import json
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from nilearn.image import load_img, math_img, smooth_img

from pet_py.util import to_path


@dataclass
class StandardizedData:
    suv_path: Path
    suvr_path: Path
    reference_region_mean: np.float64
    template_dir: Path


def get_suvbw_scale_factor(sidecar: str | Path) -> float:
    with open(sidecar) as s:
        sidecar_json = json.load(s)
    series_time = datetime.strptime(sidecar_json.get("SeriesTime").split(".")[0], "%H%M%S")
    radiopharmaceutical_start_time = datetime.strptime(
        sidecar_json.get("RadiopharmaceuticalStartTime").split(".")[0], "%H%M%S"
    )
    half_life = sidecar_json.get("RadionuclideHalfLife")
    injected_dose_bq = sidecar_json.get("RadionuclideTotalDose")
    weight_kg = sidecar_json.get("PatientWeight")
    decay_time = series_time - radiopharmaceutical_start_time
    decay_dose = injected_dose_bq * 2 ** (-decay_time.seconds / half_life)
    return weight_kg * 1000 / decay_dose


def compute_suv_volumes(
    input_volume_path: str | Path,
    suv_scale_factor: float,
    template_dir: str | Path,
    reference_region_mask_file_name: str | Path,
    fwhm: int = 7.0,
    pons_vermis_reference_region: bool = True,
) -> StandardizedData:
    input_volume_path = to_path(input_volume_path)
    template_dir = to_path(template_dir)
    input_volume = load_img(input_volume_path)
    reference_region = load_img(template_dir.joinpath(reference_region_mask_file_name))
    reference_region_mask = math_img("np.where(mask == 0.0, np.nan, 1.0)", mask=reference_region)
    out_path_suv_image = input_volume_path.parents[0].joinpath(input_volume_path.name.replace(".nii.gz", "-SUV.nii.gz"))
    smoothed_volume = smooth_img(input_volume, fwhm=fwhm)
    suv_image = math_img(f"sv * {suv_scale_factor}", sv=smoothed_volume)
    suv_image.to_filename(out_path_suv_image)
    if pons_vermis_reference_region:
        suv_reference_region_array = (
            math_img("suvimg * rr_mask", suvimg=suv_image, rr_mask=reference_region_mask).get_fdata().flatten()
        )
        suv_reference_region = suv_reference_region_array.compress(~np.isnan(suv_reference_region_array))
        part_val = suv_reference_region.size // 2
        suv_reference_region.partition(part_val)
        suv_reference_region_mean = np.mean(suv_reference_region[part_val:])
    else:
        suv_reference_region_mean = np.nanmean(
            math_img("suvimg * rr_mask", suvimg=suv_image, rr_mask=reference_region_mask).get_fdata()
        )
    out_path_suvr_image = input_volume_path.parents[0].joinpath(
        input_volume_path.name.replace(".nii.gz", "-SUVR.nii.gz")
    )
    math_img(f"suvimg / {suv_reference_region_mean}", suvimg=suv_image).to_filename(out_path_suvr_image)
    return StandardizedData(out_path_suv_image, out_path_suvr_image, suv_reference_region_mean, template_dir)


def get_roi_suvr_means(data_class: StandardizedData)->str:
    """
    calculate suvr means (mean of roi in suv / reference region). composite becomes metaroi stat.
    """
    roi_masks = list(data_class.template_dir.glob("*-mask.nii.gz"))
    suv_image = load_img(data_class.suv_path)
    out_file = data_class.suv_path.parent.joinpath(data_class.suv_path.name.split(".")[0] + ".csv")
    mrn = re.search(r"(?<=sub-)\d{7,}", data_class.suv_path.name).group()
    mean_vals = np.zeros(len(roi_masks))
    columns = []
    for i, roi in enumerate(roi_masks):
        columns.append(re.sub(r".nii.gz|-mask|\,|\\", "", roi.name).replace("harvaroxford", "harvox").replace("cortical", "cort"))
        roi_mask = math_img("np.where(mask == 0.0, np.nan, 1.0)", mask=load_img(roi))

        mean_vals[i] = (
            np.nanmean(math_img("suv_image * roi_mask", suv_image=suv_image, roi_mask=roi_mask).get_fdata())
            / data_class.reference_region_mean
        )

    df = pd.DataFrame(np.expand_dims(mean_vals, axis=0), columns=columns)
    df["mrn"] = mrn
    cols = ["mrn"]
    cols.extend(df.columns[df.columns != "mrn"])
    df.loc[:, cols].to_csv(out_file, index=False)
    return out_file
