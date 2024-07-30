import json
import re
from itertools import compress
from pathlib import Path, os
from warnings import warn

import nibabel as nib
import pandas as pd
from fsl.data.atlases import LabelAtlas, getAtlasDescription, hasAtlas, rescanAtlases
from fsl.data.image import Image
from fsl.utils.image.resample import resampleToPixdims, resampleToReference
from nilearn.image import load_img, math_img
from pydicom import Sequence

MNI = 4


def to_path(p: str | Path) -> Path:
    if not isinstance(p, Path):
        return Path(p)
    return p


def check_for_keys_present(keys_to_check: list[str], path: str | Path):
    pth = Path(path)
    search_string = "|".join(keys_to_check)
    with open(pth) as p:
        r = p.read()
    return len(re.findall(search_string, r))


def get_sequence_dict(seq: Sequence):
    return {
        x.name.replace(" ", ""): (
            x.value if not isinstance(x.value, Sequence) else get_sequence_dict(x.value)
        )
        for x in seq[0].elements()
    }


def set_s_q_form_code(input_volume_path: str | Path):
    nii = nib.load(input_volume_path)
    if nii.header.get("sform_code") != MNI or nii.header.get("qform_code") != MNI:
        nii.set_sform(nii.affine, code=MNI)
        nii.set_qform(nii.affine, code=MNI)
        nii.to_filename(input_volume_path)


def mask_volume(
    input_volume_path: str | Path,
    mask_volume_path: str | Path,
    out_dir: str,
) -> str:
    input_volume_path = to_path(input_volume_path)
    basename = input_volume_path.name
    final_name = basename.replace(".nii.gz", "_mask_applied.nii.gz")
    final_path = os.path.join(out_dir, final_name)
    input_volume = load_img(input_volume_path)
    mask_vol = load_img(mask_volume_path)
    masked_volume = math_img("iv * mv", iv=input_volume, mv=mask_vol)
    masked_volume.to_filename(final_path)
    return final_path


def resample_pix_dims(
    input_volume_path: str | Path, out_dir: str, suffix: tuple
) -> str:
    input_volume_path = to_path(input_volume_path)
    input_image = Image(input_volume_path)
    resampled_data = resampleToPixdims(input_image, [1, 1, 1])
    ni = nib.Nifti1Image(*resampled_data)
    basename = input_volume_path.name.replace(*suffix)
    out_file = os.path.join(out_dir, basename)
    ni.to_filename(out_file)
    return out_file


def resample_to_ref(input_volume: Image, lname: str, reference: str | Path) -> Path:
    reference = to_path(reference)
    resampled_data = resampleToReference(input_volume, Image(reference))
    nii = nib.Nifti1Image(*resampled_data)
    out_file = reference.parent.joinpath(lname + "-mask.nii.gz")
    nii.to_filename(out_file)
    return out_file


def generate_masks(
    ref_img: str | Path,
    atlases: tuple[str] = (
        "mni",
        "harvardoxford-subcortical",
        "harvardoxford-cortical",
    ),
):
    rescanAtlases()
    if not all(_atlases := [not hasAtlas(atlas) for atlas in atlases]):
        warn("Unable to find atlas: ", *compress(atlases, _atlases), stacklevel=1)
    label_atlases = [LabelAtlas(getAtlasDescription(atl), 1.0) for atl in atlases]
    for lat in label_atlases:
        for label in lat.desc.labels:
            label_name = "".join([label.name.replace(" ", "_"), "_", lat.desc.atlasID])
            resample_to_ref(lat.get(label), label_name, ref_img)


def append_to_csv(
    sidecar_path: str | Path,
    csv_path: str | Path,
    keys_to_add: tuple[str] = (
        "PatientAge",
        "PatientSex",
        "AccessionNumber",
        "AdditionalPatientHistory",
        "AcquisitionDateTime",
        "ReasonForStudy",
    ),
    missing_val: str = " ",
):
    df = pd.read_csv(csv_path)
    with open(sidecar_path) as s:
        jsn = json.load(s)
        for key in keys_to_add:
            rfs = jsn.get(key, missing_val)
            if key in (
                "ReasonForStudy",
                "AdditionalPatientHistory",
            ):
                rfs = re.sub(r"[^\w\s.,]|\n|\t", "", rfs).strip()
            df[key] = rfs

    if "PatientAge" in keys_to_add and len(df[df.PatientAge == missing_val]):
        df["PatientAge"] = (
            (
                pd.to_datetime(jsn.get("AcquisitionDateTime"))
                - pd.to_datetime(jsn.get("PatientBirthDate"))
            )
            .to_timedelta64()
            .astype("timedelta64[Y]")
            .item()
        )
    df.to_csv(csv_path, index=False)


# def format_csv_vals(csv_path: str | Path, format_dict):
#     df = pd.read_csv(csv_path)
#     for k in format_dict.keys():
#         if callable(format_dict[k]):
#             u = df[k].apply(format_dict.pop(k))
#             df.drop(columns=[k], inplace=True)
#             df[k] = u
#     df = df.astype(format_dict)
#     df.to_csv(csv_path)
