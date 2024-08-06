import pet_py
from zipfile import ZipFile
from pathlib import Path
from shutil import rmtree
import re
import argparse
from more_itertools import unique, strip
from warnings import warn
from multiprocessing import Pool
from pydicom.misc import is_dicom
from pandas import Series, concat, read_csv

parser = argparse.ArgumentParser()
parser.add_argument("--template_dir", required=False)
parser.add_argument("--param_file", required=False)
parser.add_argument("--zip_file", required=False)
parser.add_argument("--nifti_path", required=False)
parser.add_argument("data_path")


def check_dcm(path: Path):
    if path.is_file() and not is_dicom(path):
        path.unlink()
        return str(path)
    return ""


def process_zipfile(zip_path: str, data_path: Path):
    # unzip file to zipdest
    # check each file if dicom
    # follow notebook from-formatted-dicoms
    zip_dir = data_path.joinpath(
        pet_py.util.to_path(zip_path).name.removesuffix(".zip")
    )
    zip_parent = pet_py.util.to_path(zip_path).parent
    header_keys = data_path.parent.joinpath("additionalheaderkeys.txt")
    bids_outdir = data_path.joinpath("bids")
    derivatives = bids_outdir.joinpath("derivatives")
    derivatives.mkdir(parents=True, exist_ok=True)
    data_csv = derivatives.joinpath("datapaths.csv")

    with ZipFile(zip_path) as myzip:
        myzip.extractall(str(data_path))

    with Pool(8) as p:
        warn(
            f"Removing extrenuous files: {list(strip(unique(p.map(check_dcm, zip_dir.rglob('*[0-9]'))), lambda x: x == ''))}"
        )

    for study_path in zip_dir.glob("[0-9]*"):
        if not study_path.is_dir():
            continue
        for series_path in filter(
            lambda x: x.is_dir()
            and re.search(
                r"legs|mip|roi",
                x.name,
                flags=re.IGNORECASE,
            )
            is None,
            study_path.iterdir(),
        ):  # actual 2D images are ignored in dcm2niix config
            try:
                pet_py.data.dcm2niix(series_path)
            except StopIteration:
                continue

    for nii in zip_dir.rglob("*.nii.gz"):
        pet_py.data.to_bids(nii, bids_outdir)

    nii = Series(
        bids_outdir.joinpath(f"sub-{zip_dir.name}").rglob("*.nii.gz"),
        name="nifti_paths",
    )
    nii = nii[~nii.astype(str).str.contains(r"derivatives")]
    nii.astype(str).str.replace("nii.gz", "json").apply(
        pet_py.data.edit_json_sidecar,
        args=(
            header_keys,
            zip_parent,
        ),
    )
    nii = nii.apply(lambda x: x.relative_to(data_path)).astype(str)
    sidecars = nii.str.replace("nii.gz", "json").rename("sidecar_paths")
    df = concat([nii, sidecars], axis=1)
    if data_csv.is_file():
        df = concat([read_csv(data_csv), df], ignore_index=True).drop_duplicates(
            subset="nifti_paths", keep="last"
        )
    df.to_csv(data_csv, index=False)

    rmtree(str(zip_dir), ignore_errors=True)


def runner(
    nifti_path: Path,
    bids_dir: Path,
    templates: Path,
    param_file: Path,
):
    # todo change from using zip parent in edit_json_data
    derivatives = bids_dir.joinpath(
        "derivatives",
        re.search(r"(sub-\d{7,})_", nifti_path.name).group(1),
        re.search(r"(?<=sdit-)(.*)\.nii.gz", nifti_path.name).group(1),
    )
    derivatives.mkdir(parents=True, exist_ok=True)
    sidecar = str(nifti_path).replace("nii.gz", "json")

    robfov_out = pet_py.data.robustfov(nifti_path, derivatives, 300, ("PT", "cropped_pet"))
    pixdim = pet_py.util.resample_pix_dims(
        robfov_out, derivatives, ("cropped_pet", "pixdim")
    )
    reg_to_mni, _ = pet_py.registration.elastix_registration(
        templates.joinpath("MNI152_PET_1mm_coreg_smoothed.nii.gz"),
        pixdim,
        derivatives,
        ("pixdim", "elastix"),
        str(param_file),
    )
    stripped_deepbet = pet_py.data.skull_strip_deepbet(
        reg_to_mni, str(derivatives), ("elastix", "stripped_deepbet")
    )
    # stripped_deepbet = next(derivatives.rglob("*stripped_deepbet.nii.gz"))
    suv_scaler = pet_py.analysis.get_suvbw_scale_factor(sidecar)
    suv_data = pet_py.analysis.compute_suv_volumes(
        stripped_deepbet,
        suv_scaler,
        templates,
        next(templates.glob("ponsvermis_1mm*mask.nii.gz")),
    )
    csv_out = pet_py.analysis.get_roi_suvr_means(suv_data)
    pet_py.util.append_to_csv(sidecar, csv_out)

    # for obj in derivatives.glob("*"):
    #     if re.search(r"SUV|SUVR", obj.name):
    #         continue
    #     elif obj.is_dir():
    #         rmtree(str(obj), ignore_errors=True)
    #     else:
    #         obj.unlink(missing_ok=True)

if __name__ == "__main__":
    args = parser.parse_args()
    data_path = Path(args.data_path)
    if args.nifti_path is None:
        process_zipfile(args.zip_file, data_path)
    else:
        bids_dir = data_path.joinpath("bids-example2")
        param_file = pet_py.util.to_path(args.param_file)
        template_dir = pet_py.util.to_path(args.template_dir)
        nifti_path = pet_py.util.to_path(args.nifti_path)
        runner(nifti_path, bids_dir, template_dir, param_file)
