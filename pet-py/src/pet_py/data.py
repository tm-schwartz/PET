import json
import re
import shutil
import subprocess
import zipfile
from functools import partial
from itertools import starmap
from pathlib import Path, os
from warnings import warn

import pydicom as pydcm
from deepbet import run_bet
from more_itertools import consume
from nipype.interfaces.dcm2nii import Dcm2niix
from nipype.interfaces.fsl.utils import RobustFOV

from pet_py.util import check_for_keys_present, get_sequence_dict, to_path


def copy_data(source: str | Path, dest: Path):
    if not dest.parents[0].is_dir():
        dest.parents[0].mkdir(parents=True)
    shutil.copyfile(source, dest)


def replc(m: re.Match) -> str:
    if m.group().isspace():
        return "_"
    return ""


def format_source_data(data_dir: str | Path) -> Path:
    data_dir = to_path(data_dir)
    file_map = []
    out_dir = data_dir.joinpath("stagingdir")
    for path in data_dir.rglob("*"):
        if (
            path.is_file()
            and not re.match(".*(DICOMDIR|Viewer)", str(path))
            and pydcm.misc.is_dicom(path)
        ):
            with pydcm.dcmread(path, defer_size="1 KB") as tag_data:
                high_level_siuid = tag_data.get("StudyInstanceUID")
                if isinstance(
                    seq_or_uid := tag_data.get(
                        "RequestAttributesSequence", high_level_siuid
                    ),
                    pydcm.Sequence,
                ):
                    study_instance_uid = seq_or_uid[0].get(
                        "StudyInstanceUID", high_level_siuid
                    )
                else:
                    study_instance_uid = high_level_siuid
                file_name = ".".join(
                    [
                        tag_data.get("Modality", "NONE"),
                        tag_data.get("SOPInstanceUID", "NONE"),
                        str(tag_data.get("InstanceNumber", path.name)),
                        str(tag_data.get("ImageIndex", 1)),
                        "dcm",
                    ]
                )
                image_type = "_".join(tag_data.get("ImageType", "NONE")).replace(
                    " ", "_"
                )

                series_data = re.sub(
                    r"[^\w]",
                    replc,
                    "_".join(
                        [
                            tag_data.get("SeriesDescription", " "),
                            image_type,
                            tag_data.get("StationName", "NONE").replace(" ", "_"),
                            tag_data.get("Modality", "NONE"),
                        ]
                    ),
                )

                if not study_instance_uid:
                    print(path)
                source = str(path.resolve())
                dest = out_dir.joinpath(study_instance_uid, series_data, file_name)
                file_map.append((source, dest))
                if (jsn := to_path(f"{source}.json")).is_file():
                    file_map.append((jsn, dest.with_suffix(".json")))
    # ready for multithreading
    consume(starmap(copy_data, file_map))
    return out_dir


def dcm2niix(series_path: str | Path):
    series_path = to_path(series_path)
    folder = series_path.name
    with pydcm.dcmread(next(series_path.glob("*.dcm")), defer_size="1 KB") as dcm:
        modality = dcm.get("Modality", "NONE")
    if series_path.name.count(modality) > 1:
        folder = series_path.name.split("_")
        folder.reverse()
        i = folder.index(modality)
        _ = folder.pop(i)
        folder.reverse()
        folder = "_".join(folder)

    converter = Dcm2niix()
    converter.inputs.source_dir = str(series_path)
    converter.inputs.output_dir = str(series_path)
    converter.inputs.compression = 9
    converter.inputs.ignore_deriv = True
    converter.inputs.merge_imgs = True
    converter.inputs.anon_bids = False
    converter.inputs.out_filename = f"sub-%i_sdit-{folder}_ac-%g_dt-%t_{modality}"
    converter.run()


def to_bids(file_path: str | Path, out_dir: str | Path):
    file_path = to_path(file_path)
    out_dir = to_path(out_dir)
    modality = file_path.name.split(".")[0].split("_")[-1]
    subdir = out_dir.joinpath(file_path.name.split("_")[0], modality)
    subdir.mkdir(parents=True, exist_ok=True)
    dest = subdir.joinpath(file_path.name)
    shutil.move(file_path, dest, copy_function=shutil.copyfile)
    if (
        jsn := file_path.parent.joinpath(file_path.name.replace("nii.gz", "json"))
    ).is_file():
        shutil.move(jsn, subdir.joinpath(jsn.name), copy_function=shutil.copyfile)


def generate_json_zip_path(path: str | Path) -> Path:
    pth = Path(path)
    basename = pth.name
    mrn = basename.split("_")[0].split("-")[-1]
    with open(pth) as p:
        jsn = json.load(p)
    study_instance_uid = jsn.get("StudyInstanceUID", "NOSTUDYIUID")
    dicom_folder = re.match(".*(?<=sdit-)(.+)_ac.*", basename).group(1)
    return Path(mrn).joinpath(study_instance_uid, dicom_folder)


def _edit_json_sidecar_dcm(
    sidecar: str | Path,
    keys_file: Path,
    dicom: str | Path,
    skip_existing_keys: bool = False,
):
    keys_to_add = keys_file.read_text().split()
    nkeys = len(keys_to_add)
    if skip_existing_keys and check_for_keys_present(keys_to_add, sidecar) == nkeys:
        return f"Skipping {sidecar}, keys exist"
    dcm = pydcm.dcmread(dicom)

    with open(sidecar, "r+") as s:
        sidecar_dict = json.load(s)
        try:
            radio_seq = get_sequence_dict(dcm["RadiopharmaceuticalInformationSequence"])
            sidecar_dict.update(radio_seq)
        except Exception:
            warn(
                f"unable to add RadiopharmaceuticalInformationSequence to {sidecar} from {dcm.filename}"
            )
        for key in keys_to_add:
            if dcm.get(key) is None and sidecar_dict.get(key) is None:
                print(f"Missing {key} for sidecar {sidecar}")
                continue
            sidecar_dict[key] = dcm.get(key).strip()
        s.seek(0)
        s.truncate()
        json.dump(sidecar_dict, s, indent=1)
    return


def _edit_json_sidecar(
    sidecar: str | Path,
    keys_file: Path,
    zip_path: str | Path,
    skip_existing_keys: bool = False,
):
    keys_to_add = keys_file.read_text().split()
    nkeys = len(keys_to_add)
    if skip_existing_keys and check_for_keys_present(keys_to_add, sidecar) == nkeys:
        return f"Skipping {sidecar}, keys exist"
    json_zip_path = generate_json_zip_path(sidecar)
    search_regex = os.path.join(json_zip_path, ".*\d\.json")
    mrn = json_zip_path.parts[0]
    search_func = partial(re.match, search_regex)
    with zipfile.ZipFile(
        zip_path.joinpath(mrn + ".zip")
    ) as zip_archive, zip_archive.open(
        sorted(
            filter(search_func, zip_archive.namelist()),
            key=lambda m: int(m.split("/")[-1].split(".")[-2]),
        )[0]
    ) as jsn_dumped:
        loaded = json.load(jsn_dumped)
        update_sidecr = {k: loaded[k] for k in keys_to_add if k in loaded}
    if len(list(update_sidecr.keys())) != nkeys:
        for key in keys_to_add:
            if key not in update_sidecr:
                print(f"Missing {key} for sidecar {sidecar}")
    with open(sidecar, "r+") as s:
        sidecar_dict = json.load(s)
        sidecar_dict.update(update_sidecr)
        s.seek(0)
        s.truncate()
        json.dump(sidecar_dict, s, indent=1)
    return


def edit_json_sidecar(
    sidecar: str | Path,
    keys_file: str | Path,
    zip_path: str | Path | None = None,
    dicom: str | Path | None = None,
    skip_existing_keys: bool = False,
):
    sidecar = to_path(sidecar)
    if re.search(r"Eq_1", sidecar.name):
        sidecar_str = re.sub(r"_Eq_1", "", str(sidecar))
        folder = sidecar_str.removesuffix(".json").split("_")[-1]
        sidecar_og = re.sub(r"(?<=/)1(?=/)", folder, sidecar_str)
        shutil.copyfile(sidecar_og, str(sidecar))
    keys_fl = to_path(keys_file)
    if zip_path:
        zip_path = to_path(zip_path)
        _edit_json_sidecar(sidecar, keys_fl, zip_path, skip_existing_keys)
    elif dicom:
        dcm = to_path(dicom)
        _edit_json_sidecar_dcm(sidecar, keys_fl, dcm, skip_existing_keys)
    else:
        shd = "Should not be reached!"
        raise shd


def robustfov(
    input_volume: str | Path,
    out_dir: str,
    brain_size: int = 170,
    suffix: tuple[str, str] | list[str, str] = ("pet", "cropped_pet"),
) -> str:
    input_volume = to_path(input_volume)
    basename = input_volume.name
    out_file = os.path.join(out_dir, basename.replace(*suffix))
    out_transform = out_file.replace("nii.gz", "tfm.txt")
    if "hn" in basename.lower():
        shutil.copy(input_volume, out_file)
    else:
        robustfov_obj = RobustFOV()
        robustfov_obj.inputs.in_file = str(input_volume)
        robustfov_obj.inputs.out_roi = out_file
        robustfov_obj.inputs.brainsize = brain_size
        robustfov_obj.inputs.out_transform = out_transform
        robustfov_obj.inputs.output_type = "NIFTI_GZ"
        robustfov_obj.run()
    return out_file


def skull_strip_deepbet(
    input_volume_path: str | Path,
    out_dir: str,
    suffix: tuple[str, str] | list[str, str],
    threshold: float = 0.9,
    dilation: int = -2,
) -> str:
    input_volume_path = to_path(input_volume_path)
    basename = input_volume_path.name
    out_file = os.path.join(out_dir, basename.replace(*suffix))
    tiv_path = out_file.replace("nii.gz", "TIV.csv")
    run_bet(
        [str(input_volume_path)],
        [out_file],
        tiv_paths=[tiv_path],
        threshold=threshold,
        n_dilate=dilation,
        no_gpu=True,
    )
    return out_file


def skull_strip_mri_synthstrip(
    input_volume_path: str | Path,
    out_dir: str,
    suffix: tuple[str, str] | list[str, str],
    freesurfer_path: str = os.environ.get("FREESURFER_HOME"),
    singularity_path: str | None = None,
) -> str:
    input_volume_path = to_path(input_volume_path)
    basename = input_volume_path.name
    out_file = os.path.join(out_dir, basename.replace(*suffix))
    if singularity_path:
        subprocess.run(
            f"singularity exec {singularity_path} mri_synthstrip -i {input_volume_path!s} -o {out_file}",
            shell=True,
            check=False,
        )
    elif freesurfer_path:
        synthstrip = next(to_path(freesurfer_path).glob("bin/mri_synthstrip"))
        subprocess.run(
            f"{synthstrip} -i {input_volume_path!s} -o {out_file}",
            env={"FREESURFER_HOME": freesurfer_path},
            shell=True,
            check=False,
        )
    else:
        missing_path = "Need freesurfer_path or singularity_path"
        raise missing_path
    return out_file
