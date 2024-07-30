import json
import pickle
from importlib import resources as impresources
from pathlib import Path, os
from shutil import copyfile

import pandas as pd

import pet_py

input_image_path = next(
    impresources.files(pet_py).parents[1].glob("resources/*input_img.nii.gz")
)
input_json_path = next(
    impresources.files(pet_py).parents[1].glob("resources/*input_img.json")
)
resampled_image_path = next(
    impresources.files(pet_py).parents[1].glob("resources/*pixdim*.nii.gz")
)
synthstrip_image_path = next(
    impresources.files(pet_py).parents[1].glob("resources/*stripped-synthstrip*.nii.gz")
)
registered_image_path = next(
    impresources.files(pet_py).parents[1].glob("resources/*elastix*.nii.gz")
)
pickled_data_obj = (
    impresources.files(pet_py)
    .parents[1]
    .joinpath("resources", "standardized_data.pickle")
)
csv_path = next(impresources.files(pet_py).parents[1].glob("resources/*SUV.csv"))
templates_path = Path().home().joinpath("DATASET-20230414/mnitemplates")
fixed_volume_path = next(
    impresources.files(pet_py)
    .parents[1]
    .glob("resources/MNI152_PET_1mm_coreg_stripped_smoothed.nii.gz")
)
param_file = str(
    next(impresources.files(pet_py).parents[1].glob("resources/param27.txt"))
)


def test_function():
    sidecar = json.loads(input_json_path.read_text())
    sd = "_".join(
        [
            *sidecar["SeriesDescription"].split(),
            *sidecar["ImageType"],
            sidecar["StationName"],
        ]
    )
    assert (
        str(pet_py.data.generate_json_zip_path(input_json_path))
        == os.path.join(sidecar["PatientID"], sidecar["StudyInstanceUID"], sd)
    )  # str(Path("10342574").joinpath("1.2.124.113532.80.22161.20176.20160715.85837.521385965", "PET_AC_3DWB_ORIGINAL_PRIMARY_vumcct"))


class TestMods:
    def test_rfov(self, tmp_path):
        out_str = pet_py.data.robustfov(
            input_image_path, tmp_path, ("pet_input_img", "croppet_input_img.nii.gz")
        )
        assert Path(out_str).is_file()

    def test_resample_pix_dim(self, tmp_path):
        input_img = next(input_image_path.parent.glob("*croppet_input_img.nii.gz"))
        out_str = pet_py.util.resample_pix_dims(
            input_img, tmp_path, ("croppet", "pixdim")
        )
        assert Path(out_str).is_file()

    def test_skull_strip_deepbet(self, tmp_path):
        input_img = str(next(input_image_path.parent.glob("*pixdim_input_img.nii.gz")))
        out_str = pet_py.data.skull_strip_deepbet(
            input_img, tmp_path, ("pixdim", "stripped-deepbet"), 0.9, -2
        )
        assert Path(out_str).is_file()

    def test_skull_strip_synthstrip(self, tmp_path):
        input_img = str(next(input_image_path.parent.glob("*pixdim_input_img.nii.gz")))
        out_str = pet_py.data.skull_strip_mri_synthstrip(
            input_img, tmp_path, ("pixdim", "stripped-synthstrip")
        )
        assert Path(out_str).is_file()

    def test_registration(self, tmp_path):
        input_img = str(next(input_image_path.parent.glob("*stripped-synthstrip*")))
        reg = pet_py.registration.elastix_registration(
            fixed_volume_path,
            input_img,
            tmp_path,
            ("stripped-synthstrip", "elastix"),
            param_file,
        )
        assert Path(reg[0]).is_file()


def test_compute_suv_volumes():
    suv_scale_factor = pet_py.analysis.get_suvbw_scale_factor(input_json_path)
    out_data = pet_py.analysis.compute_suv_volumes(
        registered_image_path,
        suv_scale_factor,
        templates_path,
        next(templates_path.glob("*ponsvermis_1mm_bin*.nii.gz")),
    )

    with open(pickled_data_obj, "wb") as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(out_data, f, pickle.HIGHEST_PROTOCOL)

    assert Path(out_data.suv_path).is_file()
    assert Path(out_data.suvr_path).is_file()
    assert out_data.reference_region_mean > 0.0


def test_get_roi_suvr_means():
    sd = pickle.loads(pickled_data_obj.read_bytes())
    out_file = pet_py.analysis.get_roi_suvr_means(sd)
    assert Path(out_file).is_file()


def test_append_to_csv():
    copied = csv_path.parent.joinpath("temp.csv")
    copyfile(csv_path, copied)
    pet_py.util.append_to_csv(input_json_path, copied)
    r = pd.read_csv(copied)
    assert "PatientAge" in r.columns
    assert "PatientSex" in r.columns
    assert "AccessionNumber" in r.columns
    assert "AdditionalPatientHistory" in r.columns
    assert "AcquisitionDateTime" in r.columns
    assert "ReasonForStudy" in r.columns
    assert len(
        r.loc[
            :,
            [
                "PatientAge",
                "PatientSex",
                "AccessionNumber",
                "AdditionalPatientHistory",
                "AcquisitionDateTime",
                "ReasonForStudy",
            ],
        ].dropna()
    )


# def test_resample_to_ref(tmp_path):
#         input_img = Image(Path().home().joinpath("DATASET-20230414/mnitemplates/original_templates/Composite.nii.gz"))
#         out_str = pet_py.util.resample_to_ref(input_img, "Composite", fixed_volume_path)
#         assert Path(out_str).is_file()
