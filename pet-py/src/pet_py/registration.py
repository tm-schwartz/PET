from pathlib import Path, os

import itk

from pet_py.util import set_s_q_form_code, to_path


def elastix_registration(
    fixed_volume_path: str,
    moving_volume_path: str | Path,
    out_dir: str,
    suffix: tuple[str, str] | list[str, str],
    param_file: str | None = None,
) -> tuple[str, itk.ParameterObject]:
    moving_volume_path = to_path(moving_volume_path)
    basename = moving_volume_path.name
    out_file = os.path.join(out_dir, basename.replace(*suffix))
    itk_meta_dir = Path(out_dir).joinpath(out_file.replace(".nii.gz", "-elastix_data"))
    if not itk_meta_dir.exists():
        itk_meta_dir.mkdir()
    fixed_image = itk.imread(fixed_volume_path, itk.F)
    moving_image = itk.imread(moving_volume_path, itk.F)
    parameter_object = itk.ParameterObject.New()
    pms = [
        parameter_object.GetDefaultParameterMap(pm, 3)
        for pm in ["translation", "rigid", "affine"]
    ]
    parameter_object.AddParameterMap(pms[0])
    parameter_object.AddParameterMap(pms[1])
    parameter_object.AddParameterMap(pms[2])
    parameter_object.SetParameter(0, "DefaultPixelValue", "0")
    for i in range(3):
        parameter_object.SetParameter(i, "CompressResultImage", "true")
    if param_file:
        parameter_object.AddParameterFile(param_file)
    intermediate_image, result_transform_parameters = itk.elastix_registration_method(
        fixed_image,
        moving_image,
        output_directory=str(itk_meta_dir),
        parameter_object=parameter_object,
        # log_file_name="elastix.log")
        log_to_console=False,
    )

    result_image = itk.transformix_filter(
        moving_image, transform_parameter_object=result_transform_parameters
    )
    itk.imwrite(result_image, out_file)
    set_s_q_form_code(out_file)
    return (
        out_file,
        result_transform_parameters,
    )
