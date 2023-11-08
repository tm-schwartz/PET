using PackageCompiler
using ZipArchives
using Query
using JSON3
using Chain
using Glob
using Dates
using PyCall
using Conda
using NIfTI

const KGTOGRAMS = 1000

function initializepythonenv()
    Conda.pip_interop(true)
    Conda.pip("install", "simpleitk")
    Conda.pip("install", "deepbet")
    Conda.pip("install", "dcm2niix")
end


function generatejsonzippath(path)
    # mrn/studyuid/seriesdescription imagetype stationname/
    pathbasename = basename(path)
    mrn = @chain pathbasename begin
        split("_")
        first
        split("-")
        last
    end
    studyinstanceuid = @chain path JSON3.read get(_, "StudyInstanceUID")
    dicomfolder = match(r"(\d{0,}[A-Z].+)_ac", pathbasename).captures |> only
    return joinpath(mrn, studyinstanceuid, dicomfolder)
end

function checkforkeyspresent(keystocheck, path)
    searchstring = join(keystocheck, '|')
    return parse(Int, readchomp(`grep -Ec $searchstring $path`))
end

function editjsonsidecar(skipexistingkeys=false)
    sidecar = readlines(ENV["SIDECAR_LIST"])[parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
    keystoadd = readlines("derivatives/additionalheaderkeys.txt")
    nkeys = length(keystoadd)

    if skipexistingkeys && checkforkeyspresent(keystoadd, sidecar) == nkeys
        return "Skipping $sidecar, keys exist"
    end

    jsonzippath = generatejsonzippath(sidecar)
    searchregex = (Regex ∘ joinpath)(jsonzippath, ".*json")
    mrn = (first ∘ splitpath)(jsonzippath)
    ziparchive = (zip_open_filereader ∘ joinpath)(ENV["FDG_ZIP_PATH"], mrn * ".zip")
    keystoadd = readlines("derivatives/additionalheaderkeys.txt")

    headervalues = @chain ziparchive begin
        zip_names
        _[occursin.(searchregex, _)]
        sort
        first
        zip_readentry(ziparchive, _, String)
        JSON3.read
        get.(Ref(_), keystoadd, missing)
    end

    sidecardict = @chain sidecar read(_, String) JSON3.read copy

    open(sidecar, "w") do io
        JSON3.pretty(io, merge(sidecardict, Dict(Symbol(k) => coalesce(v, "missing") for (k, v) in zip(keystoadd, headervalues))))
    end

    indxmissing = BitArray(ismissing.(headervalues))

    if (m = sum(indxmissing)) > 0
        return "Added $(nkeys - m)/$nkeys to $sidecar. Missing $(keystoadd[indxmissing])"
    end

    return "Added all keys to sidecar:$sidecar"
end

function smoothvolume(inputvolume, outdir)
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, "pet" => "8mm_pet"))
    badgerfsl = `/data/h_vmac/vmac_imaging/fsl_v6.0.5.1.sif`
    cmd = `fslmaths $inputvolume -s 3.4 $outfile`
    if isfile(replace(string(badgerfsl), "`"=>""))
        run(`singularity exec $badgerfsl $cmd`)
    else
        run(cmd)
    end
    return outfile
end

function getsuvbwscalefactor(sidecarpath)
    json = (JSON3.read ∘ read)(sidecarpath)
    suvbwscalefactor = let
        seriestime = Time(json.SeriesTime, dateformat"HHMMSS")
        radiopharmaceuticalstarttime = Time(split(json.RadiopharmaceuticalStartTime, '.') |> first, dateformat"HHMMSS")
        halflife = Second(floor(json.RadionuclideHalfLife))
        injecteddosebq = json.RadionuclideTotalDose
        weightkg = json.PatientWeight
        decaytime = Second(seriestime - radiopharmaceuticalstarttime)
        decaydose = injecteddosebq * 2^(-decaytime / halflife)
        weightkg * KGTOGRAMS / decaydose
    end
    return suvbwscalefactor
end

function computesuvvolume(inputvolume, suvscalefactor, suffix)
    pet = niread(inputvolume)
    outpath = replace(inputvolume, suffix)
    niwrite(outpath, NIVolume(pet.header, pet.extensions, pet * suvscalefactor))
    return outpath
end

function setsformqform(inputvolume)
    nii = niread(inputvolume)
    if nii.header.sform_code != 4 || nii.header.qform_code != 4
        nii.header.sform_code = 4
        nii.header.qform_code = 4
        niwrite(inputvolume, nii)
    end
end

function loadsidecar(sidecar)
    return JSON3.read(sidecar)
end

function robustfov(inputvolume, outdir)
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, "pet" => "cropped_pet"))
    if occursin(r"HN"i, bn)
        cp(inputvolume, outfile)
    else
        badgerfsl = `/data/h_vmac/vmac_imaging/fsl_v6.0.5.1.sif`
        cmd = `robustfov -i $inputvolume -r $outfile`
        if isfile(replace(string(badgerfsl), "`"=>""))
            run(`singularity exec $badgerfsl $cmd`)
        else
            run(cmd)
        end
    end
    return outfile
end

function skullstrip(inputvolume, outdir, suffix)
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, suffix))
    tivpath = joinpath(outdir, replace(outfile, "nii.gz" => "TIV.csv"))
    deepbet = pyimport("deepbet")
    deepbet.run_bet(PyObject([inputvolume]), PyObject([outfile]), tiv_paths=PyObject([tivpath]), threshold=0.9, n_dilate=-2, no_gpu=true)
    return outfile
end

py"""
import SimpleITK as sitk


def register(
    fixed_image_path, moving_image_path, final_moving_image, final_transform_path
):
    fixed_image = sitk.ReadImage(fixed_image_path, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(moving_image_path, sitk.sitkFloat32)

    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image,
        moving_image,
        sitk.Euler3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY,
    )

    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    registration_method.SetInterpolator(sitk.sitkLinear)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=100,
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10,
    )
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Setup for the multi-resolution framework.
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    registration_method.SetInitialTransform(initial_transform)

    final_transform = registration_method.Execute(
        sitk.Cast(fixed_image, sitk.sitkFloat32),
        sitk.Cast(moving_image, sitk.sitkFloat32),
    )
    if final_transform_path != None:
        sitk.WriteTransform(final_transform, final_transform_path)
    resampled_moving = sitk.Resample(
        moving_image,
        fixed_image,
        final_transform,
        sitk.sitkLinear,
        0.0,
        moving_image.GetPixelID(),
    )

    sitk.WriteImage(resampled_moving, final_moving_image)


def resample_moving(
    fixed_image_path, moving_image_path, final_moving_image, resample_with_transform
):
    fixed_image = sitk.ReadImage(fixed_image_path, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(moving_image_path, sitk.sitkFloat32)
    transform = sitk.ReadTransform(resample_with_transform)
    # apply transform
    resampled_moving = sitk.Resample(
        moving_image,
        fixed_image,
        transform,
        sitk.sitkLinear,
        0.0,
        moving_image.GetPixelID(),
    )

    sitk.WriteImage(resampled_moving, final_moving_image)
"""

function register(fixedvolumepath, movingvolumepath, outdir, suffix, tfm=false)
    bn = basename(movingvolumepath)
    finalmovingvolumepath = joinpath(outdir, replace(bn, suffix))
    if tfm
        finaltransformpath = joinpath(outdir, replace(finalmovingvolumepath, "nii.gz" => "tfm"))
        py"register"(fixedvolumepath, movingvolumepath, finalmovingvolumepath, finaltransformpath)
    else
        py"register"(fixedvolumepath, movingvolumepath, finalmovingvolumepath, PyObject(nothing))
    end
    setsformqform(finalmovingvolumepath)
    return finalmovingvolumepath
end

function resample(fixedvolumepath, movingvolumepath, outdir, suffix)
    bn = basename(movingvolumepath)
    finalmovingvolumepath = joinpath(outdir, replace(bn, suffix))
    resamplewithtransform = joinpath(outdir, replace(bn, "nii.gz" => "tfm"))
    py"resample_moving"(fixedvolumepath, movingvolumepath, finalmovingvolumepath, resamplewithtransform)
    setsformqform(finalmovingvolumepath)
end


function compilepackage()
    ENV["JULIA_CPU_TARGET"] = "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"
    cmd = `julia --project=PetProcessing --trace-compile=pettrace.jl -e "push!(LOAD_PATH, \"PetProcessing\"); using PetProcessing; exit()"`
    run(cmd)
    create_sysimage(["PetProcessing"]; sysimage_path="PetProcessing-img.so", precompile_statements_file="pettrace.jl")
end
