#include("PetProcessing.jl")
#using .PetProcessing
include("utils.jl")
using PyCall
function setup(inputvolumepath, derivatives; sidecar=nothing)
    if occursin(r"fusion"i, inputvolumepath)
        println("Not processing fusion $inputvolumepath")
        exit()
    end

    if isnothing(sidecar)
        sidecar = replace(inputvolumepath, "nii.gz" => "json")
    end
    subject = split(inputvolumepath |> basename, "_") |> first
    subderivatives = joinpath(derivatives,
        subject,
        sidecar |> basename |> splitext |> first)
    if occursin("Eq_1", inputvolumepath)
        sidecar = replace(sidecar, "pet_Eq_1.json" => "pet.json")
    end
    if !isdir(subderivatives)
        mkpath(subderivatives)
    end
    return (subderivatives, sidecar)
end

function runSUV(inputvolumepath, subderivatives, templates, sidecar, CTpattern="SOFT_BRAIN"; suffix="")
    itk = pyimport("itk")
    suvscalefactor = getsuvbwscalefactor(sidecar)
    originalt1template = joinpath(templates, "original_templates/mni_icbm152_t1_tal_nlin_asym_55_ext.nii.gz")
    strippedt1template = joinpath(templates, "stripped-t1.nii.gz")
    parameter_object = itk.ParameterObject.New()
    parameter_object.AddParameterFile("/Users/schwartz/projects/PET-LT/bidsdir/derivatives/param27 copy.txt")
    try
        m = match.(r".*PT(.?).nii.gz", inputvolumepath)
        if isempty(suffix) && !isnothing(m) && !isempty(m.captures |> only)
            suffix = m.captures |> only
        end
        ac = match(r"ac-(\d{1,})_dt", inputvolumepath).captures |> only
        ctpath = glob("*$CTpattern*$ac*nii.gz", joinpath(inputvolumepath |> dirname |> dirname, "CT" * suffix)) |> only

        RESAMPLEDCT = resamplepixdims(ctpath, subderivatives, "_CT" => "_RESAMPLED_CT")
        RESAMPLEDPT = resamplepixdims(inputvolumepath, subderivatives, "_PT" => "_RESAMPLED_PT")
        (AFFINEPETCT, AFFINEPETCTTRM) = elastixregistration(RESAMPLEDCT, RESAMPLEDPT, subderivatives, "RESAMPLED" => "affine_composed_PT", "/Users/schwartz/projects/PET-LT/bidsdir/derivatives/elastix-parameter-files/affinetransform.txt")
        # NEW
        pet_elastix_ct = joinpath(subderivatives, "PET_ELASTIX_CT")
        mkpath(pet_elastix_ct)
        fixed_image = itk.imread(RESAMPLEDCT, itk.F)
        moving_image= itk.imread(RESAMPLEDPT, itk.F)
        (ELASTIXPETCT, ELASTIXPETCTPARAMOBJ) = itk.elastix_registration_method(fixed_image, moving_image, initial_transform_parameter_object=AFFINEPETCTTRM, parameter_object=parameter_object, output_directory=pet_elastix_ct, log_to_console=true)

        #=
            possibly run affine pet->ct then transformix before transformixpet->ctaffine->ctmrelastix
        =#
        
        (AFFINECTMR, AFFINECTMRTRM) = elastixregistration(originalt1template, RESAMPLEDCT, subderivatives, "RESAMPLED" => "affine_composed_CT", "/Users/schwartz/projects/PET-LT/bidsdir/derivatives/elastix-parameter-files/affinetransform.txt")
        CTBINARYMASK = windowimage(RESAMPLEDCT, 400, 1000, true)
        CTMASKSTRIPPED = skullstrip(CTBINARYMASK, subderivatives, "MASK" => "MASK_stripped")
        CTWINDOWED = windowimage(RESAMPLEDCT, 40, 80)
        CTBRAIN = maskvolume(CTWINDOWED, CTMASKSTRIPPED, subderivatives)
        fixed_image = itk.imread(strippedt1template, itk.F) # change to correct stripped template
        moving_image = itk.imread(CTBRAIN, itk.F)
        ct_elastix_mr = joinpath(subderivatives, "CT_ELASTIX_MR")
        mkpath(ct_elastix_mr)
        (ELASTIXCTMR, ELASTIXCTMRPARAMOBJ) = itk.elastix_registration_method(fixed_image, moving_image, initial_transform_parameter_object=AFFINECTMRTRM, parameter_object=parameter_object, output_directory=ct_elastix_mr, log_to_console=true)
        outelastixctmr = joinpath(subderivatives, "CT_ELASTIX_MR_OUT.nii.gz")
        itk.imwrite(ELASTIXCTMR, outelastixctmr)
        STRIPPEDPET = skullstrip(RESAMPLEDPT, subderivatives, "RESAMPLED" => "STRIPPED")
        moving_image = itk.imread(STRIPPEDPET, itk.F)
        pet_elastix_transformix_affine_ct = joinpath(subderivatives, "PET_ELASTIX_transformix_affine_CT")
        mkpath(pet_elastix_transformix_affine_ct)
        tfmx_pet_affine_ct = itk.transformix_filter(moving_image, transform_parameter_object=ELASTIXPETCTPARAMOBJ, output_directory=pet_elastix_transformix_affine_ct)
        outtransformixpet = joinpath(subderivatives, "PET_TRANSFORMIX_ct_affine_MR_OUT.nii.gz")
        itk.imwrite(tfmx_pet_affine_ct, outtransformixpet)
        trnsformix_pet = joinpath(subderivatives, "trnsformixpet")
        mkpath(trnsformix_pet)
        TRANSFORMIXPET = itk.transformix_filter(tfmx_pet_affine_ct, transform_parameter_object=ELASTIXCTMRPARAMOBJ, output_directory=trnsformix_pet, log_to_console=true)
        outtransformixpet = joinpath(subderivatives, "PET_TRANSFORMIX_MR_OUT.nii.gz")
        itk.imwrite(TRANSFORMIXPET, outtransformixpet)
        suvvolume = computesuvvolume(outtransformixpet, suvscalefactor, "TRANSFORMIX_MR_OUT" => "SUV")
        csv = getmeans(suvvolume, templates)
        addinfotocsv(sidecar, csv)

        zipname = subderivatives |> basename
        # run(`zip -9rTm $(joinpath(subderivatives, "$zipname-intermediatefiles.zip")) $subderivatives -x  \*.csv \*SUV\*.nii.gz `)
    catch
        logfile = replace(basename(inputvolumepath), ".nii.gz" => "_log.txt")
        exc = current_exceptions()
        open(joinpath(subderivatives, logfile), "w") do f
            showerror(f, exc)
            println(f, ARGS)
        end
        println("Failed on $inputvolumepath")
    end
end


inputvolumepath, derivatives, templates = ARGS

(outdir, sidecar) = setup(inputvolumepath, derivatives)

runSUV(inputvolumepath, outdir, templates, sidecar)