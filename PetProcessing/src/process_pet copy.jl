#include("PetProcessing.jl")
#using .PetProcessing
include("utils.jl")

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

    suvscalefactor = getsuvbwscalefactor(sidecar)

    try
        m = match.(r".*PT(.?).nii.gz" , inputvolumepath)
        if isempty(suffix) && !isnothing(m) && !isempty(m.captures|>only)
            suffix = m.captures |> only
        end
        ac = match(r"ac-(\d{1,})_dt", inputvolumepath).captures |> only
        ctpath = glob("*$CTpattern*$ac*nii.gz", joinpath(inputvolumepath |> dirname |> dirname, "CT" * suffix)) |> only
        resampledpet = resamplepixdims(inputvolumepath, subderivatives, "PT$suffix.nii"=>"resampled_PT.nii")  # add `.` to get last instance of PT
        resampledct = resamplepixdims(ctpath, subderivatives, "CT$suffix.nii"=>"resampled_CT.nii")
        petaffinect = rigidregistration(resampledct, resampledpet, subderivatives, "resampled_PT" => "CT_affine_PT")
        ctbinarymask = windowimage(resampledct, 400, 1000, true)
        ctmaskstripped = skullstrip(ctbinarymask, subderivatives, "MASK"=>"MASK_stripped")
        ctwindowed = windowimage(resampledct, 40, 80)
        ctbrain = maskvolume(ctwindowed, ctmaskstripped, subderivatives)
        petstripped = skullstrip(petaffinect, subderivatives, "CT_affine_PT.nii"=>"affine_stripped_PT.nii")
        (petelastixct, _) = elastixregistration(ctbrain, petstripped, subderivatives, "affine_stripped"=>"elastix", "/Users/schwartz/projects/PET-LT/bidsdir/derivatives/param27.txt")
        ctrigidt1 = rigidregistration(joinpath(templates, "stripped-t1.nii.gz"), ctbrain, subderivatives, "_WINDOWED_mask_applied"=>"_t1_affine_ct")
        (_, paramobj) = elastixregistration(joinpath(templates, "stripped-t1.nii.gz"), ctrigidt1, subderivatives, "_affine_"=> "_elastix_", "/Users/schwartz/projects/PET-LT/bidsdir/derivatives/param27.txt")
        ctrigidt1tfm = replace(ctrigidt1, "nii.gz"=>"affine.txt")
        petrigidt1 = applytransform(ctrigidt1, petelastixct, ctrigidt1tfm, "affine", subderivatives, "_elastix"=>"_t1_rigid_pet")
        petelastixt1 = transformix(petrigidt1, paramobj, subderivatives, "_rigid_"=>"_elastix_")
        suvvolume = computesuvvolume(petelastixt1, suvscalefactor, "t1_elastix_pet"=> "SUV")
        csv = getmeans(suvvolume, templates)
        addinfotocsv(sidecar, csv)

        zipname = subderivatives |> basename
        cd(subderivatives)
        run(`zip -9rTm $(zipname * "-intermediatefiles.zip") \* -x  \*.csv \*SUV\*.nii.gz `)
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