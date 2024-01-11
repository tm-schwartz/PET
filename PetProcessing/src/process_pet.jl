#include("PetProcessing.jl")
#using .PetProcessing
include("utils.jl")

function runSUV(inputvolume, derivatives, templates; sidecar=nothing, CT=false, CTpattern="SOFT_BRAIN", croppet=true)
    if occursin(r"fusion"i, inputvolume)
        println("Not processing fusion $inputvolume")
        exit()
    end

    if isnothing(sidecar)
        sidecar = replace(inputvolume, "nii.gz" => "json")
    end
    subject = split(inputvolume |> basename, "_") |> first
    subderivatives = joinpath(derivatives,
        subject,
        sidecar |> basename |> splitext |> first)
    if occursin("Eq_1", inputvolume)
        sidecar = replace(sidecar, "pet_Eq_1.json" => "pet.json")
    end


    suvscalefactor = getsuvbwscalefactor(sidecar)

    if !isdir(subderivatives)
        mkpath(subderivatives)
    end

    try
        if occursin("_PT.", inputvolume)
           s = "_PT"
        else
            s = "_pet"
        end
        suffix = s =>"_resampled"

        if croppet
            suffix = s =>"_cropped_PT"
            inputtoresample = robustfov(inputvolume, subderivatives; suffix)
            suffix="cropped"=>"resampled"
        else
            inputtoresample = inputvolume
        end

        resampled = resamplepixdims(inputtoresample, subderivatives, suffix)

        if CT
            ctpath = only(glob("*$CTpattern*nii.gz", joinpath(inputvolume |> dirname |> dirname, "CT")))
            pettoct = rigidregistration(ctpath, resampled, subderivatives, "resampled" => "affine")

        end
        affinereg = rigidregistration(joinpath(templates,
                "MNI152_PET_1mm_coreg_smoothed.nii.gz"),
            resampled,
            subderivatives,
            "resampled" => "affine")

        strippedvol = skullstrip(affinereg, subderivatives, "affine" => "stripped")

        (registeredpet, paramobj) = elastixregistration(joinpath(templates,
                "MNI152_PET_1mm_coreg_stripped_smoothed.nii.gz"),
            strippedvol,
            subderivatives,
            "stripped" => "mni152", "/data/h_vmac/schwart/fdg_pet/PetProcessing/src/param27.txt")

        #smoothedvol = smoothvolume(registeredpet, subderivatives; σ = 2.97)
        suffix = "mni152" => "mni152_SUV"
        suvvolume = computesuvvolume(registeredpet, suvscalefactor, suffix)

        csv = getmeans(suvvolume, templates)
        print(csv)
        addinfotocsv(sidecar, csv)
        zipname = subderivatives |> basename
        cd(subderivatives)
        run(`zip -9rTm $(zipname * "-intermediatefiles.zip") . -x  \*.csv \*suv\*.nii.gz \*SUV\*.nii.gz`)

        ## run(`zip -9rTm $(joinpath(subderivatives, "$zipname-intermediatefiles.zip")) $subderivatives -x  \*.csv \*suv_pet.nii.gz `)
        
    catch
        logfile = replace(basename(inputvolume), ".nii.gz" => "_log.txt")
        exc = current_exceptions()
        open(joinpath(subderivatives, logfile), "w") do f
            showerror(f, exc)
            println(f, ARGS)
        end
        println("Failed on $inputvolume")
    end
end

if length(ARGS) == 1 && only(ARGS) ∈ ("-h", "--help")
    msg = """Usage:
    --initpy: Initialize Python environment and install dependencies.
    --suv <inputvolume> <derivatives directory> <template directory> : Compute SUV volume and calculate mean SUV for ROI.
    --json <sidecar path> <keysfile> <zippath> : Add header keys from DICOM in Zip archive to JSON sidecar."""
    println(msg)
    exit()
elseif first(ARGS) == "--initpy"
    initializepythonenv()
    exit()
elseif first(ARGS) == "--suv"
    if length(ARGS) == 5 && ARGS[2] == "-S"
        inputvolume, derivatives, templates, sidecar = abspath.(ARGS[2:end])
        runSUV(inputvolume, derivatives, templates; sidecar=sidecar)
    elseif length(ARGS) == 5 && ARGS[2] == "--CT"
        inputvolume, derivatives, templates = abspath.(ARGS[3:end])
        runSUV(inputvolume, derivatives, templatesl; CT=true, CTpattern="SOFT_BRAIN")
    else
        inputvolume, derivatives, templates = abspath.(ARGS[2:end])
        runSUV(inputvolume, derivatives, templates)
    end
elseif first(ARGS) == "--json"
    sidecar, keysfile, zippath = abspath.(ARGS[2:end])
    editjsonsidecar(sidecar, keysfile; zippath=zippath)
end