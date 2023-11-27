#include("PetProcessing.jl")
#using .PetProcessing
include("utils.jl")

function runSUV(inputvolume, derivatives, templates; sidecar=nothing)
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

    zipname = subderivatives |> basename 

    suvscalefactor = getsuvbwscalefactor(sidecar)

    if !isdir(subderivatives)
        mkpath(subderivatives)
    end

    try
        croppedvol = robustfov(inputvolume, subderivatives)
        resampled = resamplepixdims(croppedvol, subderivatives, "cropped" => "resampled")
        affinereg = rigidregistration(joinpath(templates,
                "MNI152_PET_1mm_coreg_smoothed.nii.gz"),
            resampled,
            subderivatives,
            "resampled" => "affine")

        strippedvol = skullstrip(affinereg, subderivatives, "affine" => "stripped")

        registeredpet = elastixregistration(joinpath(templates,
                "MNI152_PET_1mm_coreg_stripped_smoothed.nii.gz"),
            strippedvol,
            subderivatives,
            "stripped" => "mni152")

        #smoothedvol = smoothvolume(registeredpet, subderivatives; σ = 2.97)
        suvvolume = computesuvvolume(registeredpet, suvscalefactor, "pet" => "suv_pet")

        getmeans(suvvolume, templates)

        run(`zip -9rTm $(joinpath(subderivatives, "$zipname-intermediatefiles.zip")) $subderivatives -x  \*.csv \*suv_pet.nii.gz `)

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
    if length(ARGS) == 5
        inputvolume, derivatives, templates, sidecar = abspath.(ARGS[2:end])
        runSUV(inputvolume, derivatives, templates; sidecar=sidecar)
    else
        inputvolume, derivatives, templates = abspath.(ARGS[2:end])
        runSUV(inputvolume, derivatives, templates)
    end
elseif first(ARGS) == "--json"
    sidecar, keysfile, zippath = abspath.(ARGS[2:end])
    editjsonsidecar(sidecar, keysfile, zippath)
end
