include("utils.jl")

# robustfov -> skullstrip -> register

inputvolume, derivatives, templates = abspath.(ARGS)

if occursin(r"fusion"i, inputvolume)
    println("Not processing fusion $inputvolume")
    exit()
end

sidecar = replace(inputvolume, "nii.gz" => "json")

subject = split(inputvolume |> basename, "_") |> first

zipname = sidecar |> basename |> splitext |> first 

suvscalefactor = getsuvbwscalefactor(sidecar)

subderivatives = joinpath(derivatives, subject)

if !isdir(subderivatives)
    mkpath(subderivatives)
end

try
    croppedvol = robustfov(inputvolume, subderivatives)
    resampled = resamplepixdims(croppedvol, subderivatives, "cropped" => "resampled")
    affinereg = rigidregistration(joinpath(templates, "MNI152_PET_1mm_coreg_smoothed.nii.gz"), resampled, subderivatives, "resampled" => "affine")

    strippedvol = skullstrip(affinereg, subderivatives, "affine" => "stripped")

    registeredpet = elastixregistration(joinpath(templates, "MNI152_PET_1mm_coreg_stripped_smoothed.nii.gz"), strippedvol, subderivatives, "stripped" => "mni152")

    smoothedvol = smoothvolume(registeredpet, subderivatives; σ=2.97)


    suvvolume = computesuvvolume(smoothedvol, suvscalefactor, "pet" => "suv_pet")

    getmeans(suvvolume, templates)

    #rm(croppedvol)
    #rm(intermediatepet)
    #rm(strippedvol)
    #rm(registeredpet)
    #rm(smoothedvol)
    run(`zip -9rTm $(joinpath(subderivatives, "$zipname-intermediatefiles.zip")) $subderivatives -x \*.csv \*suv_pet.nii.gz `)
    #println(suvvolume)

catch
    logfile = replace(basename(inputvolume), ".nii.gz" => "_log.txt")
    exc = current_exceptions()
    open(joinpath(subderivatives, logfile), "w") do f
        showerror(f, exc)
    end
    println("Failed on $inputvolume")
end