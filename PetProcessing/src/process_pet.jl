include("utils.jl")

# robustfov -> skullstrip -> register

inputvolume, derivatives, templates = abspath.(ARGS)

if occursin(r"fusion"i, inputvolume)
    println("Not processing fusion $inputvolume")
    exit()
end

sidecar = replace(inputvolume, "nii.gz" => "json")

subject = split(inputvolume |> basename, "_") |> first

subderivatives = joinpath(derivatives, subject)

if !isdir(subderivatives)
    mkpath(subderivatives)
end

try
    croppedvol = robustfov(inputvolume, subderivatives)
    intermediatepet = register(joinpath(templates, "MNI152_PET_1mm.nii"), croppedvol, subderivatives, "cropped" => "intermediatereg")
    # cropped_pet.nii.gz is input to skullstrip
    strippedvol = skullstrip(intermediatepet, subderivatives, "intermediatereg" => "stripped")

    registeredpet = register(joinpath(templates, "MNI152_T1_1mm_Brain.nii.gz"), strippedvol, subderivatives, "stripped" => "mni152")

    smoothedvol = smoothvolume(registeredpet, subderivatives)

    suvscalefactor = getsuvbwscalefactor(sidecar)

    suvvolume = computesuvvolume(smoothedvol, suvscalefactor, "pet" => "suv_pet")
    rm(croppedvol)
    rm(intermediatepet)
    rm(strippedvol)
    rm(registeredpet)
    rm(smoothedvol)
    println(suvvolume)

catch e
    logfile = replace(basename(inputvolume), ".nii.gz"=>"_log.txt")
    open(joinpath(subderivatives, logfile), "w") do f
        write(f, string(e))
    end
    println("Failed on $inputvolume")
end