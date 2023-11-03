include("utils.jl")

# robustfov -> skullstrip -> register

inputvolume, derivatives = abspath.(ARGS)

if occursin(r"fusion"i, inputvolume)
    println("Not processing fusion $inputvolume")
    exit()
end


mnitemplates = joinpath(dirname(derivatives), "mnitemplates")

sidecar = replace(inputvolume, "nii.gz" => "json")

subject = split(inputvolume |> basename, "_") |> first

subderivatives = joinpath(derivatives, subject)

if !isdir(subderivatives)
    mkpath(subderivatives)
end

try
    croppedvol = robustfov(inputvolume, subderivatives)
    intermediatepet = register(joinpath(mnitemplates, "MNI152_PET_1MM.nii"), croppedvol, subderivatives, "cropped" => "intermediatereg")
    # cropped_pet.nii.gz is input to skullstrip
    strippedvol = skullstrip(intermediatepet, subderivatives, "intermediatereg" => "stripped")

    registeredpet = register(joinpath(mnitemplates, "MNI152_T1_1mm_Brain.nii.gz"), strippedvol, subderivatives, "stripped" => "mni152")

    smoothedvol = smoothvolume(registeredpet, subderivatives)

    suvscalefactor = getsuvbwscalefactor(sidecar)

    suvvolume = computesuvvolume(smoothedvol, suvscalefactor, "pet" => "suv_pet")

    println(suvvolume)

catch e
    logfile = replace(basename(inputvolume), ".nii.gz"=>"_log.txt")
    open(joinpath(subderivatives, logfile), "w") do f
        write(f, string(e))
    end
    println("Failed on $inputvolume")
end