using PackageCompiler
using ZipArchives
using Query
using JSON3
using Chain
using Glob
using Dates

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
    sidecar = readlines(ENV["SIDECAR_LIST"])[parse(Int,ENV["SLURM_ARRAY_TASK_ID"])]
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

function getsuvbwscalefactor(path)
    json = JSON3.read(path)
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


function compilepackage()
    ENV["JULIA_CPU_TARGET"] = "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"
    cmd = `julia --project=PetProcessing --trace-compile=pettrace.jl -e "push!(LOAD_PATH, \"PetProcessing\"); using PetProcessing; exit()"`
    run(cmd)
    create_sysimage(["PetProcessing"]; sysimage_path="PetProcessing-img.so", precompile_statements_file="pettrace.jl")
end
