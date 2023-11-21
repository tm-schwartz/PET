using ZipArchives
using Query
using JSON3
using Chain
using Glob
using Dates
using PyCall
using Conda
using NIfTI
using Pkg
using Statistics: mean
using DataFrames
using CSV
const KGTOGRAMS = 1000

function initializepythonenv()
    Conda.pip_interop(true)
    Conda.add("python=3.11.6")
    reqtxt = joinpath(Pkg.project().path |> dirname, "requirements.txt")
    Conda.pip("install -r", reqtxt)
    Pkg.build("PyCall")
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
    dicomfolder = match(r"(?<=sdit-)(.+)_ac", pathbasename).captures |> only
    return joinpath(mrn, studyinstanceuid, dicomfolder)
end

function checkforkeyspresent(keystocheck, path)
    searchstring = join(keystocheck, '|')
    return parse(Int, readchomp(`grep -Ec $searchstring $path`))
end

function editjsonsidecar(sidecar, keysfile, zippath, skipexistingkeys = false)
    keystoadd = readlines(keysfile)
    nkeys = length(keystoadd)

    if skipexistingkeys && checkforkeyspresent(keystoadd, sidecar) == nkeys
        return "Skipping $sidecar, keys exist"
    end

    jsonzippath = generatejsonzippath(sidecar)
    searchregex = (Regex ∘ joinpath)(jsonzippath, ".*json")
    mrn = (first ∘ splitpath)(jsonzippath)
    ziparchive = (zip_open_filereader ∘ joinpath)(zippath, mrn * ".zip")

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
        JSON3.pretty(io,
            merge(sidecardict,
                Dict(Symbol(k) => coalesce(v, "missing")
                     for (k, v) in zip(keystoadd, headervalues))))
    end

    indxmissing = BitArray(ismissing.(headervalues))

    if (m = sum(indxmissing)) > 0
        return "Added $(nkeys - m)/$nkeys to $sidecar. Missing $(keystoadd[indxmissing])"
    end

    return "Added all keys to sidecar:$sidecar"
end

function resamplepixdims(inputvolume, outdir, suffix)
    #from fsl.data import image
    #from fsl.utils.image import resample
    fslimage = pyimport("fsl.data.image")
    fslresample = pyimport("fsl.utils.image.resample")
    i = fslimage.Image(inputvolume)
    o = fslresample.resampleToPixdims(i, [1, 1, 1])
    ni = NIVolume(first(o))
    setaffine(ni.header, last(o))
    bn = @chain inputvolume basename replace(_, suffix)
    outfile = joinpath(outdir, bn)
    niwrite(outfile, ni)
    return outfile
end

function resampletoref(inputvolume, lname, reference)
    fslresample = pyimport("fsl.utils.image.resample")
    fslimage = pyimport("fsl.data.image")
    nib = pyimport("nibabel")
    o = fslresample.resampleToReference(inputvolume, fslimage.Image(reference))
    vol = nib.Nifti1Image(o...)
    nib.save(vol, joinpath(dirname(reference), lname * "-mask.nii.gz"))
    # ni = NIVolume(first(o))
    # setaffine(ni.header, last(o))
    # return ni
end

function smoothvolume(inputvolume, outdir; σ = 3.4)
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, "pet" => "8mm_pet"))
    badgerfsl = `/data/h_vmac/vmac_imaging/fsl_v6.0.5.1.sif`
    # cmd = `fslmaths $inputvolume -s 3.4 $outfile`
    cmd = `fslmaths $inputvolume -s $σ $outfile`
    if isfile(replace(string(badgerfsl), "`" => ""))
        run(`singularity exec $badgerfsl $cmd`)
    else
        run(cmd)
    end
    return outfile
end

function getsuvbwscalefactor(sidecarpath)
    json = (JSON3.read ∘ read)(sidecarpath)
    suvbwscalefactor = let seriestime = Time(json.SeriesTime, dateformat"HHMMSS")
        radiopharmaceuticalstarttime = Time(split(json.RadiopharmaceuticalStartTime, '.') |>
                                            first,
            dateformat"HHMMSS")
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
    outpath = joinpath(dirname(inputvolume), replace(basename(inputvolume), suffix))
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
        if isfile(replace(string(badgerfsl), "`" => ""))
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
    deepbet.run_bet([inputvolume],
        [outfile],
        tiv_paths = [tivpath],
        threshold = 0.9,
        n_dilate = -2,
        no_gpu = true)
    return outfile
end

function rigidregistration(fixedvolumepath, movingvolumepath, outdir, suffix)
    bn = basename(movingvolumepath)
    finalmovingvolumepath = replace(bn, suffix)
    finaltransform = replace(finalmovingvolumepath, "nii.gz" => "affine.txt")
    outaffinedir = finaltransform |> splitext |> first
    dipy_align_affine = joinpath(Conda.BINDIR, "dipy_align_affine")
    out_dir = joinpath(outdir, outaffinedir * "_affine")
    run(`$dipy_align_affine $fixedvolumepath $movingvolumepath --transform "affine" --out_dir $out_dir --out_moved $finalmovingvolumepath --out_affine $finaltransform --progressive`)
    return joinpath(out_dir, finalmovingvolumepath)
end

"""
see https://github.com/InsightSoftwareConsortium/ITKElastix/tree/main/examples for itk examples. could probably replace
dipy in `rigidregistration` with itk. Maybe example-18 MONAI_affine_elastix_nonlinear or just equivalent 
"center of mass, translation, rigid body and full affine registration" as performed by dipy_align_affine --progressive
"""
function elastixregistration(fixedvolumepath, movingvolumepath, outdir, suffix)
    bn = basename(movingvolumepath)
    finalmovingvolume = joinpath(outdir, replace(bn, suffix))
    itk_meta_dir = joinpath(outdir, replace(finalmovingvolume, ".nii.gz" => "-elastix_data"))
    if !isdir(itk_meta_dir)
        mkpath(itk_meta_dir)
    end
    itk = pyimport("itk")
    # Load images with itk floats (itk.F). Necessary for elastix
    fixed_image = itk.imread(fixedvolumepath, itk.F)

    moving_image = itk.imread(movingvolumepath, itk.F)
    result_image, result_transform_parameters = itk.elastix_registration_method(fixed_image,
        moving_image,
        output_directory = itk_meta_dir,
        log_file_name = "elastix.log")

    itk.imwrite(result_image, finalmovingvolume)
    setsformqform(finalmovingvolume)
    return finalmovingvolume
end

function generatemasks(refimg, atlases = ["mni", "harvardoxford-subcortical", "harvardoxford-cortical"])
    fslatlases = pyimport("fsl.data.atlases")
    fslatlases.rescanAtlases()
    t = [fslatlases.hasAtlas(atlas) for atlas in atlases]
    if !all(t)
        error("Cannot find atlas: $(atlases[t])")
    end
    labelatlases = [fslatlases.LabelAtlas(fslatlases.getAtlasDescription(at), 1.0)
                    for at in atlases]
    for lat in labelatlases
        for label in lat.desc.labels
            lname = replace(label.name, " " => "_") * "_" * lat.desc.atlasID
            resampletoref(lat.get(label), lname, refimg)
        end
    end
end

function getmeans(suvimgpath, templatedir)
    roimasks = glob("*mask.nii.gz", templatedir)
    suvimg = niread(suvimgpath)
    outfile = @chain suvimgpath basename split(_, ".") first
    mrn = match(r"(?<=sub-)(\d{7,})", suvimgpath).match
    nroi = length(roimasks)
    rowdata = Array{NamedTuple{(:mrn, :label, :mean), Tuple{String, String, Float16}}}(undef,
        nroi)
    for roi in 1:nroi
        maskfile = roimasks[roi]
        label = @chain maskfile basename replace(_,
            ".nii.gz" => "",
            " " => "_",
            "-mask" => "",
            "harvardoxford" => "harvox",
            "cortical" => "cort")
        masknii = niread(maskfile)
        mask = (voxel -> voxel == 1.0 ? voxel : missing).(masknii)
        meanval = (mean ∘ skipmissing)(suvimg .* mask)
        rowdata[roi] = (; mrn, label, mean = meanval)
    end
    df = DataFrame(vcat(rowdata))
    CSV.write(joinpath(dirname(suvimgpath), "$outfile-mean-suv.csv"),
        unstack(df, :label, :mean))
end
