using ZipArchives
using Query
using JSON3
using Chain
using Glob
using Dates
using PyCall
using DICOM
using Conda
using NIfTI
using Pkg
using Statistics: mean
using DataFrames
using ImageFiltering: KernelFactors, imfilter
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

function editjsonsidecar(sidecar, keysfile; zippath=nothing, dicom=nothing, skipexistingkeys=false)
    if !isnothing(zippath)
        _editjsonsidecar(sidecar, keysfile, zippath, skipexistingkeys)
    elseif !isnothing(dicom)
        _editjsonsidecardcm(sidecar, keysfile, dicom, skipexistingkeys)
    else
        error("Should not be reached!")
    end
end

function _editjsonsidecardcm(sidecar, keysfile, dicom, skipexistingkeys=false)
    keystoadd = readlines(keysfile)
    nkeys = length(keystoadd)

    if skipexistingkeys && checkforkeyspresent(keystoadd, sidecar) == nkeys
        return "Skipping $sidecar, keys exist"
    end

    dcm = dcm_parse(dicom)


    headervalues = map(keystoadd) do k
        if isnothing(dcm[k]) && occursin("Radiopharmaceutical", k)
            return @chain dcm[:RadiopharmaceuticalInformationSequence] only _[k]
        end

        return dcm[k]
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

function _editjsonsidecar(sidecar, keysfile, zippath, skipexistingkeys=false)
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
end

function getsuvbwscalefactor(sidecarpath)
    json = (JSON3.read ∘ read)(sidecarpath)
    trimmedtime = replace(json.SeriesTime, r"\.\d{1,}" => "")
    suvbwscalefactor = let
        seriestime = Time(trimmedtime, dateformat"HHMMSS")
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

function robustfov(inputvolume, outdir; suffix="pet"=>"cropped_pet")
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, suffix))
    if occursin(r"HN"i, bn)
        cp(inputvolume, outfile; force=true)
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

function maskvolume(inputvolumepath, maskvolumepath, outdir)
    bn = basename(inputvolumepath)
    finalname = replace(bn, ".nii.gz" => "_mask_applied.nii.gz")
    inputvolume = niread(inputvolumepath)
    maskvolume = niread(maskvolumepath)
    finaldata = inputvolume .* maskvolume
    finalwritepath = joinpath(outdir, finalname)
    header = inputvolume.header
    header.scl_inter = 0.0f0  # zero out because slope/intercept auto applied when not using .raw property                                
    niwrite(finalwritepath, NIVolume(header, inputvolume.extensions, finaldata))
    return finalwritepath
end

function skullstrip(inputvolume, outdir, suffix, threshold=0.9, dilation=-2)
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, suffix))
    #tivpath = joinpath(outdir, replace(outfile, "nii.gz" => "TIV.csv"))
    tivpath = replace(outfile, "nii.gz" => "TIV.csv")
    deepbet = pyimport("deepbet")
    deepbet.run_bet([inputvolume],
        [outfile],
        tiv_paths=[tivpath],
        threshold=threshold,
        n_dilate=dilation,
        no_gpu=true)
    return outfile
end

function skullstripwb(inputvolume, outdir, suffix)
    bn = basename(inputvolume)
    outfile = joinpath(outdir, replace(bn, suffix))
    badgersynthstrip = `/data/h_vmac/vmac_imaging/synthstrip.1.5.sif`
    cmd = `mri_synthstrip -i $inputvolume -o $outfile`
    if isfile(replace(string(badgersynthstrip), "`" => ""))
        run(`singularity exec $badgersynthstrip $cmd`)
    else
        run(`/Applications/freesurfer/dev/bin/$cmd`)
    end
    #run(`singularity exec /data/h_vmac/vmac_imaging/synthstrip.1.5.sif mri_synthstrip -i $inputvolume -o $outfile`)
    #run(`/Applications/freesurfer/dev/bin/mri_synthstrip -i $inputvolume -o $outfile`)
    return outfile
end

function applytransform(fixedvolumepath, movingvolumepath, transformfile, transform)
    dipy_apply_transform = joinpath(Conda.BINDIR, "dipy_apply_transform")
    run(`$dipy_apply_transform $fixedvolumepath $movingvolumepath $transformfile --transform_type $transform`)
end

function applytransform(fixedvolumepath, movingvolumepath, transformfile, transform, outdir)
    bn = basename(movingvolumepath)
    outtransformdir = replace(bn, ".nii.gz" => "")
    out_dir = joinpath(outdir, outtransformdir * "_" * transform)
    dipy_apply_transform = joinpath(Conda.BINDIR, "dipy_apply_transform")
    run(`$dipy_apply_transform $fixedvolumepath $movingvolumepath $transformfile --transform_type $transform --out_dir $out_dir`)
end

function applytransform(fixedvolumepath, movingvolumepath, transformfile, transform, outdir, suffix)
    bn = basename(movingvolumepath)
    finalmovingvolumepath = replace(bn, suffix)
    outtransformdir = replace(finalmovingvolumepath, ".nii.gz" => "")
    out_dir = joinpath(outdir, "applytransform_" * transform)
    dipy_apply_transform = joinpath(Conda.BINDIR, "dipy_apply_transform")
    run(`$dipy_apply_transform $fixedvolumepath $movingvolumepath $transformfile --transform_type $transform --out_dir $out_dir --out_file $finalmovingvolumepath`)
    return joinpath(out_dir, finalmovingvolumepath)
end

function rigidregistration(fixedvolumepath, movingvolumepath)
    dipy_align_affine = joinpath(Conda.BINDIR, "dipy_align_affine")
    run(`$dipy_align_affine $fixedvolumepath $movingvolumepath --transform "affine" --progressive`)
    return joinpath(out_dir, finalmovingvolumepath)
end

#=
function rigidregistration(fixedvolumepath, movingvolumepath, outdir)
    bn = basename(movingvolumepath)
    outtransformdir = replace(bn, ".nii.gz" => "")
    out_dir = joinpath(outdir, outtransformdir * "_affine")
    dipy_align_affine = joinpath(Conda.BINDIR, "dipy_align_affine")
    run(`$dipy_align_affine $fixedvolumepath $movingvolumepath --transform "affine"  --out_dir $out_dir --progressive`)
    return joinpath(out_dir, finalmovingvolumepath)   this is bugged here
end
=#

function rigidregistration(fixedvolumepath, movingvolumepath, outdir, suffix, progressive=false)
    bn = basename(movingvolumepath)
    finalmovingvolumepath = replace(bn, suffix)
    finaltransform = replace(finalmovingvolumepath, "nii.gz" => "affine.txt")
    outaffinedir = "affine-registration"
    out_dir = joinpath(outdir, outaffinedir)
    dipy_align_affine = joinpath(Conda.BINDIR, "dipy_align_affine")
    cmdstr = `$dipy_align_affine $fixedvolumepath $movingvolumepath --transform "affine" --out_dir $out_dir --out_moved $finalmovingvolumepath --out_affine $finaltransform`
    if progressive
        cmdstr = `$cmdstr --progressive`
    end
    run(cmdstr)
    return joinpath(out_dir, finalmovingvolumepath)
end

"""
see https://github.com/InsightSoftwareConsortium/ITKElastix/tree/main/examples for itk examples. could probably replace
dipy in `rigidregistration` with itk. Maybe example-18 MONAI_affine_elastix_nonlinear or just equivalent 
"center of mass, translation, rigid body and full affine registration" as performed by dipy_align_affine --progressive
"""
function elastixregistration(fixedvolumepath, movingvolumepath, outdir, suffix, paramfile)
    bn = basename(movingvolumepath)
    finalmovingvolume = joinpath(outdir, replace(bn, suffix))
    itk_meta_dir = joinpath(outdir,
        replace(finalmovingvolume, ".nii.gz" => "-elastix_data"))
    if !isdir(itk_meta_dir)
        mkpath(itk_meta_dir)
    end
    itk = pyimport("itk")
    # Load images with itk floats (itk.F). Necessary for elastix
    fixed_image = itk.imread(fixedvolumepath, itk.F)

    moving_image = itk.imread(movingvolumepath, itk.F)
    parameter_object = itk.ParameterObject.New()
    parameter_object.AddParameterFile(paramfile)
    result_image, result_transform_parameters = itk.elastix_registration_method(fixed_image,
        moving_image,
        output_directory=itk_meta_dir,
        parameter_object=parameter_object,
        #log_file_name="elastix.log")
        log_to_console=false)

    itk.imwrite(result_image, finalmovingvolume)
    setsformqform(finalmovingvolume)
    return (finalmovingvolume, result_transform_parameters)
end

function transformix(movingvolumepath, result_transform_parameters, outdir, suffix)
    bn = basename(movingvolumepath)
    finalmovingvolume = joinpath(outdir, replace(bn, suffix))
    itk_meta_dir = joinpath(outdir,
        replace(finalmovingvolume, ".nii.gz" => "-transformix_data"))
    if !isdir(itk_meta_dir)
        mkpath(itk_meta_dir)
    end
    itk = pyimport("itk")
    # Load images with itk floats (itk.F). Necessary for elastix
    moving_image = itk.imread(movingvolumepath, itk.F)
    result_image_transformix = itk.transformix_filter(
        moving_image,
        result_transform_parameters, output_directory=itk_meta_dir,
        log_file_name="transformix.log")
    itk.imwrite(result_image_transformix, finalmovingvolume)
    setsformqform(finalmovingvolume)
    return finalmovingvolume
end


function generatemasks(refimg,
    atlases=["mni", "harvardoxford-subcortical", "harvardoxford-cortical"])
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

function calculatemetaroistatistic(suvimagpath, templatedir, σ=2.97)
    _metaroimask = glob("Composite-mask*", templatedir) |> only |> niread;
    _referenceregionmask = glob("ponsvermis_1mm_bin-mask*", templatedir) |> only |> niread;
    metaroimask = (voxel -> voxel == 1.0 ? voxel : missing).(vox(_metaroimask, :, :, :))
    referenceregionmask = (voxel -> voxel == 1.0 ? voxel : missing).(vox(_referenceregionmask, :, :, :))
    suvimage = niread(suvimagpath)
    suvmetaroimean =  (mean ∘ skipmissing)(suvimage .* metaroimask)
    suvreferenceregionmean =  @chain (suvimage .* referenceregionmask) skipmissing collect partialsort(_, 1:(length(_)÷2), rev=true) mean
    return suvmetaroimean/suvreferenceregionmean
end

function getmeans(suvimgpath, templatedir, σ=2.97)
    gausskern = KernelFactors.gaussian((σ, σ, σ))
    roimasks = filter(glob("*mask.nii.gz", templatedir)) do m
        !occursin("ponsvermis_1mm_bin-mask", m)
    end
    suvimg = niread(suvimgpath)
    outfile = @chain suvimgpath basename split(_, ".") first
    mrn = match(r"(?<=sub-)(\d{7,})", suvimgpath).match
    nroi = length(roimasks) + 1  # account for adding metaroi statistic after loop
    rowdata = Array{NamedTuple{(:mrn, :label, :mean),Tuple{String,String,Float16}}}(undef,
        nroi)
    for roi in 1:nroi - 1
        maskfile = roimasks[roi]
        label = @chain maskfile basename replace(_,
            ".nii.gz" => "",
            " " => "_",
            "-mask" => "",
            "harvardoxford" => "harvox",
            "cortical" => "cort")
        masknii = niread(maskfile)
        maskwithmissing = (voxel -> voxel == 1.0 ? voxel : missing).(vox(masknii, :, :, :))
        maskedsuvimg = (suvimg .* masknii)
        #smoothedmasked = imfilter(maskedsuvimg, gausskern)
        #niwrite(joinpath(homedir(), "projects/PET-LT/bidsdir/derivatives/roiouts-noct/smoothedout$roi.nii.gz"), NIVolume(suvimg.header, suvimg.extensions, smoothedmasked))
        #meanval = (mean ∘ skipmissing)(smoothedmasked .* maskwithmissing)
        meanval = (mean ∘ skipmissing)(maskedsuvimg .* maskwithmissing)
        rowdata[roi] = (; mrn, label, mean=meanval)
    end

    rowdata[nroi] = (; mrn, label = "metaroistat", mean = calculatemetaroistatistic(suvimgpath, templatedir))

    df = DataFrame(vcat(rowdata))
    outcsv = joinpath(dirname(suvimgpath), "$outfile-mean-suv.csv")
    CSV.write(outcsv,
        unstack(df, :label, :mean))
    return outcsv
end

# /Users/schwartz/projects/PET-LT/

function defaultstructure(datadir)
    dicoms = []
    tostructure = abspath(datadir)
    for (parent, _, files) in walkdir(tostructure)
        if !isempty(files)
            paths = filter(joinpath.(parent, files)) do f
                !occursin("DICOMDIR", f) && !occursin("Viewer", f) && DICOM.isdicom(f)
            end
            push!(dicoms, paths...)
        end
    end

    arr = Array{NamedTuple{(:studyinstanceuid, :modalitysopinstanceuid, :ogpath),Tuple{String,String,String}}}(undef, length(dicoms))

    for (i, path) in enumerate(dicoms)
        dcm = dcm_parse(path)
        if isnothing(dcm.RequestAttributesSequence) || !isnothing(dcm.RequestAttributesSequence[1].StudyInstanceUID) || isempty(something(dcm.RequestAttributesSequence[1].StudyInstanceUID, []))
            siu = dcm.StudyInstanceUID
        else
            siu = dcm.RequestAttributesSequence[1].StudyInstanceUID
        end
        if siu isa Vector
            println(path)
        end
        msiu = dcm.Modality * "." * dcm.SOPInstanceUID * "." * (isempty(dcm.InstanceNumber) ? basename(path) : string(dcm.InstanceNumber))
        if isnothing(siu)
            println(path)
        end
                arr[i] = (; studyinstanceuid=siu, modalitysopinstanceuid=msiu, ogpath=path)
    end
   
    mkpath("stagingdir")

    for row in arr
        mkpath(joinpath("stagingdir", row.studyinstanceuid))
        dest = joinpath("stagingdir", row.studyinstanceuid, row.modalitysopinstanceuid)
        if isfile(dest)
            @warn "DUPLICATE: $(row.ogpath) $dest"
            open("log.txt", "a") do l
                println(l, "$(row.ogpath) $dest")
            end
            continue
        end
        cp(row.ogpath, dest)
    end

end

# defaultstructure("PET-LT-DATA"); process_data_dir("stagingdir", "bidsdir", "")

function windowimage(inputvolumepath, windowcenter, windowwidth, makemask=false)
    imgmin = windowcenter - windowwidth ÷ 2
    imgmax = windowcenter + windowwidth ÷ 2
    inputvolume = niread(inputvolumepath)
    windowimage = vox(inputvolume, :)
    windowimage[windowimage.<imgmin] .= imgmin
    windowimage[windowimage.>imgmax] .= imgmax
    if makemask
        windowimage[(windowimage.>=100).|(windowimage.<0)] .= 0
        windowimage[windowimage.>0] .= 1
        suffix = ".nii.gz" => "_BINARY_MASK.nii.gz"
    else 
        suffix = ".nii.gz" => "_WINDOWED.nii.gz"
    end
    header = inputvolume.header
    header.scl_inter = 0.0f0  # zero out because we already applied slope/intercept by accessing data w/ `vox`                                
    outpath = replace(inputvolumepath, suffix)
    niwrite(outpath, NIVolume(header, inputvolume.extensions, reshape(windowimage, size(inputvolume))))
    return outpath
end

function addinfotocsv(sidecar, csv)
    df = DataFrame(CSV.File(csv))
    js = JSON3.read(read(sidecar))
    df[:, "PatientAge"] .= js.PatientAge
    df[:, "AccessionNumber"] .= parse(Int,js.AccessionNumber)
    df[:, "AdditionalPatientHistory"] .= something(js.AdditionalPatientHistory, " ") |> strip |> (s->replace(s, '\n'=>' '))
    df[:, "AcquisitionDateTime"] .= Date(match(r"(.*:\d{2})", js.AcquisitionDateTime).captures |> only, dateformat"YYYY-m-dTHH:MM:SS")
    df[:, "ReasonForStudy"] .= @chain js.ReasonForStudy begin 
    something(_, missing)
    if !ismissing(_) 
        strip  
        replace(_, "\n"=> " ")
    else
         missing
    end
end
    l =[:mrn, :AccessionNumber, :AcquisitionDateTime, :PatientAge, :AdditionalPatientHistory, :ReasonForStudy]
    CSV.write(csv, df[:, Cols(l..., Not(l))]; transform=(col, val) -> something(val, missing))
end
