module FdgPet
export process_data_dir
using DICOM
using NIfTI
using Distributed
using UUIDs
using JSON3
using Makie, CairoMakie

const excluded = r"^((?!vumcct).)*$|mip|leg[s]?|nac|screen|ds\_store|fused|statistic[s]?"i

subfiles = startswith("sub-")
exclude(input::Vector{String}) = filter(fn -> !occursin(excluded, fn), input)
exclude(input::String) = occursin(excluded, input)


function peek_tag(dicom_path, tag)
    dcm = dcm_parse(dicom_path)
    if haskey(dcm, tag)
        return dcm[tag"ImageIndex"]
    else
        return missing
    end
end

function dcm_fix_name(dir)
    paths = joinpath.(dir, readdir(dir))
    dicoms = filter(DICOM.isdicom, paths)

    updated_paths = map(dicoms) do dicom
        ext = splitext(dicom)[end]
        updated_fn = dicom
        has_imin = occursin(r"imin-\d{3}", dicom)
        imin = peek_tag(dicom, "ImageIndex")
        if has_imin == false & isempty(ext) || ext ≠ ".dcm"
            updated_fn = updated_fn * "-" * (ismissing(imin) ? "1" : string(imin)) * ".dcm"
        elseif has_imin == false & ext == ".dcm"
            updated_fn = splitext(updated_fn)[1] * "-" * (ismissing(imin) ? "1" : string(imin)) * ext
        elseif has_imin & ext ≠ ".dcm"
            updated_fn *= ".dcm"
        end
        return updated_fn
    end
    for (path, newpath) in zip(dicoms, updated_paths)
        mv(path, newpath)
    end
end

function process_descriptor(descriptor::AbstractString)
    re = r"((?<=\W)|^)(\w{1,})(?=\W?)"
    matches = findall(re, descriptor)
    return join([descriptor[m] for m in matches], "_")
end

function process_descriptor(descriptor::Vector)
    desc = join(descriptor, "_")
    return replace(desc, " " => "_")  # account for spaces in vector elements
end

function dcm_dir_to_series(dir)
    dicom_files = DICOM.find_dicom_files(dir)
    series_dict = Dict()
    for dicom in dicom_files
        dcm = dcm_parse(dicom)
        seriesdescr = replace(dcm[tag"SeriesDescription"], " " => "_", r"\W" => "")
        imagetype = process_descriptor(dcm[tag"ImageType"])
        if :StationName in propertynames(dcm)
            stationname = process_descriptor(dcm[tag"StationName"])
        else
            stationname = "NONE"
        end
        seriesdata = join([seriesdescr, imagetype, stationname], "_")
        series_path = joinpath(dir, seriesdata)
        dicom_currentpath = joinpath(dir, dicom)
        if haskey(series_dict, series_path)
            push!(series_dict[series_path], dicom_currentpath)
        else
            series_dict[series_path] = String[dicom_currentpath]
        end
    end
    return series_dict
end

function create_softlinks(target, img_dir)
    link = joinpath(img_dir, basename(target))
    if islink(link)
        lst = splitext(link)
        uuid = UUIDs.uuid4()
        link = string(lst[1], uuid, lst[2])
    end
    symlink(target, link)
end


function dcm2niix(input_dir::String, pyenv::String)
    command = `sh -c "source $pyenv; dcm2niix -a y -m y -b y -ba n -z y -f sub-%i_sdit-%f_ac-%g_dt-%t_pet $input_dir"`
    try
        output = string(readlines(Cmd(command))...)
    catch e
        display(e)
    end
end

function _slicepng(input_path)
    png_name = (first ∘ splitext ∘ first ∘ splitext ∘ basename)(input_path)
    png_location = joinpath(dirname(input_path), "png-dir")
    try
        nii = niread(input_path)
        if ndims(nii) < 3
            return
        end
        if !isdir(png_location)
            mkdir(png_location)
        end
        nii_size = size(nii)
        middle_slice = nii_size[1] ÷ 2
        resolution = nii_size[end-1:end]
        scene = Scene(resolution=resolution)
        campixel!(scene)
        heatmap!(scene, nii.raw[middle_slice, :, :], colormap=:grayC)
        save(joinpath(png_location, png_name * ".png"), scene)
    catch e
        display(e)
        display(input_path)
        return
    end
end

function slicepng(data_dir, img_link_dir)
    studypaths = filter(isdir, readdir(data_dir, join=true))
    seriespaths = exclude(filter(isdir, append!(map(sp -> readdir(sp, join=true), studypaths)...)))
    niftipaths = append!(map(seriespaths) do sp
        return filter(np -> occursin(r"\.nii\.gz$", np), readdir(sp, join=true))
    end...)
    cleaned_niftis = exclude(niftipaths)
    Distributed.pmap(_slicepng, cleaned_niftis)
    map(seriespaths) do spath
        pth = joinpath(spath, "png-dir")
        if isdir(pth)
            fl = readdir(pth, join=true)
            for target in fl
                create_softlinks(target, img_link_dir)
            end
        end
    end
    return
end

function tobids(data_dir, bidsdir, modalityfolder)
    niftifiles = readlines(Cmd(`find $data_dir -name 'sub*gz' -or -name 'sub*json'`))
    if isempty(niftifiles)
        return
    end
    subdir = niftifiles |> first |> basename |> (x -> split(x, "_")) |> first |> (x -> joinpath(bidsdir, x, modalityfolder))
    if !isdir(subdir)
        mkpath(subdir)
    end
    map(niftifiles) do f
        dest = joinpath(subdir, basename(f))
        try
            mv(f, dest)
        catch y
            display(y)
            display(f)
        end
    end
    return
end

function process_data_dir(data_dir, bidsdir, modalityfolder, pyenv)
    contents = readdir(data_dir)
    studies = filter(d -> isdir(joinpath(data_dir, d)), contents)
    img_links = joinpath(data_dir, "image_links")
    mkdir(img_links)
    for study in studies
        study_path = joinpath(data_dir, study)
        series_dict = dcm_dir_to_series(study_path)
        for (seriespath, dicompath) in series_dict
            absolute_seriespath = joinpath(data_dir, seriespath)
            mkpath(absolute_seriespath)
            for dicom in dicompath
                fn = splitdir(dicom)[2]
                original_path = joinpath(data_dir, study, fn)
                new_path = joinpath(data_dir, seriespath, fn)
                mv(original_path, new_path)
                json_file = string(original_path, ".json")
                if isfile(json_file)
                    new_json_path = string(new_path, ".json")
                    mv(json_file, new_json_path)
                end
            end
            if exclude(absolute_seriespath)
                continue
            else
                dcm_fix_name(absolute_seriespath)
                dcm2niix(absolute_seriespath, pyenv)
            end
        end
    end
    niftifiles = readlines(Cmd(`find $data_dir -name 'sub*gz' -or -name 'sub*json'`))
    if isempty(niftifiles)
        return
    else
        slicepng(data_dir, img_links)
        tobids(data_dir, bidsdir, modalityfolder)
    end
end

end # module