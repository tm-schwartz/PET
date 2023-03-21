module FdgPet
export process_data_dir
using DICOM

function peek_tag(dicom_path, tag)
    dcm = dcm_parse(dicom_path)
    if haskey(dcm, tag)
        return dcm[tag"ImageIndex"]
    else
        return missing
    end
end

function dcm_fix_name(dir)
    endswithdcm = endswith(".dcm")
    paths = joinpath.(dir, readdir(dir))
    dicoms = filter(DICOM.isdicom, paths)
    updated_paths = map(dicoms) do dicom
        ext = splitext(dicom)[end]
        updated_fn = dicom
        has_imin = occursin(r"imin-\d{3}", dicom)
        imin = peek_tag(dicom, "ImageIndex")
        if (has_imin == false) & (isempty(ext) || ext ≠ ".dcm")
            updated_fn = updated_fn * "-" * (ismissing(imin) ? "001" : string(imin)) * ".dcm"
        elseif (has_imin == false) & (ext == ".dcm")
            updated_fn = splitext(updated_fn)[1] * "-" * (ismissing(imin) ? "001" : string(imin)) * ext
        elseif has_imin & (ext ≠ ".dcm")
            updated_fn *= ".dcm"
        end
        return updated_fn
    end
    for (path, newpath) in zip(dicoms, updated_paths)
        mv(path, newpath)
    end
end

function dcm_dir_to_series(dir)
    dicom_files = DICOM.find_dicom_files(dir)
    series_dict = Dict()
    for dicom in dicom_files
        dcm = dcm_parse(dicom)
        seriesuid = dcm[tag"SeriesInstanceUID"]
        seriesdescr = replace(dcm[tag"SeriesDescription"], " "=>"_")
        seriesdata = string(seriesuid, "_", seriesdescr)
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


function dcm2niix(input_dir)
    command = `dcm2niix -a y -m y -b y -ba n -z o $input_dir`
    output = string(readlines(Cmd(command))...)
    display(output)
end 

function process_data_dir(data_dir)
    contents = readdir(data_dir)
    studies = filter(d -> isdir(joinpath(data_dir, d)), contents)
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
            dcm_fix_name(absolute_seriespath)
        end
    end
end


end # module


# med2image -i 1.2.840.113619.2.80.45947970.24763.1190911157.146.4.4_Reformatted/PT.1.2.840.113619.2.80.45947970.24763.1190911157.147-1.dcm
# -d dicom-results/middle-slice -o sample --outputFileType jpg --sliceToConvert m --verbosity 3 --convertOnlySingleDICOM