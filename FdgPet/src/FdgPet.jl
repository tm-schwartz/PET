module FdgPet
export process_data_dir
using DICOM

function dcm_add_ext(dir)
    endswithdcm = endswith(".dcm")
    filetype = "DICOM medical imaging data"
    for dicom in readdir(dir)
        path = joinpath(dir, dicom)
        command = `file -b $path`
        output = readlines(Cmd(command))[1]
        if !endswithdcm(dicom) && output == filetype
            mv(path, string(path, ".dcm"))
        end
    end
end

function dcm_dir_to_series(dir)
    dicom_files = DICOM.find_dicom_files(dir)
    series_dict = Dict()
    for dicom in dicom_files
        dcm = dcm_parse(dicom)
        series_path = joinpath(dir, dcm[tag"SeriesInstanceUID"])
        dicom_cpath = joinpath(dir, dicom)
        if haskey(series_dict, series_path)
            push!(series_dict[series_path], dicom_cpath)
        else
            series_dict[series_path] = String[dicom_cpath]
        end
    end
    return series_dict
end

function process_data_dir(data_dir)
    contents = readdir(data_dir)
    studies = filter(d -> isdir(joinpath(data_dir, d)), contents)
    println(studies)
    for study in studies
        study_path = joinpath(data_dir, study)
        series_dict = dcm_dir_to_series(study_path)
        println(study_path)
        for (seriespath, dicompath) in series_dict
            absolute_seriespath = joinpath(data_dir, seriespath)
            mkpath(absolute_seriespath)
            println(absolute_seriespath)
            for dicom in dicompath
                fn = splitdir(dicom)[2]
                original_path = joinpath(data_dir, study, fn)
                new_path = joinpath(data_dir, seriespath, fn)
                mv(original_path, new_path)
                json_file = string(original_path, ".json")
                println(json_file)
                if isfile(json_file)
                    new_json_path = string(new_path, ".json")
                    mv(json_file, new_json_path)
                end
            end
            dcm_add_ext(absolute_seriespath)
        end
    end
end

end # module
