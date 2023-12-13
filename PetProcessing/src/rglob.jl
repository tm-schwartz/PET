function rglob(dir, PATTERN)
    matches = String[]

    for entry in readdir(dir)
        entry_path = joinpath(dir, entry)

        if isdir(entry_path)
            # If the entry is a directory, recursively call the function
            subdirectory_matches = rglob(entry_path, PATTERN)
            append!(matches, subdirectory_matches)
        elseif occursin(PATTERN, entry)
            # If the entry matches the pattern, add it to the matches
            push!(matches, entry_path)
        end
    end

    return matches
end