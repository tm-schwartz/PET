using DataFrames, CSV, Query, Glob, Chain, Dates, JSON3, DataFramesMeta
scans = glob("sub-*/pet/*.nii.gz", "")
s = scans .|> splitpath .|> parts -> parts[1:2:end]
subs = @. @chain s _|>first split(_, "-")|>last parse(Int, _)
files = s .|> last
subscans = DataFrame("sub" => subs, "scan" => files)
subscans[:, "accession"] = @chain subscans.scan match.(r"ac-(\d{1,})_dt", _) broadcast(x -> x.captures,
    _) only.(_) parse.(Int, _)
gd = groupby(subscans, [:sub, :accession])

function procdt(procfile)
    if procfile == "-1"
        return Dates.unix2datetime(Dates.UNIXEPOCH)
    end
    return @chain procfile match.(r"dt-(\d{1,})_pet.?", _).captures only DateTime(string(_),
        dateformat"yyyymmddHHMMS")
end

function procgrp(g)
    if nrow(g) == 1
        return g.scan
    else
        trans = @chain g.scan match.(Regex("(sub-.*?_.*?tran.*?pet.?\\.nii\\.gz)", "i"), _) cat(_...,
            dims = 1) filter(m -> !isnothing(m), _) broadcast(m -> m.captures |> only, _) filter(m -> !occursin(r"2D"i,
                m),
            _)
        axial = @chain g.scan match.(Regex("(sub-.*?_.*?ax.*?pet.?\\.nii\\.gz)", "i"), _) cat(_...,
            dims = 1) filter(m -> !isnothing(m), _) broadcast(m -> m.captures |> only, _) filter(m -> !occursin(r"2D"i,
                m),
            _)
        orig = @chain g.scan match.(Regex("(sub-.*?_.*?original.*?pet.?\\.nii\\.gz)", "i"),
            _) cat(_...,
            dims = 1) filter(m -> !isnothing(m), _) broadcast(m -> m.captures |> only, _) filter(m -> !occursin(r"2D"i,
                m),
            _)
        headneck = @chain g.scan match.(Regex("(sub-.*?_.*?hn.*?pet.?\\.nii\\.gz|sub-.*?_.*?brain.*?pet.?\\.nii\\.gz)",
                "i"),
            _) cat(_...,
            dims = 1) filter(m -> !isnothing(m), _) broadcast(m -> m.captures |> only, _) filter(m -> !occursin(r"2D"i,
                m),
            _)
        topush = ""
        for re in (headneck, orig, axial, trans)
            if !isempty(re)
                if length(re) > 1
                    origprim = filter(m -> occursin(r"original.*primary"i, m), re)
                    primary = filter(m -> occursin(r"primary"i, m), re)
                    ax = filter(m -> occursin(r"ax"i, m), re)
                    sag = filter(m -> occursin(r"sag"i, m), re)
                    if !isempty(origprim)
                        topush = origprim
                    elseif !isempty(primary)
                        topush = primary
                    elseif !isempty(ax)
                        topush = ax
                    else
                        topush = sag
                    end
                else
                    topush = re |> only
                end
                if any(occursin.("Eq_1", re)) && !occursin("Eq_1", topush)
                    totest = replace(topush, ".nii.gz" => "_Eq_1.nii.gz")
                    if occursin.(totest, re)
                        topush = totest
                    end
                end

                if topush isa Vector
                    topush = first(topush)
                end
                topush = repeat([string(topush)], nrow(g))
                break
            end
        end
        if isempty(topush)
            println(g.scan)
            return repeat(["-1"], nrow(g))
        end

        return topush
    end
end

df = combine(gd, :scan, procgrp);
df[:, :procdt] = map(procdt, df.x1);
rename!(df, :x1 => :procfile)
udf = df[:, Not(:scan)] |> unique
udf[:, :reasonforstudy] = @with udf @byrow begin
    if :procfile == "-1"
        return "-1"
    end
    jsnfile = joinpath("sub-" * string(:sub), "pet", :procfile)
    jsnfile = replace(jsnfile, "nii.gz" => "json")
    if occursin("Eq_1", jsnfile)
        jsnfile = replace(jsnfile, "_Eq_1" => "")
    end
    if stat(jsnfile).size > 0
        open(jsnfile, "r") do f
            jsn = JSON3.read(f)
            reason = get(jsn, "ReasonForStudy")
            stri = strip(reason)
            replace(stri, '\n' => ' ')
        end
    else
        return "unspecified"
    end
end

CSV.write("/nobackup/h_vmac/PET/DATASET-20230414/derivatives/processfiles-test.csv", udf)