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
    return @chain procfile match.(r"dt-(\d{1,})_pet", _).captures only DateTime(string(_),
        dateformat"yyyymmddHHMMS")
end

function procgrp(g)
    if nrow(g) == 1
        return g.scan
    else
        c = join(g.scan)
        o = "original"
        a1 = "ax"
        a2 = "axial"
        t1 = "tran"
        t2 = "trans"
        t3 = "transaxial"
        hn = "hn"
        orig = match(Regex(".*(sub-\\d{1,}.*\\Q$o\\E.*?gz).*", "i"), c)
        axial = match(Regex(".*(sub-\\d{1,}.*(?>\\Q$a1\\E|\\Q$a2\\E).*?gz).*", "i"), c)
        trans = match(Regex(".*(sub-\\d{1,}.*(?>\\Q$t1\\E|\\Q$t2\\E|\\Q$t3\\E).*?gz).*",
                "i"),
            c)
        headneck = match(Regex(".*(sub-\\d{1,}.*\\Q$hn\\E.*?gz).*", "i"), c)
        if !isnothing(orig) && !occursin(r"2D"i, string(orig.captures))
            topush = orig.captures |> only
            if occursin("Eq_1", c) && !occursin("Eq_1", topush)
                totest = replace(topush, ".nii.gz" => "_Eq_1.nii.gz")
                if occursin(totest, c)
                    topush = totest
                end
            end
            return repeat([string(topush)], nrow(g))
        elseif !isnothing(axial) && !occursin(r"2D"i, string(axial.captures))
            topush = axial.captures |> only
            if occursin("Eq_1", c) && !occursin("Eq_1", topush)
                totest = replace(topush, ".nii.gz" => "_Eq_1.nii.gz")
                if occursin(totest, c)
                    topush = totest
                end
            end
            return repeat([string(topush)], nrow(g))
        elseif !isnothing(trans) && !occursin(r"2D"i, string(trans.captures))
            topush = trans.captures |> only
            if occursin("Eq_1", c) && !occursin("Eq_1", topush)
                totest = replace(topush, ".nii.gz" => "_Eq_1.nii.gz")
                if occursin(totest, c)
                    topush = totest
                end
            end
            return repeat([string(topush)], nrow(g))
        elseif !isnothing(headneck) && !occursin(r"2D"i, string(headneck.captures))
            topush = headneck.captures |> only
            if occursin("Eq_1", c) && !occursin("Eq_1", topush)
                totest = replace(topush, ".nii.gz" => "_Eq_1.nii.gz")
                if occursin(totest, c)
                    topush = totest
                end
            end
            return repeat([string(topush)], nrow(g))
        else
            println(c)
            sub = split(c, "_") |> first
            return repeat(["-1"], nrow(g))
        end
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
            stri=strip(reason)
            replace(stri, '\n'=>' ')
        end
    else
        return "unspecified"
    end
end

CSV.write("/nobackup/h_vmac/PET/DATASET-20230414/derivatives/processfiles.csv", udf)