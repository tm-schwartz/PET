using Tar
using JSON3
using DelimitedFiles


function constructargs(bidsjson)
    """need study and series dir for regex"""
    j = JSON3.read(read(bidsjson, String))
    study = j[:StudyInstanceUID]
    mrn = j[:PatientID]
    seriesdir = match(r"(?<=sdit-).*(?=_ac.*)", bidsjson).match
    pattern = Regex(".*/$study/$seriesdir/.*json")
    outdir = ((x -> chopsuffix(x, ".json")) ∘ basename)(bidsjson)
    return (pattern, outdir, mrn,)
end

function tarextractor(pattern, tarball, path)
    toextract = first(filter(hdr -> occursin(pattern, hdr.path), Tar.list(`zcat $tarball`))).path
    Tar.extract(x -> (x.path == toextract), `zcat $tarball`, path)
end

function tarextract(datadir, tarballdir)
    """
    expects datadir to be top level data dir
   """
    cmd = `find $datadir -name 'sub*.json'`
    jsonpaths = Cmd(cmd)
    for (pattern, outdir, mrn) ∈ constructargs.(jsonpaths)
        tarextractor(pattern, joinpath(datadir, "derivatives", outdir), joinpath(tarballdir, mrn * ".tar.gz"))
    end
end