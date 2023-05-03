using Tar
using JSON3
using DelimitedFiles
using Distributed


function constructargs(datadir, tarballdir, bidsjson)
    """need study and series dir for regex"""
    j = JSON3.read(read(bidsjson, String))
    study = j[:StudyInstanceUID]
    mrn = j[:PatientID]
    seriesdir = match(r"(?<=sdit-).*(?=_ac.*)", bidsjson).match
    pattern = Regex(".*/$study/$seriesdir/.*json")
    outdir = ((x -> chopsuffix(x, ".json")) ∘ basename)(bidsjson)
    return (pattern, joinpath(abspath(tarballdir), mrn * ".tar.gz"), joinpath(abspath(datadir), "derivatives", outdir))
end

function tarextractor(args)
    function _tarextractor(pattern, tarball, path)
        try
            toextract = first(filter(hdr -> occursin(pattern, hdr.path), Tar.list(`gunzip -dc $tarball`))).path
            Tar.extract(x -> (x.path == toextract), `gunzip -dc $tarball`, path)
        catch e
            display(e)
            display(tarball)
        end
    end

    _tarextractor(args...)
end

function tarextract(datadir, tarballdir)
    """
    expects datadir to be top level data dir
   """

    cmd = `find $datadir -name "sub*json"`
    jsonpaths = readlines(Cmd(cmd))
    Distributed.pmap(tarextractor, constructargs.(datadir, tarballdir, jsonpaths))
end