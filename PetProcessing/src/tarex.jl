using Tar

function tarextractor(pattern, tarball, path)
    toextract = first(filter(hdr->occursin(pattern, hdr.path), Tar.list(`gzcat $tarball`))).path       
    Tar.extract(x->(x.path == toextract), `gzcat $tarball`, path)
end
