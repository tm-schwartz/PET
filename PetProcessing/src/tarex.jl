using Tar

function tarextractor(pattern, tarball, path)
    toextract = first(filter(hdr->occursin(pattern, hdr.path), Tar.list(`zcat $tarball`))).path       
    Tar.extract(x->(x.path == toextract), `zcat $tarball`, path)
end
