module PetProcessing

export
    # bids
    process_data_dir,
    
    # tarex
    tarextract,
    
    # utils
    initializepythonenv,
    editjsonsidecar,
    getsuvbwscalefactor,
    robustfov,
    resamplepixdims,
    rigidregistration,
    skullstrip,
    elastixregistration,
    smoothvolume,
    computesuvvolume,
    getmeans

include("bids.jl")
include("tarex.jl")
include("utils.jl")

end # module