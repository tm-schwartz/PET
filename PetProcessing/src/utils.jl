using PackageCompiler

function compilepackage()
    cmd = `julia --project=PetProcessing --trace-compile=pettrace.jl -e "push!(LOAD_PATH, \"PetProcessing\"); using PetProcessing; exit()"`
    run(cmd)
    create_sysimage(["PetProcessing"]; sysimage_path="PetProcessing-img.so", precompile_statements_file="pettrace.jl")
end