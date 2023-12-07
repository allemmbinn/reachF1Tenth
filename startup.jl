# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

const TARGET_FOLDER = "result"
const RESULTS_FILE = "results.csv"

function main()
    if !isdir(TARGET_FOLDER)
        mkdir(TARGET_FOLDER)
    end
    global io = open(joinpath(TARGET_FOLDER, RESULTS_FILE), "w")

    println("Running NLN benchmarks...")

    # # Robertson benchmark (ROBE21)
    # println("###\nRunning Robertson benchmark\n###")
    # include("models/Robertson/robertson_benchmark.jl")

    # # Coupled Van der Pol benchmark (CVDP22)
    # println("###\nRunning Van der Pol benchmark\n###")
    # include("models/VanDerPol/vanderpol_benchmark.jl")

    # # Laub-Loomis benchmark (LALO20)
    # println("###\nRunning Laub-Loomis benchmark\n###")
    # include("models/LaubLoomis/laubloomis_benchmark.jl")

    # # Lotka-Volterra tangential crossing benchmark (LOVO21)
    # println("###\nRunning Lotka-Volterra tangential crossing benchmark\n###")
    # include("models/LotkaVolterra/lotka_volterra_benchmark.jl")

    # # Spacecraft benchmark (SPRE22)
    # println("###\nSpacecraft benchmark\n###")
    # include("models/Spacecraft/spacecraft_benchmark.jl")

    # F1Tenth
    println("###\nF1Tenth benchmark\n###")
    include("models/F1Tenth/f1tenth_benchmark.jl")

    print(io, "\n")
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
