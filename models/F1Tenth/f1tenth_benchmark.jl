using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "F1Tenth"
cases = [""]
SUITE[model] = BenchmarkGroup()

include("f1tenth.jl")
validation = []

LazySets.deactivate_assertions()


println("Running Different discritisation of Input")

prob = f1vehicle()
alg = TMJets(;abstol=1e-7, orderT=5, orderQ=1, adaptive=false)

# warm-up run
sol = solve(prob, T=3.0, alg=alg)

# verify that specification holds
solz = overapproximate(sol, Zonotope)
vlim = 7;
deltalim = 0.41;
property = (ρ(ey1, solz) < vlim) && (ρ(-1*ey1, solz) < vlim) && (ρ(ey2, solz) < deltalim) && (ρ(-1*ey2, solz) < deltalim)
push!(validation, Int(property))

# # benchmark
# SUITE[model][cases[1]] = @benchmarkable solve($prob, T=3.0, alg=$alg)

# # ==============================================================================
# # Execute benchmarks and save benchmark results
# # ==============================================================================

# # tune parameters
# tune!(SUITE)

# # run the benchmarks
# results = run(SUITE, verbose=true)

# # # return the sample with the smallest time value in each test
# # println("minimum time for each benchmark:\n", minimum(results))

# # # return the median for each test
# # println("median time for each benchmark:\n", median(results))

# # export runtimes
# runtimes = Dict()
# for (i, c) in enumerate(cases)
#     t = median(results[model][c]).time * 1e-9
#     runtimes[c] = t
# end

# if !@isdefined io
#     io = stdout
# end

# for (i, c) in enumerate(cases)
#     print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c])\n")
# end


# ==============================================================================
# Plot
# ==============================================================================

fig = plot()
plot!(fig, solz,  vars=(1, 2), lw=0.0, color=:blue,
    tickfont=font(5, "Times"), guidefontsize=5,
    xlab=L"x",
    ylab=L"y",
    xtick= -1:5:100, ytick= -1:5:100,
    xlims=(-2, 120), ylims=(-2, 120),
    bottom_margin=8mm, left_margin=2mm, right_margin=8mm, top_margin=3mm,
    size=(1000, 1000))
hline!(fig, [3.7], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, "F1Tenth_reachability.png")
