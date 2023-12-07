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

# ----------------------------------------
# Case 1
# ----------------------------------------
println("Running Case 1: a = 1m/s^2; steering angle = 0.2 rad")
prob = f1vehicle()
alg = TMJets(;abstol=1e-7, orderT=5, orderQ=1, adaptive=false)

# warm-up run
sol = solve(prob, T=7.0, alg=alg)

# verify that specification holds
solz = overapproximate(sol, Zonotope)

property = (ρ(ey1, solz) < 4.95) && (ρ(ey2, solz) < 4.95)
push!(validation, Int(property))

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob, T=10.0, alg=$alg)

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=true)

# return the sample with the smallest time value in each test
println("minimum time for each benchmark:\n", minimum(results))

# return the median for each test
println("median time for each benchmark:\n", median(results))

# export runtimes
runtimes = Dict()
for (i, c) in enumerate(cases)
    t = median(results[model][c]).time * 1e-9
    runtimes[c] = t
end

if !@isdefined io
    io = stdout
end

for (i, c) in enumerate(cases)
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c])\n")
end


# ==============================================================================
# Plot
# ==============================================================================

fig = plot()
plot!(fig, solz,  vars=(1, 2), lw=0.0, color=:blue,
    tickfont=font(5, "Times"), guidefontsize=5,
    xlab=L"x",
    ylab=L"y",
    xtick= -2:0.5:15, ytick= -2:0.5:15,
    xlims=(-2.5, 16), ylims=(-2.5, 16),
    bottom_margin=8mm, left_margin=2mm, right_margin=8mm, top_margin=3mm,
    size=(1000, 1000))
hline!(fig, [3.7], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, "F1Tenth_reachability.png")
