using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "CVDP22"
cases = [""]
SUITE[model] = BenchmarkGroup()

include("vanderpol.jl")
validation = []

LazySets.deactivate_assertions()

# ----------------------------------------
# Case 1: μ = 1, b = [60, 80]
# ----------------------------------------

prob = vanderpolN2_b(b=interval(60.0, 80.0))
alg = TMJets(abstol=1e-10, orderT=8, orderQ=1)

# warm-up run
sol = solve(prob, T=7.0, alg=alg)

# verify that specification holds
solz = overapproximate(sol, Zonotope)
property = (ρ(ey₁, solz) < 3.7) && (ρ(ey₂, solz) < 3.7)
push!(validation, Int(property))

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob, T=7.0, alg=$alg)

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
plot!(fig, solz,  vars=(1, 2), lw=0.0, alpha=1.0, color=:blue,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"x_{1}",
    ylab=L"y_1",
    xtick=[-2.0, 0.0, 2.0], ytick=[-4.0, -2.0, 0.0, 2.0, 4.0],
    xlims=(-2.5, 2.5), ylims=(-4.05, 4.05),
    bottom_margin=8mm, left_margin=2mm, right_margin=8mm, top_margin=3mm,
    size=(1000, 1000))
hline!(fig, [3.7], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, "ARCH-COMP22-JuliaReach-CVDP22.png")
