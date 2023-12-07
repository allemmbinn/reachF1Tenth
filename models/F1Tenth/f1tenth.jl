using ReachabilityAnalysis
using ReachabilityAnalysis: TaylorModelReachSet, AbstractLazyReachSet

const ey1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] 
const ey2 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] 

@taylorize function f1vehicle!(dx, x, p, t)
    lr = 0.17145
    lf = 0.15875

    dx[1] = x[4] * cos(x[3]) - x[4] * sin(x[3])
    dx[2] = x[4] * sin(x[3]) + x[4] * cos(x[3])
    dx[3] = x[6]
    dx[4] = one(x[1]) # Constant Accelartion
    dx[5] = lr/(lr+lf) * (x[7])
    dx[6] = 1 /(lr+lf) * (x[7])
    dx[7] = zero(x[1]) # Constant Steering Angle
    return dx
end

function f1vehicle()
    # X0 = (1.25 .. 1.55) × (2.35 .. 2.45) × (1.25 .. 1.55) × (2.35 .. 2.45) × b
    X0 = Hyperrectangle([0.0, 0.0, 0.0, 0.0, 0., 0., 0.05],
                        [0.3, 0.3, 0., 0.1, 0.1, 0., 0.])
    @ivp(x' = f1vehicle!(x), dim: 7, x(0) ∈ X0)
end

# function absolute_velocity(R::AbstractLazyReachSet)
#     vx = 4 # Position of vx in state
#     vy = 5 # Position of vy in state
#     vx2 = set(overapproximate(project(R, vars=(vx,)), Interval)).dat
#     vy2 = set(overapproximate(project(R, vars=(vy,)), Interval)).dat
#     sqrt(vx2 + vy2)
# end

# function absolute_velocity(R::TaylorModelReachSet)
#     Z = overapprximate(R, Zonotope)
#     absolute_velocity(Z)
# end

# function velocity_constraint(sol)
#     all_idx = findall(x -> x == 2, location.(sol)) # attempt
#     for idx in all_idx
#         # maximum velocity measured in m/min
#         verif = all(absolute_velocity(R) < 0.055 * 60. for R in sol[idx])
#         !verif && return false
#     end
#     return true
# end