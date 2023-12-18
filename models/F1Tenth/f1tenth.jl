using ReachabilityAnalysis
using ReachabilityAnalysis: TaylorModelReachSet, AbstractLazyReachSet

vlim = 7.0;
deltalim = 0.41;
ratlim = deltalim/vlim;

# Constraint for vx
const ey1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] 
# Constraint for vx
const ey2 = [0.0, 0.0, 0.0, ratlim, 0.0, 0.0, 1.0]


@taylorize function f1vehicle!(dx, x, p, t)
    # x[9]  -> Frx
    # x[10] -> dsteer
    lr = 0.17145
    lf = 0.15875
    # m  = 3.47

    dx[1] = x[4] * cos(x[3]) - x[4] * sin(x[3])
    dx[2] = x[4] * sin(x[3]) + x[4] * cos(x[3])
    dx[3] = x[6]
    dx[4] = p[1]
    dx[7] = p[2]
    dx[5] = lr/(lr+lf) * (x[7]*dx[4] + dx[7]*x[4])
    dx[6] = 1 /(lr+lf) * (x[7]*dx[4] + dx[7]*x[4])
    return dx
end

function f1vehicle()
    # X0 = (1.25 .. 1.55) × (2.35 .. 2.45) × (1.25 .. 1.55) × (2.35 .. 2.45) × b
    X0 = Hyperrectangle([0.0, 0.0, 0.0, 1.0, 0.0, 0., 0.],
                        [0.3, 0.3, 0.001, 0.1, 0.1, 0.001, 0.001])
    U0 = ZeroSet(2);
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