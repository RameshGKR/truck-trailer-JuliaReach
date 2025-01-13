using ReachabilityAnalysis
using Plots

@taylorize function trucktrailer!(du, u, p, t)
    # truck and trailer Parameters
    L0 = 0.3375  # truck length
    M0 = 0.1     # truck distance to center of mass
    L1 = 0.3     # trailer length
    M1 = 0.06    # trailer distance to center of mass

    # state variables: x1 (trailer x), y1 (trailer y), theta0 (truck yaw), theta1 (trailer yaw)
    # control inputs (v0, dtheta0) are now treated as additional states with constant derivatives (i.e., as controls)
    x1, y1, theta0, theta1, v0, dtheta0 = u

    # intermediate variables
    beta01 = theta0 - theta1  # angle between truck and trailer
    v1 = v0 * cos(beta01) + M0 * sin(beta01) * dtheta0  # trailer velocity

    # ODEs + zero(u[.]) for type stability
    du[1] = v1 * cos(theta1) + zero(u[1])
    du[2] = v1 * sin(theta1) + zero(u[2])
    du[3] = dtheta0 + zero(u[3])
    du[4] = (v0 / L1) * sin(beta01) - (M0 / L1) * cos(beta01) * dtheta0 + zero(u[4])
    
    du[5] = zero(u[5])  # v0 remains constant
    du[6] = zero(u[6])  # dtheta0 remains constant

    return du
end

# initial conditions
U₀ = (0.0..0.0) × (0.25..0.25) × (1.57..1.57) × (1.57..1.57) × (0.1..0.3) × (0.2..0.4)  # x1, y1, theta0, theta1, v0, dtheta0

# Define initial value problem
prob = @ivp(u' = trucktrailer!(u), u(0) ∈ U₀, dim:6)

T = 5.0 # Time horizon for the reachability analysis
alg = TMJets21a(; orderT=5, orderQ=1)  #Taylor model-based

# Solve
sol = solve(prob; T=T, alg=alg)

# Overapproximate the solution using zonotopes
solz = overapproximate(sol, Zonotope)

# Plotting
fig = plot(solz; vars=(1, 2), xlab="x1 (trailer x)", ylab="y1 (trailer y)", lw=0.2, color=:blue, lab="Flowpipe")
fig = plot(sol; vars=(1, 2), xlab="x1 (trailer x)", ylab="y1 (trailer y)", lw=0.2, color=:blue, lab="Flowpipe")

plot(solz; vars=(0, 1), xlab="Time (t)", ylab="x1 (trailer x)", lw=0.2, color=:blue, lab="x1(t)")
plot(solz; vars=(0, 2), xlab="Time (t)", ylab="y1 (trailer y)", lw=0.2, color=:blue, lab="y1(t)")

plot(solz; vars=(0, 3), xlab="Time (t)", ylab="theta0 (truck yaw)", lw=0.2, color=:blue, lab="theta0(t)")
plot(solz; vars=(0, 4), xlab="Time (t)", ylab="theta1 (trailer yaw)", lw=0.2, color=:blue, lab="theta1(t)")

plot(solz; vars=(3, 4), xlab="theta0 (truck yaw)", ylab="theta1 (trailer yaw)", lw=0.2, color=:blue, lab="theta1 vs theta0")




