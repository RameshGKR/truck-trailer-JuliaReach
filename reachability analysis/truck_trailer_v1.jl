# Truck-Trailer Model using JuliaReach
using ReachabilityAnalysis, LazySets, MathematicalSystems

# Model Parameters
const L1 = 0.3  # Length of the trailer
const M0 = 0.1  # Parameter related to trailer dynamics

@taylorize function truck_trailer!(du, u, p, t)
    # States: theta1, x1, y1, theta0
    theta1, x1, y1, theta0 = u

    # Control inputs: v0 (linear velocity of truck), dtheta0 (angular velocity of truck)
    v0, dtheta0 = p

    # Derived variables
    beta01 = theta0 - theta1
    dtheta1 = v0 / L1 * sin(beta01) - M0 / L1 * cos(beta01) * dtheta0
    v1 = v0 * cos(beta01) + M0 * sin(beta01) * dtheta0

    # Define the ODEs
    du[1] = dtheta1               # d(theta1)/dt
    du[2] = v1 * cos(theta1)      # d(x1)/dt
    du[3] = v1 * sin(theta1)      # d(y1)/dt
    du[4] = dtheta0               # d(theta0)/dt

    return du
end

# Define initial set as a hyperrectangle
theta1_0 = -0.2    
x1_0 = 0.0         
y1_0 = 2.5         
theta0_0 = -0.2    

# Define initial state range
X0 = Hyperrectangle([theta1_0, x1_0, y1_0, theta0_0], [0.05, 0.1, 0.1, 0.05])

# Define control parameters (as constants for now, or use ranges)
v0 = 0.1
dtheta0 = 0.01

# Define the continuous system using @system macro
truck_trailer_system = @system(x' = truck_trailer!(x, [v0, dtheta0], 0.0), dim:4)

# Define the initial-value problem using @ivp
prob = @ivp(x' = truck_trailer_system, x(0) ∈ X0)

# Time span for reachability analysis
T = 10.0

# Perform reachability analysis
sol = solve(prob, tspan=(0.0, T), alg=TMJets21a())

# Plotting the results
using Plots

# Plot showing θ1 over time
fig1 = plot(sol; vars=(0, 1), xlab="t (s)", ylab="theta1 (rad)", lw=0.5, color=:blue)
display(fig1)

# Plot showing x1 (position) over time
fig2 = plot(sol; vars=(0, 2), xlab="t (s)", ylab="x1 (m)", lw=0.5, color=:green)
display(fig2)

# Plot showing y1 (position) over time
fig3 = plot(sol; vars=(0, 3), xlab="t (s)", ylab="y1 (m)", lw=0.5, color=:red)
display(fig3)

# Phase plot of x1 vs. y1
fig4 = plot(sol; vars=(2, 3), xlab="x1 (m)", ylab="y1 (m)", lw=0.5, color=:purple)
display(fig4)

