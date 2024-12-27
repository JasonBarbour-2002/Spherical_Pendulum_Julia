push!(LOAD_PATH, "src")
using Roots
using GLMakie
using HCubature
using StaticArrays
using LaTeXStrings
using OrdinaryDiffEq
using pendulum: spherical_pendulum!, to_cartesian, phi_dot, Energy, Energy0, stable_point, estimation_theta_thetadot, Period_theta, Period_phi, compute_Tphi
using My_Themes: custum_dark_Theme, custum_light_Theme
# Set theme to light or dark mode: 
Theme_color = "dark" 
# Theme_color = "light"
GLMakie.closeall() # Close any pervious figures

# Change here to change the theme
dark = custum_dark_Theme()
light = custum_light_Theme()
my_theme  = (Theme_color == "dark") ? dark : light
if my_theme == dark
    rod_color = "#FFFFFF"
    ball_color = "#589AFA" 
    trajectory_color = "#FF9300"
    center_color = "#CE3E00"
    Label_color = :white
else
    rod_color = :black
    ball_color = :blue
    trajectory_color = :green
    center_color = :red
    Label_color = :black
end
set_theme!(my_theme)

l = 3.0
g = 9.81
m = 1.0
p = [g, l, m]
Status = Observable("Start")
is_Energy = Observable("")

# Initialize the plot
# Almost all of the bellow is for the plot till line 140
fig = Figure(size=(1600,1000));
ax = Axis3(fig[1:4,1:4], xlabel = L"x", ylabel = L"y", zlabel = L"z", title = "Spherical Pendulum")
params_slider = (
    format= x -> string(round(x, digits = 1)),
    startvalue=0.0,
)
fig[5, 1:6] = slides = GridLayout()
sg = SliderGrid(
    slides[1, 1:2],
    (label=L"L_z", range=0:0.1:20, params_slider...),
    (label="Energy", range=0:0.1:300, params_slider...),
)
fig[5,7] = bg = GridLayout()
start, reset = bg[2:3,1] = [
    Button(fig; label=Status, width=100);
    Button(fig; label="Reset", width=100)
    ]
Label(bg[1,1], is_Energy, color = :red, width = 200, height=50, justification=:center)

Inputs = [s.value for s in sg.sliders]

Lz, E = Inputs 

# Find a corresponding θ and θ̇
# Let's try to make a reasonable guess as to what is the θ and θ̇
estimates = @lift estimation_theta_thetadot($Lz, $E, p)
θ_in , θ̇_in, E_diff = [@lift $estimates[i] for i in 1:3]

θ = @lift copy($θ_in)
θ̇ = @lift copy($θ̇_in)
ϕ = Observable(0.0)
ϕ̇ = @lift phi_dot($θ_in, $Lz, p)

Period = @lift Period_theta($Lz, $E, p)
Period_θ, θ₀, _θ̇ = [@lift $Period[i] for i in 1:3]
Period_ϕ = @lift compute_Tphi($Lz, $E, p)
Label_Period_θ = @lift latexstring("T_\\theta = " * string(round($Period_θ, digits = 5)))
Label_Period_ϕ = @lift latexstring("T_\\phi = " * string(round($Period_ϕ, digits = 5)))
Label(slides[2, 1], Label_Period_θ, color = Label_color, width = 200, height = 40, justification=:center, tellwidth = false)
Label(slides[2, 2], Label_Period_ϕ, color = Label_color, width = 200, height = 40, justification=:center, tellwidth = false)

cartesianCoordinates = @lift to_cartesian(l, $θ, $ϕ)
ball = @lift Point3f($cartesianCoordinates)
rod = @lift Point3f[(0, 0, 0), $cartesianCoordinates]
trajectory = Observable(Point3f[])

# Plot origine
scatter!(ax, 0, 0, 0, markersize = 10, color = center_color)
# Plot the rod
lines!(ax, rod, color = rod_color)
# Plot the ball
scatter!(ax, ball, markersize = 20, color = ball_color)
# Plot the trajectory
lines!(ax, trajectory, color = trajectory_color, linestyle = :dash)

xlims!(ax, -l, l)
ylims!(ax, -l, l)
zlims!(ax, -l, l)

s = 100
theta = range(-π, π, length=s)
theta_dot = range(-3π, 3π, length=s)
theta_mesh = theta * ones(s)'
theta_dot_mesh = ones(s) * theta_dot'
Energy_mesh = @lift Energy(theta_mesh, theta_dot_mesh, $Lz, p)

ax2 = Axis(fig[1:2,5:6]; xlabel = L"\theta", ylabel = L"\dot\theta", title = "Energy Contour Plot", aspect = AxisAspect(1) )
# Plot Contour plot
# Colors
con = heatmap!(ax2, theta, theta_dot, Energy_mesh, colorrange=(0, 300), interpolate = true)
# Constour lines
contour!(ax2, theta, theta_dot, Energy_mesh, levels = range(0, 300, length=30), color=:black, transparency = true)
# Ball on the plot
pen = scatter!(ax2, θ, θ̇, markersize = 20, color=ball_color)
translate!(pen, 0, 0, 1)
Colorbar(fig[1:2,7], con, label = "Energy")
xlims!(ax2, -π, π)
ylims!(ax2, -3π, 3π)

ax3 = Axis(fig[3,5:7]; xlabel = "Energy", ylabel = L"T_\theta", title = L"$T_\theta$ vs Energy Plot")
E_range = range(0, 300, length=100)
thetaPeriod = @lift [Period_theta($Lz, En, p) for En in E_range]
θ_Period = @lift first.($thetaPeriod)
thetaPeriod_now = @lift Period_theta($Lz, $E, p)[1]
scatter!(ax3, E_range, θ_Period, color = :blue)
scatter!(ax3, E, thetaPeriod_now, color = :red)
xlims!(ax3, 0, 300)
ylims!(ax3, 0, 10)

ax4 = Axis(fig[4,5:7]; xlabel = L"L_z", ylabel = L"T_\phi", title = L"$T_\phi$ vs Angular Momentum in Z Plot")
Lz_range = collect(0:0.1:20)
phiPeriod = @lift [compute_Tphi(Lz, $E, p) for Lz in Lz_range]
phiPeriod_now = @lift compute_Tphi($Lz, $E, p)
scatter!(ax4, Lz_range, phiPeriod, color = :blue)
scatter!(ax4, Lz, phiPeriod_now, color = :red)
xlims!(ax4, 0, 20)
ylims!(ax4, 0, 10)
display(fig) # Show the figure

function reset_plot()
    # Reset the plot back to it's original positions
    trajectory[] = Point3f[]
    θ[] = θ_in[]
    θ̇[] = θ̇_in[]
    ϕ[] = 0
    ϕ̇[] = phi_dot(θ_in[], Lz[], p)
end

function anim_step!(integ)
    # Evolve the system to the next step to animate
    for _ in 1:5
        step!(integ)
    end
    # Send everything to the plot
    t = integ[1]
    if t > 0 
        θ[] = ((abs(t % 2π) > π) ? ((t % 2π) - 2π) : t % π)
    else
        θ[] = ((abs(t % 2π) > π) ? ((t % 2π) + 2π) : t % π)
    end
    ϕ[] = integ[2]
    θ̇[] = integ[3]
    ϕ̇[] = integ[4]
    coo = to_cartesian(l, integ[1], integ[2])
    new_point = Point3f(coo)
    trajectory[] = push!(trajectory[], new_point)
end

function animate_pendulum(integ, is_running)
    # This function takes care of animating the pendulum
    @async while is_running[]
        if !isopen(fig.scene)
            is_running[] = false;
            break
        end
        anim_step!(integ)
        sleep(0.005)
    end
end
# Here on out it's all about buttons and sliders you can ignore this part
is_running = Observable(false)
y₀ = [θ[], ϕ[], θ̇[], ϕ̇[]]
total_time = s
# Initialize the ODE Problem 
problem = ODEProblem(spherical_pendulum!, y₀, (0.0, total_time), p)
integ = Observable{Any}(init(problem, Tsit5(); adaptive = false, dt = 0.001))
on(start.clicks) do n
    if abs(E_diff[]) > 1
        is_Energy[] = "Energy too small"
    else 
        is_Energy[] = ""
        is_running[] = !is_running[]
    end
end
on(start.clicks) do n
    if is_Energy[] == ""
        if Status[] == "Start"
            y₀ = [θ[], ϕ[], θ̇[], ϕ̇[]]
            total_time = 1000
            # Initialize the ODE Problem 
            problem = ODEProblem(spherical_pendulum!, y₀, (0.0, total_time), p)
            integ[] = init(problem, Tsit5(); adaptive = false, dt = 0.001)
            Status[] = "Pause"
        elseif Status[] == "Pause"
        Status[] = "Resume"
        elseif Status[] == "Resume"
            Status[] = "Pause"
        end 
        animate_pendulum(integ[], is_running)
    end
end

on(reset.clicks) do n; is_running[] = false; end
on(reset.clicks) do n
    Status[] = "Start"
    reset_plot()
end;