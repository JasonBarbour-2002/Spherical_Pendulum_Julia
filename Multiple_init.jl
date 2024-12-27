push!(LOAD_PATH, "src")
import Random
using GLMakie
using OrdinaryDiffEq
using StaticArrays
using pendulum: spherical_pendulum!, to_cartesian, phi_dot, Energy, Energy0, stable_point, estimation_theta_thetadot
using My_Themes: custum_dark_Theme, custum_light_Theme


Theme_color = "dark" 
# Theme_color = "light" 
GLMakie.closeall() # Close any pervious figures

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

function create_ball(x, y, radius, Lz, p, num_points)
    Random.seed!(42)
    θs = zeros(num_points)
    θ̇s = zeros(num_points)
    Es = zeros(num_points)
    for i in 1:num_points
        t = 2π * rand()
        r = rand() + rand() # to get a uniform distribution in the unit circle
        r = r > 1 ? 2 - r : r
        r *= radius
        θs[i] = x + r*cos(t)
        θ̇s[i] = y + r*sin(t)
        Es[i] =Energy(θs[i], θ̇s[i], Lz, p)
    end
    return θs, θ̇s, Es
end

# Initial conditions
l = 3.0
g = 9.81
m = 1.0
p = SVector{3, Float64}(g, l, m)
num_pendulum = 100

fig = Figure(size = (1000, 1600))
params_slider = (
    format= x -> string(round(x, digits = 1)),
)
fig[2, 1] = slides = GridLayout()
sg = SliderGrid(
    slides[1, 1:2],
    (label="Energy", range=0:0.1:250, startvalue = 10, params_slider...),
    (label="Radius", range=0.01:0.01:2, startvalue = 0.1, params_slider...),
    (label="Lz", range=0:0.1:20, startvalue = 0.0, params_slider...),
)
Status = Observable("Start")
fig[2,2] = bg = GridLayout()
start, reset = bg[2:3,1] = [
    Button(fig; label=Status, width=100);
    Button(fig; label="Reset", width=100)
    ]
Inputs = [s.value for s in sg.sliders]

init_E = Inputs[1]
radius = Inputs[2]
init_Lz = Inputs[3]

estimates = @lift estimation_theta_thetadot($init_Lz, $init_E, p)
θ_in, θ̇_in, E_diff = [@lift $estimates[i] for i in 1:3]
# ball of initial conditions centered at θ_in, θ̇_in
ball = @lift create_ball($θ_in, $θ̇_in, $radius, $init_Lz, p, num_pendulum)
θ_list, θ̇_list, E_list = [@lift $ball[i] for i in 1:3]



θ = @lift copy($θ_list)
θ̇ = @lift copy($θ̇_list)
ϕ = Observable(zero(θ_list[]))
ϕ̇ = @lift [phi_dot(θ_, $init_Lz, p) for θ_ in θ_list[]]

s = 100
theta = range(-π, π, length=s)
theta_dot = range(-3π, 3π, length=s)
theta_mesh = theta * ones(s)'
theta_dot_mesh = ones(s) * theta_dot'
Energy_mesh = @lift Energy(theta_mesh, theta_dot_mesh, $init_Lz, p)

ax = Axis(fig[1, 1]; xlabel = L"\theta", ylabel = L"\dot{\theta}", aspect=AxisAspect(1))

con = heatmap!(ax, theta, theta_dot, Energy_mesh, colorrange=(0, 300), interpolate = true)
contour!(ax, theta, theta_dot, Energy_mesh, levels = range(0, 300, length=30), color=:black, transparency = true)

pen = scatter!(ax, θ, θ̇, markersize = 5, color = :blue)
translate!(pen, 0, 0, 1)
Colorbar(fig[1, 2], con, label = "Energy")
xlims!(ax, -π, π)
ylims!(ax, -3π, 3π)
display(fig)

function reset_plot()
    θ[] = θ_list[]
    θ̇[] = θ̇_list[]
    ϕ[] = zeros(num_pendulum)
    ϕ̇[] = [phi_dot(θ_, init_Lz[], p) for θ_ in θ_list[]]
end

function anim_step!(integ_list, num_pendulum = num_pendulum)
    θ_l = zeros(num_pendulum)
    ϕ_l = zeros(num_pendulum)
    θ̇_l = zeros(num_pendulum)
    ϕ̇_l = zeros(num_pendulum)
    for (i, integ) ∈ enumerate(integ_list)
        for _ in 1:5
            step!(integ)
        end
        t = integ[1]
        if t > 0 
            θ_l[i] = ((abs(t % 2π) > π) ? ((t % 2π) - 2π) : t % π)
        else
            θ_l[i] = ((abs(t % 2π) > π) ? ((t % 2π) + 2π) : t % π)
        end 
        ϕ_l[i] = integ[2]
        θ̇_l[i] = integ[3]
        ϕ̇_l[i] = integ[4]
    end
    θ[] = θ_l
    ϕ[] = ϕ_l
    θ̇[] = θ̇_l
    ϕ̇[] = ϕ̇_l
end

function animate_pendulum(integ_list, is_running, num_pendulum = num_pendulum)
    @async while is_running[]
        if !isopen(fig.scene)
            is_running[] = false;
            break
        end
        anim_step!(integ_list, num_pendulum)
        sleep(0.01)
    end
end

is_running = Observable(true)
problems = [ODEProblem(spherical_pendulum!, [θ[][i], ϕ[][i], θ̇[][i], ϕ̇[][i]], (0.0, 100.0), p) for i in 1:num_pendulum]
integ_list = [init(problems[i], Tsit5(); adaptive=false, dt=0.001) for i in 1:num_pendulum]
on(start.clicks) do n
    if Status[] == "Start"
        Status[] = "Pause"
        is_running[] = true
        problems[:] = [ODEProblem(spherical_pendulum!, [θ[][i], ϕ[][i], θ̇[][i], ϕ̇[][i]], (0.0, 100.0), p) for i in 1:num_pendulum]
        integ_list[:] = [init(problems[i], Tsit5(); adaptive=false, dt=0.001) for i in 1:num_pendulum]
        animate_pendulum(integ_list, is_running)
    elseif Status[] == "Pause"
        Status[] = "Resume"
        is_running[] = false
    else 
        Status[] = "Pause"
        is_running[] = true
        animate_pendulum(integ_list, is_running)
    end

end

on(reset.clicks) do n
    Status[] = "Start"
    is_running[] = false
    problems[:] = [ODEProblem(spherical_pendulum!, [θ[][i], ϕ[][i], θ̇[][i], ϕ̇[][i]], (0.0, 100.0), p) for i in 1:num_pendulum]
    integ_list[:] = [init(problems[i], Tsit5(); adaptive=false, dt=0.001) for i in 1:num_pendulum]
    reset_plot()
end
