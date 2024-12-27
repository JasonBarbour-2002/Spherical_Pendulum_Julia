push!(LOAD_PATH, "src")
using GLMakie
using StaticArrays
using OrdinaryDiffEq
using pendulum: spherical_pendulum!, to_cartesian, Energy, phi_dot
using My_Themes: custum_dark_Theme, custum_light_Theme

# Set theme to light or dark mode: 
Theme_color = "dark" 
# Theme_color = "light"


GLMakie.closeall() # Close any pervious figures

function Angular_Momentum_z(θ, ϕ̇, p)
    g, l, m = p
    return m * l^2 * ϕ̇ * sin(θ)^2
end

dark = custum_dark_Theme()
light = custum_light_Theme()
my_theme  = (Theme_color == "dark") ? dark : light
# Change here to change the theme
with_theme(my_theme) do
    if my_theme == dark
        rod_color = "#FFFFFF"
        ball_color = "#589AFA" 
        trajectory_color = "#FF9300"
        center_color = "#CE3E00"
    else
        rod_color = :black
        ball_color = :blue
        trajectory_color = :green
        center_color = :red
    end 
    l = 3.0
    g = 9.91
    m = 1.0
    p = SVector{3, Float64}(g, l, m)
    Status = Observable("Start")

    fig = Figure(size=(1600,1000));
    ax = Axis3(fig[1,1:3])
    ax2 = Axis(fig[1,4:6]; xlabel = L"\theta", ylabel = L"$\dot\theta$", title = "Phase space")

    params_slider = (
        format= x -> string(round(x/π, digits = 2)) * "π",
        startvalue=0.0,
    )
    sg = SliderGrid(
        fig[2,1:6],
        (label=L"\theta", range=-π:0.01:π, params_slider...),
        (label=L"\dot\theta", range=-2π:0.01:2π, params_slider...),
        (label=L"\phi", range=0:0.01:π, params_slider...),
        (label=L"\dot\phi", range=0:0.01:2π, params_slider...)
    )

    fig[2,7] = bg = GridLayout()
    start, reset = bg[1:2,1] = [Button(fig; label=Status, width=100); Button(fig; label="Reset", width=100)]
    spherical = [s.value for s in sg.sliders]
    θ, θ̇, ϕ, ϕ̇ = spherical
    cart = @lift to_cartesian(l, $θ, $ϕ)
    ball = @lift Point3f($cart)
    rod = @lift Point3f[(0, 0, 0), $cart]
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
    # Phase space
    Lz = @lift Angular_Momentum_z($θ, $ϕ̇, p)
    s = 100
    theta = range(-π, π, length=s)
    theta_dot = range(-3π, 3π, length=s)
    theta_mesh = theta * ones(s)'
    theta_dot_mesh = ones(s) * theta_dot'
    E = @lift Energy(theta_mesh, theta_dot_mesh, $Lz, p)
    teta = @lift copy($θ)
    tetadot = @lift copy($θ̇)

    con = heatmap!(ax2, theta, theta_dot, E, colorrange=(0, 300), interpolate = true)
    contour!(ax2, theta, theta_dot, E, levels = range(0, 300, length=30), color=:black, transparency = true)
    pen = scatter!(ax2, teta, tetadot, markersize = 20, color=ball_color)
    translate!(pen, 0, 0, 10)
    Colorbar(fig[1,7], con, label = "Energy")
    xlims!(ax2, -π, π)
    ylims!(ax2, -3π, 3π)

    display(fig)

    function reset_plot()
        cart[] = to_cartesian(l, θ[], ϕ[])
        trajectory[] = Point3f[]
        teta[] = θ[]
        tetadot[] = θ̇[]
    end

    function anim_step!(integ)
        for i in 1:5
            step!(integ)
        end
        t = integ[1]
        if t > 0 
            teta[] = ((abs(t % 2π) > π) ? ((t % 2π) - 2π) : t % π)
        else
            teta[] = ((abs(t % 2π) > π) ? ((t % 2π) + 2π) : t % π)
        end
        tetadot[] = integ[3]
        coo = to_cartesian(l, integ[1], integ[2])
        cart[] = coo
        new_point = Point3f(coo)
        trajectory[] = push!(trajectory[], new_point)
    end

    function animate_pendulum(integ, is_running)
        @async while is_running[]
            if !isopen(fig.scene)
                is_running[] = false;
                break
            end
            anim_step!(integ)
            sleep(0.005)
        end
    end

    is_running = Observable(false)

    θ, θ̇, ϕ, ϕ̇ = spherical
    y₀ = [θ[], ϕ[], θ̇[], ϕ̇[]]
    total_time = s
    # Initialize the ODE Problem 
    problem = ODEProblem(spherical_pendulum!, y₀, (0.0, total_time), p)
    integ = Observable{Any}(init(problem, Tsit5(); adaptive = false, dt = 0.001))

    on(start.clicks) do n; is_running[] = !is_running[]; end
    on(start.clicks) do n
        if Status[] == "Start"
            θ, θ̇, ϕ, ϕ̇ = spherical
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

    on(reset.clicks) do n; is_running[] = false; end
    on(reset.clicks) do n
        Status[] = "Start"
        reset_plot()
    end
end;