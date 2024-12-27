module pendulum
    using Roots
    using OrdinaryDiffEq
    using HCubature
    using LinearAlgebra
    using StaticArrays

    function spherical_pendulum!(du, u, p, t)
        # The equations of motion for the spherical pendulum
        θ, ϕ, θ̇, ϕ̇ = u
        g, l, m = p
        du[1] = θ̇
        du[2] = ϕ̇
        du[3] = ϕ̇^2 * sin(θ) * cos(θ) - g/l * sin(θ)
        du[4] = -2 * θ̇ * ϕ̇ * cot(θ)
        return nothing
    end

    function to_cartesian(r::Float64, θ::Float64, ϕ::Float64)
        # Convert spherical coordinates to cartesian coordinates
        X = r * sin(θ) * cos(ϕ)
        Y = r * sin(θ) * sin(ϕ)
        Z = r * cos(θ)
        return X, Y, -Z
    end
    function to_cartesian(r, θ, ϕ)
        # Convert spherical coordinates to cartesian coordinates
        X = r .* sin.(θ) .* cos.(ϕ)
        Y = r .* sin.(θ) .* sin.(ϕ)
        Z = r .* cos.(θ)
        return X, Y, - Z
    end

    function phi_dot(θ::Float64, Lz::Float64, p)
        # Compute the angular velocity in the ϕ direction
        g, l, m = p
        return Lz / (m * l^2 * sin(θ)^2)
    end

    function Energy(θ::Matrix{T}, θ̇::Matrix{T}, Lz::Float64, p) where T
        # Compute the energy of the system
        g, l, m = p
        Term1 = 0.5m .* l^2 .* θ̇.^2
        Term2 = Lz^2 ./ (2m * l^2 * sin.(θ).^2)
        Term3 = m * g * l * (1 .- cos.(θ))
        return Term1 + Term2 + Term3
    end
    function Energy(θ::Float64, θ̇::Float64, Lz::Float64, p)
        # Compute the energy of the system
        g, l, m = p
        Term1 = 0.5m * l^2 * θ̇^2
        Term2 = Lz^2 / (2m * l^2 * sin(θ)^2)
        Term3 = m * g * l * (1 - cos(θ))
        return Term1 + Term2 + Term3
    end
    function Energy(θ::T, θ̇::U, Lz::V, p) where {T<:Real, U<:Real, V<:Real}
        # Compute the energy of the system with the correct types
        return Energy(convert(Float64, θ), convert(Float64, θ̇), convert(Float64, Lz), p)
    end

    function Energy0(E::Float64, θ, θ̇, Lz, p)
        # Compute the difference between the energy of the system and the given energy
        return Energy(θ, θ̇, Lz, p) - E
    end
    function Energy0(E::T, θ, θ̇, Lz, p) where T<:Real
        # Compute the difference between the energy of the system and the given energy
        return Energy0(convert(Float64, E), θ, θ̇, Lz, p)
    end
    
    function stable_point(θ::Float64, Lz::Float64, p)
        # Derivative of the equations of motion with respect to θ given θ̇ = 0
        g, l, m = p
        return Lz^2 / (m^2 * l^4 * sin(θ)^3) * cos(θ) - g/l * sin(θ) 
    end
    
    function estimation_theta_thetadot(Lz::Float64, E::Float64, p)
        # Given the energy and the angular momentum, compute a possible value for θ and θ̇
        # First check if the Energy is big enough
        min_θ = find_zero(θ -> stable_point(θ, Lz, p), 0.1) + eps()
        min_E = Energy(min_θ, 0.0, Lz, p)
        # If the Energy is too small there is no point in computing this 
        if E < min_E
            return NaN, NaN, Inf64
        end
        g, l, m = p
        θ = min_θ
        PE = Lz^2 / (2m * l^2 * sin(θ)^2) + m * g * l * (1 - cos(θ))
        KE = E - PE
        # KE = m l² θ̇²
        θ̇ = sqrt(2KE / (m*l^2))
        new_E = Energy(θ, θ̇, Lz, p)
        return θ, θ̇, new_E - E
    end
    
    function Period_theta(Lz, E, p)
        # Compute the period of the system in the θ direction
        g, l, m = p
        # The following is just picking the start and end of the path 
        # for different edge cases
        if Lz == 0.0
            θ̇ = 0.0
            if E >= 2 * m * g * l
                # θ̇ ≠ 0
                KE = E - 2 * m * g * l
                θ₀ = -π
                θₜ = π
                θ̇ = sqrt(abs(2KE / (m*l^2)))
            elseif E == 0.0
                # Something is wrong
                return NaN, NaN, NaN
            else 
                # Here θ̇ = 0 at some point. we can find the possible positions in θ
                θₜ = find_zero(θ -> Energy0(E, θ, θ̇, Lz, p), [eps(), π-eps()])
                θ₀ = find_zero(θ -> Energy0(E, θ, θ̇, Lz, p), [-π+eps(), -eps()])
            end
        else
            # Here we have Lz ≠ 0
            # We can find the possible positions in θ
            θ̇ = 0.0
            θmin = find_zero(θ -> stable_point(θ, Lz, p), 0.1) + eps()
            Emin = Energy(θmin, θ̇, Lz, p)
            if E < Emin
                return NaN, NaN, NaN
            end
            θ₀ = find_zero(θ -> Energy0(E, θ, θ̇, Lz, p), (0.0001, θmin))
            θₜ = find_zero(θ -> Energy0(E, θ, θ̇, Lz, p), (θmin, π-0.001))
        end 
        if θ₀ > θₜ
            θ₀, θₜ = θₜ, θ₀
        end
        # Now we can compute the period by integrating the function 1/sqrt(2(E - Ueff(θ, Lz, p)))
        function Ueff(θ, Lz, p)
            g, l, m = p
            if Lz == 0.0
                return m * g * l * (1 - cos(θ))
            end
            return m * g * l * (1 - cos(θ)) + Lz^2 / (2m * l^2 * sin(θ)^2)
        end
        integrand(θ) = 1 / sqrt(abs(2 * (E - Ueff(θ, Lz, p))))
        # Here the integrand should be 1 / sqrt(2 * (E - Ueff(θ, Lz, p))
        # But because of numerical error, the difference between E and Ueff can 
        # be a very small negative number, which will cause problems with the sqrt
        # So we take the absolute value of the difference
        quadrature = hquadrature(integrand, θ₀, θₜ; rtol=1e-4)
        # Here I am multiplying by 2 because this is just the half 
        # trajectory.
        return 2 * l * sqrt(m) * quadrature[1], θ₀, θ̇
    end
    function Period_phi(T,  θ₀, θ̇, Lz, p)
        # Compute the period of the system in the ϕ direction
        # By running the simulation for the θ period, then check how far has
        # the system moved in the ϕ direction and compute the period
        # in ϕ accordingly
        if isnan(θ₀) || isnan(θ̇) || isnan(T) || Lz == 0.0
            return NaN
        end
        ϕ₀ = 0.0
        ϕ̇ = phi_dot(θ₀, Lz, p)
        y₀ = [θ₀, ϕ₀, θ̇, ϕ̇]
        problem = ODEProblem(spherical_pendulum!, y₀, (0.0, T), p)
        sol = solve(problem, Tsit5())
        phi_end = sol[2, end]
        return T * 2π/phi_end
    end
    function compute_Tphi(Lz, E, p)
        # Compute the period of the system in the ϕ direction
        if Lz == 0.0
            return NaN
        end
        T_θ = Period_theta(Lz, E, p)
        if isnan(T_θ[1])
            return NaN
        end
        T_θ, θ₀, θ̇ = T_θ
        T_ϕ = Period_phi(T_θ, θ₀, θ̇, Lz, p)
        return T_ϕ
    end
end # module pendulum
