
export efr_timestep!
import IncompressibleNavierStokes.ExplicitRungeKuttaMethod
import IncompressibleNavierStokes:apply_bc_u!
import IncompressibleNavierStokes:apply_bc_temp!
import IncompressibleNavierStokes:momentum!
import IncompressibleNavierStokes:project!

create_stepper(::ExplicitRungeKuttaMethod; setup, psolver, u, temp, t, n = 0) =
    (; setup, psolver, u, temp, t, n)

function efr_timestep!(
    method::ExplicitRungeKuttaMethod,
    stepper,
    Δt;
    θ = nothing,
    cache,
    filter = nothing,
    relax = nothing,
    project_after_filter::Bool = false,
)
    (; setup, psolver, u, temp, t, n) = stepper
    (; closure_model, temperature) = setup
    (; A, b, c) = method
    # Read pre-computed cache in case of differential filter
    (; ustart, ku, p, tempstart, ktemp, diff) = cache
    lu_filter_mats = hasproperty(cache, :lu_filter_mats) ? cache.lu_filter_mats : nothing

    nstage = length(b)
    m = closure_model

    # Save starting time
    tstart = t
    copyto!(ustart, u)
    isnothing(temp) || copyto!(tempstart, temp)

    for i = 1:nstage
        apply_bc_u!(u, t, setup)
        isnothing(temp) || apply_bc_temp!(temp, t, setup)
        momentum!(ku[i], u, temp, t, setup)

        if !isnothing(temp)
            ktemp[i] .= 0
            convection_diffusion_temp!(ktemp[i], u, temp, setup)
            temperature.dodissipation && dissipation!(ktemp[i], diff, u, setup)
        end

        # Add closure term
        isnothing(m) || (ku[i] .+= m(u, θ, closure_stuff, setup))

        # Stage time
        t = tstart + c[i] * Δt

        # Combine stage results
        u .= ustart
        for j = 1:i
            @. u += Δt * A[i, j] * ku[j]
        end
        if !isnothing(temp)
            temp .= tempstart
            for j = 1:i
                @. temp += Δt * A[i, j] * ktemp[j]
            end
        end

        apply_bc_u!(u, t, setup)
        project!(u, setup; psolver, p)
    end

    # === EFR logic ===
    if !isnothing(filter)
        if !isnothing(filter.filter_radius)
            u_filtered = apply_filter(filter, stepper.u, stepper.setup, lu_filter_mats)
        else
            u_filtered = apply_filter(filter, stepper.u, stepper.setup)
        end
        
        
        if !isnothing(relax)
            u_relaxed = apply_relax(relax, stepper.u, u_filtered)
            copyto!(stepper.u, u_relaxed)
        else
            copyto!(stepper.u, u_filtered)
        end
        
    end

    if project_after_filter
        apply_bc_u!(stepper.u, stepper.t, stepper.setup)
        project!(stepper.u, stepper.setup; psolver = stepper.psolver, p = cache.p)
    end
    # Final BCs
    apply_bc_u!(u, t, setup)
    isnothing(temp) || apply_bc_temp!(temp, t, setup)
    # copy of the full timestep! logic, then at the end:
    create_stepper(method; setup, psolver, u, temp, t, n = n + 1)

end

