module Hopf
    
    using Plots
    include("utils/hopfs_eqn_funcs.jl")

    const spacesteps = 1000::Int64
    const L = 1.00::Float64
    const h = (L/spacesteps)::Float64

    const c = 1.00::Float64

    const endtime = 0.10::Float64
    const dt = abs(h/c)::Float64#1e-4::Float64
    const timesteps = Int64(div(endtime,dt)+1)::Int64

    const γ = c*dt/h

    #################################################################################
    function transport_eqn()
        x = LinRange(0, 2*L, spacesteps)
        t = LinRange(0, endtime, timesteps)
        u_init = u0.(x)

        u = Lax_Wandroff_trans(u_init, spacesteps, timesteps; γ)

        u_an = trans_anal(u0, x, t; c)

        anim = @animate for n in 1:10:timesteps
            p1 = plot(x, u[:, n], ylim =(-0.1, 1.1))
            p2 = plot!(x, u_an[:, n], ylim =(-0.1, 1.1))
        end
        gif(anim, "figures/transport_eqn/advection_equation_solution_new.gif", fps=7)
        
    end
    
    ##################################################################################

    function hopfs_eqn()
        x = LinRange(0, L, spacesteps)
        t = LinRange(0, endtime, timesteps)
        u_init = sin.(π * x)

        u = Lax_Wandroff_Hopf(u_init, spacesteps, timesteps; γ)
        anim = @animate for n in 1:timesteps
            p1 = plot(x, u[:, n])
            # p = plot(p1,p2, layout= (1,2), size = (800,500), suptitle = "\nAnalitical and numerical solution at time t = $(round(time[n], digits = 4))")
        end
        gif(anim, "figures/transport_eqn/hopf_equation_solution.gif", fps=7)
    end

    ##################################################################################
end

Hopf.transport_eqn()
Hopf.hopfs_eqn()