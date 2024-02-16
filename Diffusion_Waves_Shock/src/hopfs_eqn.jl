module Hopf
    
    using Plots
    include("./utils/hopfs_eqn_funcs.jl")

    const spacesteps = 200::Int64
    const L = 2.00::Float64
    const h = (L/spacesteps)::Float64

    const c = 2.00::Float64

    const endtime = 1.00::Float64

    #################################################################################
    function transport_eqn()
        dt = abs(h/(2*c))::Float64#1e-4::Float64
        timesteps = Int64(div(endtime,dt)+1)::Int64
        γ = c*dt/h
        
        x = LinRange(0, L, spacesteps)
        t = LinRange(0, endtime, timesteps)
        u_init = step_func.(x)

        u = Lax_Wandroff_trans(u_init, spacesteps, timesteps; γ)
        u_an = trans_anal(step_func, x, t; c)
        
        difference = zeros(timesteps)
        
        for i in length(timesteps)
            difference[i] = maximum(abs.(u[i] .- u_an[i,:]))
        end
        print(maximum(difference))
        p1 = plot(xlabel = "\$x\$ (a.u.)", ylabel = "\$u(x,t)\$", ylims = (-0.1,1.4), title = "\nCFL = $(round(γ, digits = 3))")
        plot!(x, u[:, 20], label = "t = $(round(t[10],digits = 3))")
        plot!(x, u[:, 100], label = "t = $(round(t[100],digits = 3))")
        plot!(x, u[:, 200], label = "t = $(round(t[200],digits = 3))")
        savefig(p1, "./figures/transport_eqn/advection_equation_CFL05.png")

        dt = abs(h/(c))::Float64#1e-4::Float64
        timesteps = Int64(div(endtime,dt)+1)::Int64
        t = LinRange(0, endtime, timesteps)
        γ = c*dt/h
        u = Lax_Wandroff_trans(u_init, spacesteps, timesteps; γ)
        
        p2 = plot(xlabel = "\$x\$ (a.u.)", ylabel = "\$u(x,t)\$", ylims = (-0.1,1.4), title = "\nCFL = $(round(γ, digits = 3))")
        plot!(x, u[:, 10], label = "t = $(round(t[10],digits = 3))")
        plot!(x, u[:, 55], label = "t = $(round(t[50],digits = 3))")
        plot!(x, u[:, 99], label = "t = $(round(t[100],digits = 3))")
        savefig(p2, "./figures/transport_eqn/advection_equation_CFL1.png")

        savefig(plot(p1,p2,layout = (1,2), suptitle = "Transport equation for different CFL number"), "./figures/transport_eqn/advection_subplots.png")
    end
    
    ##################################################################################
    function hopf_init(x)
        return sin.(x)
    end

    function hopfs_eqn()
        dt = abs(h/(c))::Float64#1e-4::Float64
        timesteps = Int64(div(endtime,dt)+1)::Int64
        γ = c*dt/h
        x = LinRange(0, L, spacesteps)
        t = LinRange(0, endtime, timesteps)
        
        # u_init = 2*sin.(π*x)
        # u_init = step_func.(x)
        u_init = gaussian.(x)

        u = Lax_Wandroff_Hopf(u_init, spacesteps, timesteps; γ)
        # u_an = 2*sin.(π*(x .- u_init .* t'))

        # u_an = step_func.(x .- u_init .* t')
        u_an = gaussian.(x .- u_init .* t')

        anim = @animate for n in 1:timesteps
            plot(x, u[:, n],ylims=(0,2.5))
        end
        gif(anim, "./figures/transport_eqn/hopf_equation_solution.gif", fps=60)
        
        anim = @animate for n in 1:timesteps
            plot(x, u_an[:,n])
        end
        gif(anim, "./figures/transport_eqn/hopf_char_solution.gif", fps=60)
        
        p2 = plot(xlabel = "\$x\$ (a.u.)", ylabel = "\$u(x,t)\$", ylims = (-0.1,1.4), xlims = (0,1.2), title = "Hopf's equation and shock")
        plot!(x, u[:, 10], label = "t = $(round(t[10],digits = 3))")
        plot!(x, u[:, 55], label = "t = $(round(t[50],digits = 3))")
        plot!(x, u[:, 99], label = "t = $(round(t[100],digits = 3))")
        savefig(p2, "./figures/transport_eqn/hopf_CFL1.png")
    
    end
    ##################################################################################
end


Hopf.transport_eqn()
# Hopf.hopfs_eqn()        

