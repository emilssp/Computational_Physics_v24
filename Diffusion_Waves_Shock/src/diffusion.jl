module Diffusion

include("utils/diffusion_funcs.jl")
using LaTeXStrings
    const diffusivity = 1.00::Float64
    
    const endtime = 1.00::Float64
    const L = 1.00::Float64

    const dt = 1e-4::Float64
    const timesteps = Int64(div(endtime,dt)+1)::Int64

    const spacesteps = 100::Int64
    const dx = L/spacesteps
    const u0_tilda = 1::Int64
    const u0 = (u0_tilda*spacesteps/L)::Float64

    const cfl = (diffusivity*dt/dx^2)::Float64

    const D_left = 0.50::Float64
    const D_right = 1.00::Float64
    
    const subdir = "figures/diffusion/"

    ##########################################################################

    function plot_array(vals_CN, vals_an, time_arr, x, t; figname = "const_D/new_fig.png")
        plt = plot(xlabel="x(a.u.)", ylabel = "u(x,t)" )
        for time in time_arr
            scatter!(x,vals_CN[round(Int, time/dt),:], label = round(t[round(Int, time/dt)],digits=4))
            plot!(x,vals_an[round(Int, time/dt),:], color = "black", label = "")#t[round(Int, time/dt)])
        end
        scatter!(x,vals_CN[end,:], label = round(t[end],digits = 4))
        plot!(x,vals_an[end,:], color = "black", label = "")#t[round(Int, time/dt)])
        savefig(plt, joinpath(subdir,figname))
        return plt 
    end

    function main()

        x = LinRange(0,L, spacesteps)
        t = LinRange(0,endtime, timesteps)

        u_init = zeros(spacesteps)
        u_init[div(spacesteps, 2)] = u0/2
        u_init[div(spacesteps, 2)+1] = u0/2


        u_res_absorbing = CN_absorbing_BC(u_init, cfl)

        u_res_reflecting = CN_reflecting_BC(u_init, cfl)

        u_an_unbounded = unbounded.(L/2, x', t[2:end])

        u_an_absorbing = absorbing_bounds.(L/2, x', t)
        u_an_reflecting= reflecting_bounds.(L/2, x', t)
        print(size(u_an_reflecting[10,:]))


        println("Done")

        plot_array(u_res_absorbing, u_an_unbounded, [0.01,0.04,0.07], x, t; figname = "const_D/u_res_absorbing.png")
        plot_array(u_res_reflecting, u_an_unbounded, [0.01,0.04,0.07], x, t; figname = "const_D/u_res_reflecting.png")
        plot_array(u_res_absorbing, u_an_absorbing, [0.01,0.04,0.07], x, t; figname = "const_D/u_an_absorbing.png")
        plot_array(u_res_reflecting, u_an_reflecting, [0.01,0.04,0.07], x, t; figname = "const_D/u_an_reflecting.png")

        x_boundary = x[div(spacesteps, 3)]
        step_alphas = step_diffusivity.(x,x_boundary)

        absorbing_variable_D = absorbing_step(u_init, step_alphas)
        reflecting_variable_D = reflecting_step(u_init, step_alphas)
        
        plt2 = plot()
        scatter!(x,absorbing_variable_D[round(Int, 0.07/dt)])
        scatter!(x,reflecting_variable_D[round(Int, 0.07/dt)])
        vline!([x_boundary],color ="black", label="")
        savefig(plt2, joinpath(subdir, "variable_D/step_D.png"))

        alphas = variable_diffusivity(linear_func, x)

        absorbing_variable_D = absorbing_step(u_init, alphas)
        reflecting_variable_D = reflecting_step(u_init, alphas)

        plt3 = plot()
        scatter!(x,absorbing_variable_D[round(Int, 0.05/dt)])
        scatter!(x,reflecting_variable_D[round(Int, 0.05/dt)])
        plot!(x,linear_func.(x))
        savefig(plt3, joinpath(subdir,"variable_D/cont_D.png"))
    end
    function test_timestep()

        x = LinRange(0,L, spacesteps)
        
        u_init = zeros(spacesteps)
        u_init[div(spacesteps, 2)] = u0/2
        u_init[div(spacesteps, 2)+1] = u0/2

        dt_array = 0.1*dt:5*dt:2*dt
        difference_abs = zeros(size(dt_array)[1])
        difference_ref = zeros(size(dt_array)[1])
        for i in 1:size(dt_array)[1]
            α = diffusivity*dt_array[i]/dx^2
            t = 0:dt_array[i]:endtime

            u_res_absorbing = CN_absorbing_BC(u_init, α)

            u_res_reflecting = CN_reflecting_BC(u_init, α)
    
            u_an_absorbing = absorbing_bounds.(L/2, x', t)

            u_an_reflecting= reflecting_bounds.(L/2, x', t)


            difference_abs[i] = maximum(abs.(u_res_absorbing[round(Int64, 0.01/dt_array[i])] .- u_an_absorbing[round(Int64, 0.01/dt_array[i]),:]))
            difference_ref[i] = maximum(abs.(u_res_reflecting[round(Int64, 0.01/dt_array[i])] .- u_an_reflecting[round(Int64, 0.01/dt_array[i]),:]))
        end

        plt = plot( title = "Max error with absorbing BC at time t = 0.01", xlabel = "dt", ylabel = "Error")
        plot!(dt_array, difference_abs, label ="Absorbing boundary")
        plot!(dt_array, difference_ref, label ="Reflective boundary")
        savefig(plt, joinpath(subdir,"test_dt_1.png"))


    end
end

# Diffusion.main()

Diffusion.test_timestep()