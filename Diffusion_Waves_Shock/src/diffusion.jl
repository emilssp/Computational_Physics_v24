module Diffusion

include("utils/diffusion_funcs.jl")
using LaTeXStrings
    const diffusivity = 1.00::Float64
    
    const endtime = 1.00::Float64
    const L = 1.00::Float64

    const dt = 1e-4::Float64
    const timesteps = Int64(div(endtime,dt)+1)::Int64

    const spacesteps = 50::Int64
    const dx = L/spacesteps
    const u0_tilda = 1::Int64
    const u0 = (u0_tilda*spacesteps/L)::Float64

    const cfl = (diffusivity*dt/dx^2)::Float64
    
    const subdir = "figures/diffusion/"

    ##########################################################################

    function plot_array(vals_CN, vals_an, time_arr, x, t; title = "Title", figname = "const_D/new_fig.png")
        plt = plot(xlabel="x(a.u.)", ylabel = "u(x,t)", title = title )
        for time in time_arr
            scatter!(x,vals_CN[round(Int, time/dt),:], label = "t = $(round(t[round(Int, time/dt)],digits=4))")
            plot!(x,vals_an[round(Int, time/dt),:], color = "black", label = "")#t[round(Int, time/dt)])
        end
        scatter!(x,vals_CN[end,:], label = "t = $(round(t[end],digits = 4))")
        plot!(x,vals_an[end,:], color = "black", label = "")#t[round(Int, time/dt)])
        savefig(plt, joinpath(subdir,figname))
        return plt 
    end

    function main()
        # print(cfl)
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

        plt_res1=plot_array(u_res_absorbing, u_an_unbounded, [0.01,0.04,0.07], x, t; title = "\nAbsorbing boundary" , figname = "const_D/u_res_absorbing.png")
        plt_res2=plot_array(u_res_reflecting, u_an_unbounded, [0.01,0.04,0.07], x, t; title = "\nReflective boundary" , figname = "const_D/u_res_reflecting.png")
        plt_an1=plot_array(u_res_absorbing, u_an_absorbing, [0.01,0.04,0.07], x, t; title = "\nAbsorbing boundary" , figname = "const_D/u_an_absorbing.png")
        plt_an2=plot_array(u_res_reflecting, u_an_reflecting, [0.01,0.04,0.07], x, t; title = "\nReflective boundary" , figname = "const_D/u_an_reflecting.png")
        savefig(plot(plt_res1,plt_res2, layout = (1,2), suptitle = "Num. result vs. solution for unbounded case"), joinpath(subdir,"const_D/subplts1.png"))
        savefig(plot(plt_an1,plt_an2, layout = (1,2), suptitle = "Numerical result vs. exact solution"), joinpath(subdir,"const_D/subplts2.png"))


        x_boundary = x[div(spacesteps, 3)]
        D_left = 0.99
        D_right = 1.00
        step_alphas = step_diffusivity.(x,x_boundary; D_left, D_right)

        absorbing_variable_D = absorbing_step(u_init, step_alphas)
        reflecting_variable_D = reflecting_step(u_init, step_alphas)
        
        u_an_step = u_an_step_D.(x', t, 0.5, x_boundary; D_left, D_right)

        plt_step1 = plot(xlabel="x(a.u.)", ylabel = "u(x,t)", ylims = (-0.02,2.8), title = "\n \$D_- = 0.99\$, \$D_+ = 1.00\$")
        scatter!(x,absorbing_variable_D[round(Int, 0.02/dt)], label = "Absorbing BC")
        scatter!(x,reflecting_variable_D[round(Int, 0.02/dt)], label = "Reflective BC")
        vline!([x_boundary],color ="black", label="")
        plot!(x, u_an_step[round(Int, 0.02/dt),:], label = "Analytical, unbounded")
        savefig(plt_step1, joinpath(subdir, "variable_D/step_D.png"))

        D_left = 0.50
        step_alphas = step_diffusivity.(x,x_boundary; D_left, D_right)

        absorbing_variable_D = absorbing_step(u_init, step_alphas)
        reflecting_variable_D = reflecting_step(u_init, step_alphas)
        
        u_an_step = u_an_step_D.(x', t, 0.5, x_boundary; D_left, D_right)

        plt_step2 = plot(xlabel="x(a.u.)", ylabel = "u(x,t)", ylims = (-0.02,2.8), title = "\n \$D_- = 0.50\$, \$D_+ = 1.00\$")
        scatter!(x,absorbing_variable_D[round(Int, 0.02/dt)], label = "Absorbing BC")
        scatter!(x,reflecting_variable_D[round(Int, 0.02/dt)], label = "Reflective BC")
        vline!([x_boundary],color ="black", label="")
        plot!(x, u_an_step[round(Int, 0.02/dt),:], label = "Analytical, unbounded")
        savefig(plt_step2, joinpath(subdir, "variable_D/step_D.png"))

        savefig(plot(plt_step1,plt_step2,layout = (1,2), suptitle ="Num. solution vs. unbounded analytical sol."), joinpath(subdir, "variable_D/step_D_subplots.png"))
        
        alphas = variable_diffusivity(linear_func, x)

        absorbing_variable_D = absorbing_step(u_init, alphas)
        reflecting_variable_D = reflecting_step(u_init, alphas)

        plt3 = plot(xlabel="x(a.u.)", ylabel = "u(x,t)", title = "Numerical solution for continuous D(x)" )
        scatter!(x,absorbing_variable_D[round(Int, 0.02/dt)], label = "Absorbing BC")
        scatter!(x,reflecting_variable_D[round(Int, 0.02/dt)], label = "Reflective BC")
        plot!(x,linear_func.(x), label = "Diffusivity D(x)")
        savefig(plt3, joinpath(subdir,"variable_D/cont_D.png"))
    end
    function test_timestep()

        dt = 1e-4
        x = LinRange(0,L, spacesteps)
        
        u_init = zeros(spacesteps)
        u_init[div(spacesteps, 2)] = u0/2
        u_init[div(spacesteps, 2)+1] = u0/2

        α = diffusivity*dt/dx^2
        t = 0:dt:endtime

        u_res_absorbing = CN_absorbing_BC(u_init, α)

        u_res_reflecting = CN_reflecting_BC(u_init, α)

        u_an_absorbing = absorbing_bounds.(L/2, x', t)

        u_an_reflecting= reflecting_bounds.(L/2, x', t)
        
        difference_abs = zeros(length(t)-1)
        difference_ref = zeros(length(t)-1)
        
        print(size(u_res_absorbing[end]))
        print(size(u_an_absorbing[end,:]))
        for i in length(t)-1
            difference_abs[i] = maximum(abs.(u_res_absorbing[i] .- u_an_absorbing[i,:]))
            difference_ref[i] = maximum(abs.(u_res_reflecting[i] .- u_an_reflecting[i,:]))
        end
        return maximum(difference_abs), maximum(difference_ref)
    end

end

# Diffusion.main()

println(Diffusion.test_timestep())