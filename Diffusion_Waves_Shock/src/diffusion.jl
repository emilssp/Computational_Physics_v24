module Diffusion

include("utils/diffusion_funcs.jl")

    const diffusivity = 1.00::Float64
    
    const endtime = 1.00::Float64
    const L = 1.00::Float64

    const dt = 1e-4::Float64
    const timesteps = Int64(div(endtime,dt)+1)::Int64

    const spacesteps = 150::Int64

    const u0_tilda = 1::Int64
    const u0 = (u0_tilda*spacesteps/L)::Float64

    const cfl = (diffusivity*(spacesteps*spacesteps)*dt/L^2)::Float64

    const D_left = 0.50::Float64
    const D_right = 1.00::Float64
    
    const subdir = "figures/diffusion/"

    ##########################################################################

    function plot_array(vals_CN, vals_an, time_arr, x, t; figname = "const_D/new_fig.png")
        plt = plot()
        for time in time_arr
            scatter!(x,vals_CN[round(Int, time/dt),:], label = t[round(Int, time/dt)])
            plot!(x,vals_an[round(Int, time/dt),:], color = "black", label = "")#t[round(Int, time/dt)])
        end
        scatter!(x,vals_CN[end,:], label = t[end])
        plot!(x,vals_an[end,:], color = "black", label = "")#t[round(Int, time/dt)])
        savefig(plt, joinpath(subdir,figname))
        return plt 
    end

    function main()

        x = reduce(vcat, [0, cumsum(fill(L/(spacesteps-1), spacesteps-1))])
        t = LinRange(0,endtime, timesteps)

        u_init = zeros(spacesteps)
        u_init[div(spacesteps, 2)] = u0/2
        u_init[div(spacesteps, 2)+1] = u0/2


        u_res_absorbing = CN_absorbing_BC(u_init, cfl)

        u_res_reflecting = CN_reflecting_BC(u_init, cfl)

        u_an_unbounded = unbounded.(L/2, x', t[2:end])

        u_an_absorbing = absorbing_bounds.(L/2, x', t)

        u_an_reflecting= reflecting_bounds.(L/2, x', t)

        println("Done")

        plot_array(u_res_absorbing, u_an_unbounded, [0.01,0.04,0.07], x, t; figname = "const_D/u_res_absorbing")
        plot_array(u_res_reflecting, u_an_unbounded, [0.01,0.04,0.07], x, t; figname = "const_D/u_res_reflecting")
        plot_array(u_res_absorbing, u_an_absorbing, [0.01,0.04,0.07], x, t; figname = "const_D/u_an_absorbing")
        plot_array(u_res_reflecting, u_an_reflecting, [0.01,0.04,0.07], x, t; figname = "const_D/u_an_reflecting")

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

end