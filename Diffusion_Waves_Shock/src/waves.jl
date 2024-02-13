module Waves
    
    using Plots
    include("./utils/waves_funcs.jl")

    const x0 = 0.5::Float64
    const y0 = 0.5::Float64
    const σ = 0.001::Float64

    const spacesteps = 100::Int64
    const L = 1.00::Float64
    const h = (L/spacesteps)::Float64

    const c = 1.00::Float64

    const endtime = 1.00::Float64
    const dt = (h/(2*c))::Float64#1e-4::Float64
    const timesteps = Int64(div(endtime,dt)+1)::Int64

    x = LinRange(0,L, spacesteps)
    y = LinRange(0,L, spacesteps)
    time = LinRange(0,endtime, timesteps)

    β = (c*dt/h)^2
    function print_cfl()
        println(c*dt/h)
    end
    ############################# Sinus initial values #################################
    function sinusoidal_IC()

        u = zeros(spacesteps, spacesteps, timesteps)

        u_init = sin.(pi*x') .* sin.(2*pi*y) 
        u_init[:,1] .= 0
        u_init[:,end] .= 0
        u_init[1,:] .= 0
        u_init[end,:] .= 0

        u[:,:,1] = u_init
        β = (c*dt/h)^2                
        u = solve_u(u,x,y, timesteps, β)

        u_an = analytical_solution(x,y,time,spacesteps, timesteps)
        anim = @animate for n in 1:10:timesteps
            p1 = surface(x, y, u_an[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nAnalytical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
            p2 = surface(x, y, u[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nNumerical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
            p = plot(p1,p2, layout= (1,2), size = (800,500), suptitle = "\nAnalitical and numerical solution at time t = $(round(time[n], digits = 4))")
        end

        difference =0
        temp = 0

        for i in length(x)
            for j in length(y)
                for t in length(time)
                    temp = abs.(u[i,j,t] .- u_an[i,j,t])
                    if temp>difference
                        difference = temp 
                    end
                end
            end
        end
        print("Max_err = $difference")

        gif(anim, "./figures/waves/wave_equation_solution_new.gif", fps=10)
        
        subplt1 = surface(x, y, u_an[:, :, round(Int, 0.4/dt)]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nAnalytical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        subplt2 = surface(x, y, u[:, :, round(Int, 0.4/dt)]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nNumerical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        plt = plot(subplt1, subplt2, layout= (1,2), size = (800,500), suptitle = "\nAnalitical and numerical solution at time t = $(round(time[round(Int, 0.4/dt)], digits = 4))")

        savefig(plt, "./figures/waves/sine_02.png")

    
    end

    ######################### Test stability #########################################################
    function test_stability()
        dt = 0.1
        endtime = 5.00
        timesteps = Int64(div(endtime,dt)+1)::Int64

        u = zeros(spacesteps, spacesteps, timesteps)

        u_init = sin.(pi*x') .* sin.(2*pi*y) 
        u_init[:,1] .= 0
        u_init[:,end] .= 0
        u_init[1,:] .= 0
        u_init[end,:] .= 0

        u[:,:,1] = u_init

        u = solve_u(u, timesteps, β)

        anim = @animate for n in 1:timesteps
            surface(x, y, u[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nNumerical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        end

        gif(anim, "./figures/waves/wave_equation_solution_test.gif", fps=10)
    end

    function membrane_wave()

        dt = 4*h/(2*c)
        print(dt)
        print(c)
        endtime = 3.00
        timesteps = Int64(div(endtime,dt)+1)::Int64

        time = LinRange(0,endtime, timesteps)

        u = zeros(spacesteps, spacesteps, timesteps)

        u_init = gaussian.(x',y; x0, y0, σ)
        u_init[:,1] .= 0
        u_init[:,end] .= 0
        u_init[1,:] .= 0
        u_init[end,:] .= 0

        u[:,:,1] = u_init

        u = solve_u(u, x, y, timesteps, β)
        u_an = analytical_solution(x,y,time, spacesteps,timesteps)

        anim3 = @animate for n in 1:timesteps
            surface(x, y, u[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "Membrane react to deformation t = $(round(time[n], digits = 2))", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        end
        gif(anim3, "./figures/waves/wave_equation_membrane_tst.gif", fps=20)
        
        
        subplt1 = surface(x, y, u[:, :, 1]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\n \$u(0,x,y)\$", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        subplt2 = surface(x, y, u[:, :, round(Int, 0.3/dt)]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\n \$u(0.3,x,y)\$", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        subplt3 = surface(x, y, u[:, :, round(Int, 1.2/dt)]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\n \$u(1.2,x,y)\$", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
        subplt4 = surface(x, y, u[:, :, round(Int, 2/dt)]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\n \$u(2.0,x,y)\$", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")


        plt = plot(subplt1, subplt2, layout= (1,2), size = (800,500), suptitle = "\nMembrane reacts to deformation")
        plt2 = plot(subplt3, subplt4, layout= (1,2), size = (800,500))
        savefig(plt, "./figures/waves/membrane_test.png")
        savefig(plt2, "./figures/waves/membrane2_test.png")

    end

end
Waves.sinusoidal_IC()
()