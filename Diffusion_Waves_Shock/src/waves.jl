module Waves
    
using Plots
include("./utils/waves_funcs.jl")


const spacesteps = 100::Int64
const L = 1.00::Float64
const h = (L/spacesteps)::Float64

const c = 1.00::Float64

const endtime = 1.00::Float64
const dt = h/(2*c)#1e-4::Float64
const timesteps = Int64(div(endtime,dt)+1)::Int64


############################# Sinus initial values #################################

x = reduce(vcat, [0, cumsum(fill(L/(spacesteps-1), spacesteps-1))])
y = reduce(vcat, [0, cumsum(fill(L/(spacesteps-1), spacesteps-1))])

time = LinRange(0,endtime, timesteps)

β = (c*dt/h)^2

u = zeros(spacesteps, spacesteps, timesteps)

u_init = sin.(pi*x') .* sin.(2*pi*y) 
u_init[:,1] .= 0
u_init[:,end] .= 0
u_init[1,:] .= 0
u_init[end,:] .= 0

u[:,:,1] = u_init
                    
u = solve_u(u, timesteps, β)
u_an = analytical_solution(x,y,time,spacesteps, timesteps)
anim = @animate for n in 1:10:timesteps
    p1 = surface(x, y, u_an[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nAnalytical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
    p2 = surface(x, y, u[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nNumerical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
    p = plot(p1,p2, layout= (1,2), size = (800,500), suptitle = "\nAnalitical and numerical solution at time t = $(round(time[n], digits = 4))")
end

gif(anim, "figures/waves/wave_equation_solution.gif", fps=10)


######################### Test stability #########################################################

dt = 1.45h/(sqrt(2)*c)
endtime = 1.00
timesteps = Int64(div(endtime,dt)+1)::Int64

time = LinRange(0,endtime, timesteps)

u = zeros(spacesteps, spacesteps, timesteps)

u_init = sin.(pi*x') .* sin.(2*pi*y) 
u_init[:,1] .= 0
u_init[:,end] .= 0
u_init[1,:] .= 0
u_init[end,:] .= 0

u[:,:,1] = u_init

u = solve_u(u, timesteps, β)

anim = @animate for n in 1:timesteps
    p2 = surface(x, y, u[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "\n\nNumerical solution", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
end

gif(anim, "figures/waves/wave_equation_solution_test.gif", fps=10)

const x0 = 0.5::Float64
const y0 = 0.5::Float64
const σ = 0.001::Float64

dt = h/(2*c)
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

u = solve_u(u, timesteps, β)
u_an = analytical_solution(x,y,time, spacesteps,timesteps)

anim3 = @animate for n in 1:timesteps
    p3 = surface(x, y, u[:, :, n]', zlim=(-1, 1), c=:thermal, colorbar = false, title = "Membrane react to deformation t = $(round(time[n], digits = 2))", xlabel = "\$x\$ (a.u.)", ylabel = "\$y\$ (a.u.)", zlabel = "\$u(x,y,t)\$")
end
gif(anim3, "figures/waves/wave_equation_membrane.gif", fps=20)

end