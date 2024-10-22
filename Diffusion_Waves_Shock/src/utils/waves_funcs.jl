function solve_u(u::Array{Float64,3}, x, y, timesteps::Int64, β)
    center = [0.5,0.5]
    radius = 0.5
    mask = ((x' .- center[1]).^2 .+ (y .- center[2]).^2 .<= radius^2)
    u[2:end-1, 2:end-1, 2] = u[2:end-1, 2:end-1, 1] + β * (u[2:end-1, 3:end, 1] - 2*u[2:end-1, 2:end-1, 1] + u[2:end-1, 1:end-2, 1] +
                                u[3:end, 2:end-1, 1] - 2*u[2:end-1, 2:end-1, 1] + u[1:end-2, 2:end-1, 1])

    for n in 2:timesteps-1
        u[2:end-1, 2:end-1, n+1] = 2*u[2:end-1, 2:end-1, n] - u[2:end-1, 2:end-1, n-1] + β * (u[2:end-1, 3:end, n] - 2*u[2:end-1, 2:end-1, n] + u[2:end-1, 1:end-2, n] +
                                                                    u[3:end, 2:end-1, n] - 2*u[2:end-1, 2:end-1, n] + u[1:end-2, 2:end-1, n])  
        u[:,:,n+1] =u[:,:,n+1] .* mask      
    end
    return u
end



function analytical_step(x::Float64,y::Float64,t::Float64)
    return cos(sqrt(5) * π * t) * sin(π * x) * sin(2 *π * y)
end



function analytical_solution(x::LinRange{Float64, Int64}, y::LinRange{Float64, Int64}, t::LinRange{Float64, Int64}, spacesteps::Int64, timesteps::Int64)
    u_an = zeros(spacesteps,spacesteps,timesteps)
    for n in 1:timesteps
        u_an[:,:,n] = analytical_step.(x',y,t[n])
    end  
    return u_an
end



function gaussian(x::Float64, y::Float64; x0 = 0.00, y0 = 0.00, σ = 1.00)
    return exp(-((x - x0)^2 +(y - y0)^ 2) / σ)
end