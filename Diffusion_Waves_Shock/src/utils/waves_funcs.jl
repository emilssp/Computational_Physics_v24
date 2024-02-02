function solve_u(u, timesteps, β)

    u[2:end-1, 2:end-1, 2] = u[2:end-1, 2:end-1, 1] + β * (u[2:end-1, 3:end, 1] - 2*u[2:end-1, 2:end-1, 1] + u[2:end-1, 1:end-2, 1] +
                                u[3:end, 2:end-1, 1] - 2*u[2:end-1, 2:end-1, 1] + u[1:end-2, 2:end-1, 1])

    for n in 2:timesteps-1
            u[2:end-1, 2:end-1, n+1] = 2*u[2:end-1, 2:end-1, n] - u[2:end-1, 2:end-1, n-1] + β * (u[2:end-1, 3:end, n] - 2*u[2:end-1, 2:end-1, n] + u[2:end-1, 1:end-2, n] +
                                                                        u[3:end, 2:end-1, n] - 2*u[2:end-1, 2:end-1, n] + u[1:end-2, 2:end-1, n])        
    end
    return u
end



function analytical_step(x,y,t)
    return cos(sqrt(5) * π * t) * sin(π * x) * sin(2 *π * y)
end



function analytical_solution(x, y, t, spacesteps, timesteps)
    u_an = zeros(spacesteps,spacesteps,timesteps)
    for n in 1:timesteps
        u_an[:,:,n] = analytical_step.(x',y,time[n])
    end  
    return u_an
end



function gaussian(x,y; x0 = 0.00, y0 = 0.00, σ = 1.00)
    return exp(-((x - x0)^2 +(y - y0)^ 2) / σ)
    
end