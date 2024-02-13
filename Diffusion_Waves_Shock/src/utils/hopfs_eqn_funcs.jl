function gaussian(x::Float64; x0 = 0.4, σ = 0.03)
    return exp(-(x - x0)^2 / σ)
end

function shock(x)
    if x<1 return 1 - 0.5*x
    else return 0.5*x - 0.5
    end
end
function step_func(x; p = 1)
    if 0.3<x<0.5
        return p
    else return 0
    end
end

function sigmoid(x)
    return 1/(1+exp((-x+0.5)*16))
end

function Lax_Wandroff_trans(u_init, spacesteps, timesteps; γ = 1.00)
    u = zeros(spacesteps, timesteps)
    u[:,1] = u_init
    for n in 1:timesteps-1
        u[2:end-1, n+1] = @. u[2:end-1, n] - 0.5 * γ * (u[3:end, n] - u[1:end-2, n]) + 0.5 * γ^2 * (u[3:end, n] - 2 * u[2:end-1,n] + u[1:end-2, n]) 
        u[end, n] = u[1,n] = 0
    end
    return u
end

function trans_anal(f::Function, x::LinRange{Float64, Int64}, t::LinRange{Float64, Int64}; c = 1)
    arg = x .- c .* t'
    return f.(arg)

end

function u0(x)
    return max(1 - (4 * (x - 0.5))^2, 0)
end

function Lax_Wandroff_Hopf(u_init, spacesteps, timesteps; γ = 1.00)
    u = zeros(spacesteps,timesteps)
    u[:, 1] = u_init
    for n in 1:timesteps-1
        u[2:end-1, n+1] = @. u[2:end-1, n] - γ/4 * ((u[3:end, n]) ^ 2 - (u[1:end-2, n])^2) + (γ^2 / 8) *((u[3:end, n] + u[2:end-1, n]) * ((u[3:end, n])^2 - (u[2:end-1, n])^2) -(u[2:end-1, n] + u[1:end-2, n]) * ((u[2:end-1, n])^2 - (u[1:end-2, n])^2))
        u[1, n+1] = u[2, n+1]
        u[end, n+1] = u[end-1, n+1]                 
    end
    return u
end