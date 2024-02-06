
using Plots, LinearAlgebra, SpecialFunctions

function Crank_Nicolson_step(A::Tridiagonal{Float64, Vector{Float64}}, B::Tridiagonal{Float64, Vector{Float64}}, u_inint::Vector{Float64})
    return A\(B*u_inint)
end

function CN_absorbing_BC(u_init::Vector{Float64}, alpha::Float64 )
    #A * U_n+1 = B * U_n
    n_steps = size(u_init,1)-2

    A = Tridiagonal(fill(-alpha/2, n_steps-1),
                    fill(1+alpha, n_steps),
                    fill(-alpha/2, n_steps-1))

    B = Tridiagonal(fill(alpha/2, n_steps-1),
                    fill(1-alpha, n_steps),
                    fill(alpha/2, n_steps-1))


    res = Array{Array{Float64, 1}, 1}(undef, timesteps)
    res[1] = u_init
    for n in 1:timesteps-1
        res[n+1] = reduce(vcat, [0, Crank_Nicolson_step(A,B,res[n][2:end-1]),0])
    end
    return res
end 

function CN_reflecting_BC(u_init::Vector{Float64}, alpha::Float64)
    #A * U_n+1 = B * U_n
    n_steps = size(u_init,1)

    A = Tridiagonal(reduce(vcat, [fill(-alpha/2, n_steps-2), -alpha]),
                        fill(1+alpha, n_steps),
                        reduce(vcat, [-alpha, fill(-alpha/2, n_steps-2)]))
    
    B = Tridiagonal(reduce(vcat, [fill(alpha/2, n_steps-2), alpha]),
                        fill(1-alpha, n_steps),
                        reduce(vcat, [alpha, fill(alpha/2, n_steps-2)]))

    res = Array{Array{Float64, 1}, 1}(undef, timesteps)
    res[1] = u_init
    for n in 1:timesteps-1
        res[n+1] = Crank_Nicolson_step(A,B,res[n])
    end
    return res
end

function absorbing_bounds(x0::Float64, x::Float64, t::Float64; max_iter = 1000, tol = 1e-8)
    res = 0
    term = 1 
    for n in 1:max_iter
        vn_x0 = sqrt(2/L) * sin(n * pi * x0/L)
        vn_x = sqrt(2/L) * sin(n * pi * x/L)
        term = exp(-((n * pi/L)^2) * diffusivity * t) * vn_x0 * vn_x
        res += term
        if tol > abs(term) && abs(term) > 1e-15
            break
        end
    end
    return u0_tilda * res
end

function reflecting_bounds(x0::Float64, x::Float64, t::Float64; max_iter = 1000, tol = 1e-8)
    res = 1/L
    term = 1 
    for n in 1:max_iter
        vn_x0 = sqrt(2/L) * cos(n * pi * x0/L)
        vn_x = sqrt(2/L) * cos(n * pi * x/L)
        term = exp(-((n * pi/L)^2) * diffusivity * t) * vn_x0 * vn_x
        res += term
        if tol > abs(term) && abs(term) > 1e-15
            break
        end
    end
    return u0_tilda * res
end

function unbounded(x0::Float64, x::Float64, t::Float64)
    return u0_tilda/sqrt(4 * pi * diffusivity * t) * exp(-(x-x0)^2/(4* diffusivity * t))
end

function left_or_right(x,x_boundary)
    if x < x_boundary
        return true
    else x >= x_boundary
        return false
    end
end

function step_diffusivity(x::Float64,x_boundary::Float64)
    
    if left_or_right(x,x_boundary)
        return D_left*dt*spacesteps^2/L^2
    else left_or_right(x,x_boundary)
        return D_right*dt*spacesteps^2/L^2
    end
end

function make_diagonals(alpha::Vector{Float64}; method ='a')
    
    if method == 'a'
        n_steps = spacesteps-2
        
        LDA = zeros(n_steps-1) - (alpha[2:end-2] + alpha[3:end-1])/4
        MDA = ones(n_steps) + (alpha[1:end-2] + 2*alpha[2:end-1] + alpha[3:end])/4
        UDA = zeros(n_steps-1) - (alpha[2:end-2] + alpha[3:end-1])/4

        LDB = zeros(n_steps-1) + (alpha[2:end-2] + alpha[3:end-1])/4
        MDB = ones(n_steps) - (alpha[1:end-2] + 2*alpha[2:end-1] + alpha[3:end])/4
        UDB = zeros(n_steps-1) + (alpha[2:end-2] + alpha[3:end-1])/4

    end

    if method == 'r'
        n_steps = spacesteps
        left_imaginary = reduce(vcat, [alpha[1], alpha[1:end-1]])
        right_imaginary =reduce(vcat, [alpha[2:end], alpha[end]])
        reduce(vcat, [alpha[2:end]/2, alpha[end-1]/2])
        LDA = zeros(n_steps-1) - (alpha[1:end-1] + alpha[2:end])/4
        MDA = ones(n_steps) + (left_imaginary + 2*alpha + right_imaginary)/4
        UDA = zeros(n_steps-1) - (alpha[1:end-1] + alpha[2:end])/4

        LDB = zeros(n_steps-1) + (alpha[1:end-1] + alpha[2:end])/4
        MDB = ones(n_steps) - (left_imaginary + 2*alpha + right_imaginary)/4
        UDB = zeros(n_steps-1) + (alpha[1:end-1] + alpha[2:end])/4

        LDA[end]+= -(alpha[end] + alpha[end-1])/4
        UDA[1]  += -(alpha[1] + alpha[2])/4

        LDB[end]+= (alpha[end] + alpha[end-1])/4
        UDB[1]  += (alpha[1] + alpha[2])/4

    end

    return LDA, MDA, UDA, LDB, MDB, UDB
end

function absorbing_step(u_init::Vector{Float64}, alphas::Vector{Float64})

    n_steps = spacesteps - 2

    LDA, MDA, UDA, LDB, MDB, UDB = make_diagonals(alphas)
    
    A = Tridiagonal(LDA, MDA, UDA)

    B = Tridiagonal(LDB, MDB, UDB)


    res = Array{Array{Float64, 1}, 1}(undef, timesteps)
    res[1] = u_init
    for n in 1:timesteps-1
        res[n+1] = reduce(vcat, [0, Crank_Nicolson_step(A,B,res[n][2:end-1]),0])
    end
    return res
end 

function reflecting_step(u_init::Vector{Float64}, alphas::Vector{Float64})
    #A * U_n+1 = B * U_n
    n_steps = spacesteps

    LDA, MDA, UDA, LDB, MDB, UDB = make_diagonals(alphas; method = 'r')
    
    A = Tridiagonal(LDA, MDA, UDA)

    B = Tridiagonal(LDB, MDB, UDB)


    res = Array{Array{Float64, 1}, 1}(undef, timesteps)
    res[1] = u_init
    for i in 2:timesteps
        res[i] = Crank_Nicolson_step(A,B,res[i-1])
    end
    return res
end


function A_prefix_r(t, x0, x_boundary)

    erf_right = erf(x0 / sqrt(4*D_right*t))-1
    erf_left = erf(x0 / sqrt(4*D_left*t))-1

    exponent = (D_right - D_left)*(x0)^2/(4*D_right*D_left*t)
    
    denominator = 1 + erf_right + sqrt(D_left/D_right) * exp(exponent) * (1-erf_left)
    return 1/denominator
end

function A_prefix_l(t, x0, x_boundary)


    exponent = (D_right - D_left)*(x0)^2/(4*D_right*D_left*t)

    return A_prefix_r(t, x0, x_boundary)*sqrt(D_left/D_right)*exp(exponent)
    
end

function u_an_step_D(x, t, x0, x_boundary)
    # if t == 0 return 0 end
    if left_or_right(x,x_boundary)
        return u0_tilda*A_prefix_l(t, x0, x_boundary)*exp(-(x-x0)^2/(4*D_left*t))/sqrt(4*pi*D_left*t)
    else
        return u0_tilda*A_prefix_r(t, x0, x_boundary)*exp(-(x-x0)^2/(4*D_right*t))/sqrt(4*pi*D_right*t)
    end
    
end

function quadratic_func(x::Float64; C = 1.00)
    return (C*(x-1)^2 + (x-1) + 2)
end

function cubic_func(x::Float64; C = 1.00)
    return (C*(x-1)^3 + (x-1)^2 + (x-1) + 2)
end

function linear_func(x::Float64; C = 1.00)
    return C*x + 1
end

function variable_diffusivity(f, x; C = 1.00)
    return f.(x;C)
end

