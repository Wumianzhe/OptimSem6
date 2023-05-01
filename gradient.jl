using LinearAlgebra

function f(x :: Vector{Float64})
    x[1]^2 + 4*x[2]^2 + sin(6*x[1] + 7*x[2]) + 3*x[1] + 2*x[2]
end

function g(x :: Vector{Float64})
    res = zeros(2);
    res[1] = 2x[1] + 6*cos(6x[1] + 7x[2]) + 3;
    res[2] = 8x[2] + 7*cos(6x[1] + 7x[2]) + 2;
    return res;
end

function h(x :: Vector{Float64})
    res = zeros(2,2);
    res[1,1] = -36sin(6x[1] + 7x[2]) + 2;
    res[1,2] = res[2,1] = -42sin(6x[1] + 7x[2]);
    res[2,2] = -49sin(6x[1] + 7x[2]) + 8;
    return res;
end

sol = [-1.9043886748, -0.3679467219];
# sol =

function newton(f,g,h,x0 :: Vector{Float64},eps :: Float64)
    xk = x0;
    xp = x0;
    g_x = g(xk);
    while norm(g_x) > eps
        p = h(xk)\ (-g_x)
        alpha = 1.0;
        e = 0.5;
        l = 0.667;
        while (f(xk + p*alpha) - f(xk) > alpha * e * dot(g_x, p))
            alpha *= l
        end
        xp = xk;
        xk = xk + p*alpha;
        g_x = g(xk);
        println(string(norm(xk - sol) / (norm(xp - sol)^2)))
    end
    return xk;
end

print(newton(f,g,h,[0.0,-0.0], 1e-9))
