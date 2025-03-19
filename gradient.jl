using LinearAlgebra

function f(x :: Vector{Float64})
    x[1]^2 + 4*x[2]^2 + sin(6*x[1] + 7*x[2]) + 3*x[1] + 2*x[2]
end

function f2(x::Vector{Float64})
    (x[1]*x[2]-cos(x[1]))^2 + sinh(x[1]*x[2])^2;
end

function g(x :: Vector{Float64})
    res = zeros(2);
    res[1] = 2x[1] + 6*cos(6x[1] + 7x[2]) + 3;
    res[2] = 8x[2] + 7*cos(6x[1] + 7x[2]) + 2;
    return res;
end

function g2(xs::Vector{Float64})
    res = zeros(2);
    x = xs[1];
    y = xs[2];

    res[1] = 2*x*y^2 + 2*sin(x)*x*y - 2*cos(x)*y + sinh(2x*y)*y - sin(2x)
    res[2] = 2*x^2*y-2*cos(x)*x + sinh(2*x*y)*x
    return res;
end

function h(x :: Vector{Float64})
    res = zeros(2,2);
    res[1,1] = -36sin(6x[1] + 7x[2]) + 2;
    res[1,2] = res[2,1] = -42sin(6x[1] + 7x[2]);
    res[2,2] = -49sin(6x[1] + 7x[2]) + 8;
    return res;
end
function h2(xs :: Vector{Float64})
    res = zeros(2,2);
    x = xs[1];
    y = xs[2];

    res[1,1] = 2 * cos(x) * x*y + 2y^2 + 2 * cosh(2x*y) * y^2 + 4 * sin(x) * y − 2 * cos(2x)
    res[1,2] = res[2,1] = 6x*y + 4 * sinh(x*y)^2*x*y + 2 * sin(x) * x − 2 * cos(x) + sinh(2x*y)
	res[2,2] = 2x^2 + 2 * cosh(2x*y) * x^2
    return res;
end

sol = [-0.4050064442611152, -0.06937312042384142];
# sol =

function newton(f,g,h,x0 :: Vector{Float64},eps :: Float64)
    xk = x0;
    xp = x0;
    g_x = g(xk);
    while norm(g_x) > eps
        println("x: $xk\t∇ Norm: ", norm(g_x),"\tH(xk) eigenvalues: ",eigvals(h(xk)))
        p = h(xk)\ (-g_x)
        alpha = 1.0;
        # e = 0.5;
        # l = 0.67;
        # while (f(xk + p*alpha) - f(xk) > alpha * e * dot(g_x, p))
        #     alpha *= l
        # end
        xp = xk;
        xk = xk + p*alpha;
        g_x = g(xk);
        #println("Ratio: ",string(norm(xk - sol) / (norm(xp - sol)^2), "; ∇ Norm: ", norm(g_x)))
        #println("H(xk) eigenvalues: ",eigvals(h(xk)), "; R/m: ", 58.5/eigvals(h(xk))[1])
    end
    return xk;
end

print(newton(f2,g2,h2,[0.0,0.0], 1e-3))
