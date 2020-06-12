# 2D Path Planning (before 1D, over S-T graph, Motion Planning)
# --- Symbolic computations to get analytical expressions ---

using SymPy

a0, a1, a2, a3, θ0 = @syms a0 a1 a2 a3 θ0 
s = symbols("s")

# f: is a symbolic expression, s: symbol, n: number
function integral_simpson(f, s, n)
    res =  f.subs(s, 0) + f
    for i in 1:2:n
        res += 4*f.subs(s, i*s/n)
	end
    for i in 2:2:n
        res += 2*f.subs(s, i*s/n)
	end
    res *= s/(3*n)
    return res
end

# Parametric curve used: polynomial spiral
K = a0 + a1*s + a2*s^2 + a3*s^3
θ = θ0 + integrate(K, s)
x = integral_simpson(cos(θ), s, 8)
y = integral_simpson(sin(θ), s, 8)

println("\n*** Path represented as a polynomial spiral ***\n")
println("curvature  : K=$K\n")
println("heading    : θ=$θ\n")
println("cartesian x: x=$x\n")
println("cartesian y: y=$y\n")

# Bending Energy distributes curvature evenly along spiral to promote comfort
fbe = integrate(K^2, s)
println("bending energy: fbe=$fbe\n")

# We switch from parameters a0, a1, a2, a3, sf -> p0, p1, p2, p3, p4
p0, p1, p2, p3, p4 = @syms p0 p1 p2 p3 p4 

# P0=K(0), p1=K(sf/3), p2=K(2*sf/3), p3=K(sf), p4=sf
eqns=[a0-p0, a3*(p4/3)^3+a2*(p4/3)^2+a1*(p4/3)+p0-p1, a3*(2*p4/3)^3+a2*(2*p4/3)^2+a1*(2*p4/3)+p0-p2, a3*(p4)^3+a2*(p4)^2+a1*(p4)+p0-p3]
convert = solve(eqns, [a0, a1, a2, a3])

println(convert)
println("\nWe auto-magically get the same result as in -Motion Planning for Autonomous Driving with a Conformal Spatio temporal Lattice- section 3 (ICRA 2011):\n")
for (key,val) in convert
	println("$key = $val")
end


# --- Optimization problem ---

# We want to find a parametric curve
# which minimizes fbe(a0, a1, a2, a3, sf)
# with constraints:
#   initial conditions (hard constraints): K0, θ0, x0, y0
#   final conditions   (soft constraints): Kf, θf, xf, yf
#   K(sf/3) <= Kmax
#   K(2*sf/3) <= Kmax

# We may add some more conditions ...

