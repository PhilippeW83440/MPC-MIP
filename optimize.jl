#=
        project2.jl -- This is where the magic happens!

    All of your code must either live in this file, or be `include`d here.
=#

using LinearAlgebra
using Statistics
#using Plots

global constraints
global nmax

plot_simple = false
npath = 1
barrier = "inv"
barrier = "log"


function backtrack_line_search(f0, g0, c0, f, x, d, α, max_bt_iters, g_x; p=0.5, β=1e-4)
	y = f(x)
	i = 0
	while i < max_bt_iters
		if count(f0, g0, c0) > nmax-5
			return Nothing
		end
		if f(x + α*d) > y + β*α*(g_x⋅d)
			α *= p
			i += 1
		else
			break
		end
	end
	#println(α)
	return α
end


function bfgs(prob, f0, g0, c0, f, x0, feas; ϵ=1e-3, α_max=2.5, max_bt_iters=10)
	x = x0
	m = length(x0)
	Id = Matrix(I, m, m)
	Hinv = Id

	x_history = []
	push!(x_history, x)

	f_x = f(x)
	g_x = grad(f0, g0, c0, f, x, f_x)
	if g_x == Nothing
		return x, x_history
	end

	while count(f0, g0, c0) < nmax - 20 && norm(g_x) > ϵ

		prob == "path1d" && println("cost=$(f0(x))")

		# 1) Compute search direction
		d = -Hinv * g_x
		d = d./norm(d)

		#println("dsearch=$d g_x=$g_x")

		# 2) Compute step size
		α = backtrack_line_search(f0, g0, c0, f, x, d, α_max, max_bt_iters, g_x)
		if α == Nothing
			break
		end
		#println("step size α=$α")

		xp = x + α .* d

		f_xp = f(xp)
		g_xp = grad(f0, g0, c0, f, xp, f_xp)
		if g_xp == Nothing
			break
		end

		δ = xp - x
		y = g_xp - g_x

		# 3) Update Hinv
		#den = y'⋅δ + 1e-5 # to avoid divide by 0 ...
		den = max(1e-3, y'⋅δ)
		Hinv = (Id - δ*y'/den) * Hinv * (Id - y*δ'/den) + δ*δ'/den 

		x, g_x = xp, g_xp
		push!(x_history, x)
	end

	return x, x_history
end

# was used for secret1 in place of bfgs
function conj_grad(prob, f0, g0, c0, f, x0, feas; ϵ=1e-3, α_max=2.5, max_bt_iters=10)
    x = x0
	iters = 0
	iters_feas = 0

	x_history = []
	push!(x_history, x)

	g_x = Inf

	while count(f0, g0, c0) < nmax-10 && norm(g_x) > ϵ
		f_x = f(x)
		println("iters=$iters: f($x) = $f_x")
		g_x = grad(f0, g0, c0, f, x, f_x)
		if g_x == Nothing
			break
		end

		# 1) Compute search direction
		if iters > 0
			#β = max(0, dot(g_x, g_x) / (g_x ⋅ gprev_x)) # Fletcher-Reeves
			β = max(0, dot(g_x, g_x - gprev_x) / (g_x ⋅ gprev_x)) # Polak-Ribiere
			d = -g_x + β*dprev
		else
			d = -g_x
		end
		d = d./norm(d)

		# 2) Compute step size
		α = backtrack_line_search(f0, g0, c0, f, x, d, α_max, max_bt_iters, g_x)
		if α == Nothing
			break
		end
		xp = x + α .* d

		if norm(xp - x) < 1e-5
			break
		end
		x = xp
		push!(x_history, x)

		global dprev, gprev_x = d, g_x
	end

	return x, x_history
end


# Used during feasibility search
function cg(f0, g0, c0, f, x0, n, feas; α_max=2.5, max_bt_iters=10)
    x = x0
	iters = 0
	iters_feas = 0

	x_history = []
	push!(x_history, x)

	while iters < n && !feas(x)
		f_x = f(x)
		println("iters=$iters: f($x) = $f_x")
		g_x = grad(f0, g0, c0, f, x, f_x)
		if g_x == Nothing
			break
		end

		# 1) Compute search direction
		if iters > 0
			#β = max(0, dot(g_x, g_x) / (g_x ⋅ gprev_x)) # Fletcher-Reeves
			β = max(0, dot(g_x, g_x - gprev_x) / (g_x ⋅ gprev_x)) # Polak-Ribiere
			d = -g_x + β*dprev
		else
			d = -g_x
		end
		d = d./norm(d)

		# 2) Compute step size
		α = backtrack_line_search(f0, g0, c0, f, x, d, α_max, max_bt_iters, g_x)
		if α == Nothing
			break
		end
		x = x + α .* d
		push!(x_history, x)

		global dprev, gprev_x = d, g_x
		iters += 1
	end

	return x, x_history
end


function grad(f0, g0, c0, f, x, f_x; h=sqrt(eps(Float64)))
#function grad(f0, g0, c0, f, x, f_x; h=cbrt(eps(Float64)))
	n = length(x)
	gradient = zeros(n)

	for i in 1:n
		u = zeros(n)
		u[i] = 1
		if count(f0, g0, c0) > nmax-5
			return Nothing
		end
		gradient[i] = (f(x+h*u)-f_x) / h
		#gradient[i] = (f(x+0.5*h*u)-f(x-0.5*h*u)) / h
	end
	#println("gradient=$gradient")
	return gradient
end

function f_penalty(x)
	penalty = 0
	c_x = constraints(x)
	for i in 1:length(c_x)
		if c_x[i] > 0
			penalty += c_x[i]^2
		end
	end
	return penalty
end

function feasible(x)
	c_x = constraints(x)
	for i in 1:length(c_x)
		if c_x[i] > 0
			return false
		end
	end
	return true
end

# --- Inverse Barrier for Interior Point method
function f_invb(x)
	res = 0
	c_x = constraints(x)
	for i in 1:length(c_x)
		if c_x[i] > 0
			return Inf
		else
			res -= 1/(c_x[i] + 1e-10)
		end
	end
	return res
end

# --- Log Barrier for Interior Point method
function f_logb(x)
	res = 0
	c_x = constraints(x)
	for i in 1:length(c_x)
		if c_x[i] > 0
			return Inf
		elseif c_x[i] >= -1.0
			res -= log(-c_x[i] + 1e-10)
		end
	end
	return res
end

if barrier == "log"
	f_barrier = f_logb
else
	f_barrier = f_invb
end


"""
    optimize(f, g, c, x0, n, prob)

Arguments:
    - `f`: Function to be optimized
    - `g`: Gradient function for `f`
    - `c`: Constraint function for 'f'
    - `x0`: (Vector) Initial position to start from
    - `n`: (Int) Number of evaluations allowed. Remember `g` costs twice of `f`
    - `prob`: (String) Name of the problem. So you can use a different strategy for each problem. E.g. "simple1", "secret2", etc.

Returns:
    - The location of the minimum
"""
feas_scores = []
feas_counts = []

scores = []
counts = []

function optimize(f, g, c, x0, n, prob)
    x_best = x0
	x = x0
	global constraints = c
	global nmax = n

	x_history = []

	x_feas, x_hist = cg(f, g, c, f_penalty, x, 200, feasible; α_max=1.5, max_bt_iters=100)
	append!(x_history, x_hist)

	push!(feas_scores, f(x_feas))
	push!(feas_counts, count(f,g,c))
	println("mean: feasibility score = $(mean(feas_scores)), count=$(mean(feas_counts))")

	x = x_feas

	# Interior Point Method
	prob=="simple2" ? ρ=10 : ρ=100
	ρ_max = 1e9
	delta = Inf
	while count(f, g, c) < nmax - 10 && delta > 1e-5 #&& ρ < ρ_max
		if prob=="simple2"
			x_int, x_hist = bfgs(prob, f, g, c, x -> f(x) + f_barrier(x)/ρ, x, feasible; ϵ=0.01, α_max=1.4, max_bt_iters=30)
			ρ *= 5
		else
			x_int, x_hist = bfgs(prob, f, g, c, x -> f(x) + f_barrier(x)/ρ, x, feasible; ϵ=1e-2, α_max=1.4, max_bt_iters=60)
			ρ *= 10
		end
		#println("fint($x_int)=$(f(x_int))")
		delta = norm(x_int - x)
		println("delta=$delta")
		x = x_int
		append!(x_history, x_hist)
	end

	push!(scores, f(x))
	push!(counts, count(f,g,c))
	println("$prob $(barrier)barrier max_count=$(maximum(counts))")
	println("mean: interior score = $(mean(scores)), count=$(mean(counts))")
	println("max rho = $ρ")

	prob=="path1d" && println("x_history: ", x_history[1])

	global npath
	global p1
	global p2
	global p3

	if plot_simple && (prob=="simple1" || prob=="simple2") && npath <= 3
		# Just for plotting
		function simple1()
			fobj(x) = -x[1]*x[2] +2/(3*sqrt(3))
			fmin(x1, x2) = -x1*x2 +2/(3*sqrt(3))
			c1(x1, x2) = x1 + x2^2 -1
			c2(x1, x2) = -x1 - x2
			#lev = [collect(0:1:10); collect(20:10:100); collect(120:20:300); collect(400:100:1000)]
			lev = 0:0.1:30
			return fobj, fmin, c1, c2, lev
		end
		
		# Just for plotting
		function simple2()
			fobj(x) = 100.0*(x[2]-x[1]^2)^2 + (1-x[1])^2
			fmin(x1, x2) = 100.0*(x2-x1^2)^2 + (1-x1)^2
			c1(x1, x2) = (x1-1)^3 -x2 +1
			c2(x1, x2) = x1 + x2 -2
			lev = [collect(0:2:10); collect(20:10:100); collect(120:20:300); collect(400:100:1000)]
			return fobj, fmin, c1, c2, lev
		end

		if prob=="simple1"
		 	fobj, fmin, c1, c2, lev = simple1()
		elseif prob=="simple2"
		 	fobj, fmin, c1, c2, lev = simple2()
		end

		f_history = [fobj(x) for x in x_history]
		c_history = [maximum(c(x)) for x in x_history]
		#println("$prob : ", x_history)
		println("$prob : ", c_history)
		println(length(x_history))

		if npath == 1
			x1 = -3:0.01:3
			x2 = -3:0.01:3

			p1 = plot(x1, x2, fmin, st=:contourf, c=:viridis, levels=lev, xaxis=("x1", (-3,3)), yaxis=("x2", (-3,3)), title="$prob $(barrier)barrier", contour_labels = true, legend=:bottomright)
			plot!(p1, x1, x2, c1, st=:contour, levels=0:.1:50, label="c1", c=:red, linealpha=0.5)
			plot!(p1, x1, x2, c2, st=:contour, levels=0:.1:50, label="c2", c=:red, linealpha=0.5)

			p2 = plot(f_history, yaxis=:log, title="Objective for $prob with $(barrier)barrier", label="Path $npath", xlabel="iterations", ylabel="objective", reuse=false)
			#p2 = plot(f_history, title="Objective for $prob with $(barrier)barrier", label="Path $npath", xlabel="iterations", ylabel="objective", reuse=false)
			p3 = plot(c_history, title="Max constraint for $prob with $(barrier)barrier", label="Path $npath", xlabel="iterations", ylabel="max constraint", reuse=false)
		else
			plot!(p2, f_history, label="Path $npath")
			plot!(p3, c_history, label="Path $npath")
		end
		#plot!(p1, Tuple.(x_history), label="path", c = :green1, marker=2)
		col = [:green1, :blue, :cyan]
		plot!(p1, Tuple.(x_history), label="Path $npath", c=col[npath], marker=2, linealpha=0.75)
		npath += 1

		if npath >= 4
			savefig(p1, "paths_$(prob)_$(barrier)barrier.png")
			savefig(p2, "obj_$(prob)_$(barrier)barrier.png")
			savefig(p3, "constr_$(prob)_$(barrier)barrier.png")
		end
	end

    return x
end
