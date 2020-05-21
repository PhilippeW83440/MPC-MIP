using LinearAlgebra
using Statistics

global constraints


function backtrack_line_search(f, x, d, α, max_bt_iters, g_x; p=0.5, β=1e-4)
	y = f(x)
	i = 0
	while i < max_bt_iters && f(x + α*d) > y + β*α*(g_x⋅d)
		α *= p
		i += 1
	end
	#println(α)
	return α
end

function bfgs(prob, f0, g0, c0, f, x0, feas; ϵ=1e-3, α_max=2.5, max_bt_iters=10)
	x = x0
	m = length(x0)
	Id = Matrix(I, m, m)
	Hinv = Id

	f_x = f(x)
	g_x = grad(f, x, f_x)

	if prob == "secret1"
		nmax = 1750 # 1750
	elseif prob =="secret2"
		nmax = 1900 # 1850 # 1500
	else
		nmax = 2000
	end

	while count(f0, g0, c0) < nmax - 100 && norm(g_x) > ϵ
		start_counts = count(f0, g0, c0)
		# 1) Compute search direction
		d = -Hinv * g_x
		d = d./norm(d)

		# 2) Compute step size
		α = backtrack_line_search(f, x, d, α_max, max_bt_iters, g_x)

		xp = x + α .* d

		f_xp = f(xp)
		g_xp = grad(f, xp, f_xp)

		δ = xp - x
		y = g_xp - g_x

		# 3) Update Hinv
		#den = y'⋅δ + 1e-5 # to avoid divide by 0 ...
		den = max(1e-3, y'⋅δ)
		#println("BFGS den=$den")
			# # SR-1 update
			# den = ((δ - Hinv * y)'⋅y)
			# if den > 1e-5
			# 	Hinv = Hinv + (δ - Hinv * y) * (δ - Hinv * y)' / den
			# end
		if prob != "secret1"
			Hinv = (Id - δ*y'/den) * Hinv * (Id - y*δ'/den) + δ*δ'/den 
		end

		x, g_x = xp, g_xp
		counts = count(f0, g0, c0)
		bfgs_counts = counts - start_counts
		if counts + 2 * bfgs_counts > nmax
			break
		end
	end

	return x
end

function conj_grad(prob, f0, g0, c0, f, x0, feas; ϵ=1e-3, α_max=2.5, max_bt_iters=10)
    x = x0
	iters = 0
	iters_feas = 0

	g_x = Inf

	if prob == "secret1"
		nmax = 1750 # 1750
	elseif prob =="secret2"
		nmax = 1850 # 1500
	else
		nmax = 2000
	end

	while count(f0, g0, c0) < nmax - 100 && norm(g_x) > ϵ
		start_counts = count(f0, g0, c0)
		f_x = f(x)
		println("iters=$iters: f($x) = $f_x")
		g_x = grad(f, x, f_x)

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
		α = backtrack_line_search(f, x, d, α_max, max_bt_iters, g_x)
		xp = x + α .* d

		if norm(xp - x) < 1e-5
			break
		end
		x = xp

		global dprev, gprev_x = d, g_x
		counts = count(f0, g0, c0)
		cg_counts = counts - start_counts
		if counts + 2 * cg_counts > nmax
			break
		end
	end

	return x
end


function cg(f, x0, n, feas; α_max=2.5, max_bt_iters=10)
    x = x0
	iters = 0
	iters_feas = 0

	while iters < n && !feas(x)
		f_x = f(x)
		println("iters=$iters: f($x) = $f_x")
		g_x = grad(f, x, f_x)

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
		α = backtrack_line_search(f, x, d, α_max, max_bt_iters, g_x)
		x = x + α .* d

		global dprev, gprev_x = d, g_x
		iters += 1
	end

	return x
end


# ---------------------------------------------------------------------

function grad(f, x, f_x; h=sqrt(eps(Float64)))
	n = length(x)
	gradient = zeros(n)

	for i in 1:n
		u = zeros(n)
		u[i] = 1
		gradient[i] = (f(x+h*u)-f_x) / h
		#gradient[i] = (f(x+0.5*h*u)-f(x-0.5*h*u)) / h
	end
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

	if prob == "secret1"
		x_best = x0 
		fx_best = Inf
		ndim = length(x0)
		xrand = x0
		ρ = 1e2
		while count(f, g, c) < 2000 - 300
			#x_int = bfgs(prob, f, g, c, x -> f(x) + f_invb(x)/ρ, x0, feasible; ϵ=0.01, α_max=7.4, max_bt_iters=50)
			x_int = conj_grad(prob, f, g, c, x -> f(x) + f_invb(x)/ρ, x0, feasible; ϵ=0.01, α_max=7.4, max_bt_iters=50)
			fx_int = f(x_int)
			if fx_int < fx_best
				fx_best = fx_int
				x_best = x_int
			end
			#xrand = 5*randn(ndim)
			#x0 = cg(f_penalty, xrand, 200, feasible; α_max=1.5, max_bt_iters=100)
			x0 = x_int
			ρ *= 10
		end
		x = x_best
	elseif prob == "secret2"
		x0 = cg(f_penalty, x0, 200, feasible; α_max=1.5, max_bt_iters=100)
		x_best = x0 
		fx_best = Inf
		ndim = length(x0)
		ρ = 1e2
		while count(f, g, c) < 2000 - 100 #- 300
			#x_int = bfgs(prob, f, g, c, x -> f(x) + f_logb(x)/ρ, x0, feasible; ϵ=0.01, α_max=1.4, max_bt_iters=70)
			x_int = bfgs(prob, f, g, c, x -> f(x) + f_logb(x)/ρ, x0, feasible; ϵ=1/ρ, α_max=1.5, max_bt_iters=70)
			fx_int = f(x_int)
			if fx_int < fx_best
				fx_best = fx_int
				x_best = x_int
			end
			x0 = x_int
			ρ *= 2
		end
		x = x_best
	else
		x_feas = cg(f_penalty, x, 200, feasible; α_max=1.5, max_bt_iters=100)

		push!(feas_scores, f(x_feas))
		push!(feas_counts, count(f,g,c))
		println("mean: feas score= $(mean(feas_scores)), count=$(mean(feas_counts))")

		x = x_feas

		# Interior Point Method
		prob=="simple2" ? ρ=10 : ρ=100
		ρ_max = 1e9
		delta = Inf
		while count(f, g, c) < 2000 && delta > 1e-3 && ρ < ρ_max
			#x_int = bfgs(prob, f, g, c, x -> f(x) + f_invb(x)/ρ, x, feasible; ϵ=1/ρ, α_max=1.4, max_bt_iters=60)
			if prob=="simple2"
				x_int = bfgs(prob, f, g, c, x -> f(x) + f_logb(x)/ρ, x, feasible; ϵ=0.01, α_max=1.4, max_bt_iters=30)
				#x_int = conj_grad(prob, f, g, c, x -> f(x) + f_logb(x)/ρ, x, feasible; ϵ=1e-2, α_max=1.4, max_bt_iters=30)
				ρ *= 5
			else
				x_int = bfgs(prob, f, g, c, x -> f(x) + f_logb(x)/ρ, x, feasible; ϵ=1e-2, α_max=1.4, max_bt_iters=60)
				#x_int = conj_grad(prob, f, g, c, x -> f(x) + f_logb(x)/ρ, x, feasible; ϵ=1e-2, α_max=1.4, max_bt_iters=60)
				ρ *= 10
			end
			#println("fint($x_int)=$(f(x_int))")
			delta = norm(x_int - x)
			x = x_int
		end
	end

	push!(scores, f(x))
	push!(counts, count(f,g,c))
	println("mean: int score= $(mean(scores)), count=$(mean(counts)) max_count=$(maximum(counts))")
	println("$ρ")

    return x
end
