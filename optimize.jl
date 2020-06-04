#=
        project2.jl -- This is where the magic happens!

    All of your code must either live in this file, or be `include`d here.
=#

using LinearAlgebra
using Statistics
using Plots
#using Printf

global constraints_ineq
global constraints_eq
global constraints_eq_Ax_b
global nmax
global test_num = 1

visu_graph = true
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


function bfgs(prob, f0, g0, c0, f, x0; feas=nothing, ϵ=1e-4, α_max=2.5, max_bt_iters=10, method="bfgs", sherman=false)
	x = x0
	m = length(x0)
	Id = Matrix(I, m, m)
	Hinv = H = Id

	x_history = []
	push!(x_history, x)

	f_x = f(x)
	g_x = grad(f0, g0, c0, f, x, f_x)
	if g_x == Nothing
		return x, x_history
	end

	Ax_b = constraints_eq_Ax_b(x)
	if length(Ax_b) > 0
		A, b = Ax_b
		ATinv = A'
	end

	iter = 1

	while count(f0, g0, c0) < nmax - 20 && norm(g_x) > ϵ

		t1 = time_ns()

		if feas != nothing && feas(x)
			break # just used for the feasibility search phase
		end

		occursin("path1d", prob) && println("cost=$(f0(x)) counts=$(count(f0, g0, c0))")

		# 1) Compute search direction
		d = -Hinv * g_x
		if length(Ax_b) > 0
			d -= ATinv*(A*x-b) # cf Stephen Boyd book p.532
		end
		d = d./norm(d)

		#println("dsearch=$d g_x=$g_x")

		# 2) Compute step size
		bttime = @elapsed α = backtrack_line_search(f0, g0, c0, f, x, d, α_max, max_bt_iters, g_x)
		if α == Nothing
			break
		end
		#println("step size α=$α")

		xp = x + α .* d

		f_xp = f(xp)
		gradtime = @elapsed g_xp = grad(f0, g0, c0, f, xp, f_xp)
		if g_xp == Nothing
			break
		end

		δ = xp - x
		y = g_xp - g_x

		# 3) Update Hinv
		#den = y'⋅δ + 1e-5 # to avoid divide by 0 ...
		if method == "bfgs"
			den = max(1e-3, y'⋅δ)
			if sherman # Inverse via Sherman-Morisson formula
				Hinv = (Id - δ*y'/den) * Hinv * (Id - y*δ'/den) + δ*δ'/den 
			else
				den2 = max(1e-3, δ'*H*δ)
				H = H + y*y'/den - H*δ*δ'*H'/den2
				#H = H + y*y'/den - H*δ*δ'*H'/(δ'*H*δ + 1e-6) # Works fine as well

				if length(Ax_b) > 0
					# Cf Book "Convex Optimization" Ch.10
					#println("Ax-b = ", norm(A*x-b))
					rowsH, colsH = size(H)
					rowsA, colsA = size(A)
					n = rowsH + rowsA
					Heq = zeros(n, n)
					Heq[1:rowsH, 1:colsH] = H
					Heq[rowsH+1:end, 1:colsA] = A
					Heq[1:colsA, colsH+1:end] = A'
					invtime = @elapsed Heq_inv = inv(Heq)
					Hinv = Heq_inv[1:rowsH, 1:colsH]
					ATinv = Heq_inv[1:colsA, colsH+1:end]
				else
					Hinv = inv(H)
				end
			end
		elseif method == "broyden"
			den = max(1e-4, δ⋅δ')
			# More intuitive/obvious than bfgs
			# Just check that H*δ=y .... but not as good
			H = H + (y - H*δ)*δ'/den
			#H = H + (y - H*δ)*δ'/(δ'⋅δ+1e-7)
			Hinv = inv(H')
		end

		x, g_x = xp, g_xp
		push!(x_history, x)

		t2 = time_ns()

		println("Iter $(iter) BFGS=$((t2-t1)/1.0e6) ms: INV=$(invtime*1000) ms, BT=$(bttime*1000) ms, GRAD=$(gradtime*1000) ms")
		iter += 1
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

function p_quadratic(x)
	penalty = 0
	c_x = constraints_ineq(x)
	for i in 1:length(c_x)
		if c_x[i] > 0
			penalty += c_x[i]^2
		end
	end

	Ax_b = constraints_eq_Ax_b(x)
	if length(Ax_b) > 0
		A, b = Ax_b
		penalty += sum((A*x-b).^2)
	end

	#h_x = constraints_eq(x)
	#for i in 1:length(h_x)
	#	penalty += h_x[i]^2
	#end
	return penalty
end

function feasible(x)
	c_x = constraints_ineq(x)
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
	c_x = constraints_ineq(x)
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
	c_x = constraints_ineq(x)
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


function penalty_method(prob, f0, g0, c0, f, p, x, k_max; ρ=1, γ=2)
	x_history = []
	for k in 1:k_max
		println("ρ=$ρ")
		x, x_hist = bfgs(prob, f0, g0, c0, x -> f(x) + ρ*p(x), x; ϵ=1e-2, α_max=1.4, max_bt_iters=60)
		append!(x_history, x_hist)
		ρ *= γ
		if p(x) == 0
			return x
		end
	end
	return x, x_history
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

function optimize(f, g, c, h, h_Ax_b, x0, n, prob)
	x = x0
	global constraints_ineq = c
	global constraints_eq = h
	global constraints_eq_Ax_b = h_Ax_b
	global nmax = n

	x_history = []

	if occursin("penalty", prob)
		# ----------------
		# Penalty Method
		# ----------------
		x, x_history = penalty_method(prob, f, g, c, f, p_quadratic, x, 20; ρ=16, γ=2)
		x_feas = x
	else
		# ----------------------
		# Interior Point Method
		# ----------------------

		# Use BFGS during feasibility search (it is much better than CG: we have a quadratic penalty)
		x_feas, x_hist = bfgs(prob, f, g, c, p_quadratic, x; feas=feasible, α_max=1.5, max_bt_iters=100)
		append!(x_history, x_hist)

		push!(feas_scores, f(x_feas))
		push!(feas_counts, count(f,g,c))
		println("mean: feasibility score = $(mean(feas_scores)), count=$(mean(feas_counts))")

		x = x_feas

		# Interior Point Method
		prob=="simple2" ? ρ=10 : ρ=10#0
		ρ_max = 1e9
		delta = Inf
		while count(f, g, c) < nmax - 10 && delta > 1e-5 #&& ρ < ρ_max
			if prob=="simple2"
				x_int, x_hist = bfgs(prob, f, g, c, x -> f(x) + f_barrier(x)/ρ, x; ϵ=1e-2, α_max=1.4, max_bt_iters=30)
				ρ *= 5
			else
				#x_int, x_hist = bfgs(prob, f, g, c, x -> f(x) + f_barrier(x)/ρ, x; ϵ=1e-2, α_max=1.4, max_bt_iters=60)
				x_int, x_hist = bfgs(prob, f, g, c, x -> f(x) + f_barrier(x)/ρ, x; ϵ=1e-3, α_max=10.0, max_bt_iters=60)
				ρ *= 10
			end
			#println("fint($x_int)=$(f(x_int))")
			delta = norm(x_int - x)
			println("delta=$delta")
			x = x_int
			append!(x_history, x_hist)
		end
		println("max rho = $ρ")
	end

	push!(scores, f(x))
	push!(counts, count(f,g,c))
	println("$prob max_count=$(maximum(counts))")
	println("mean: final score = $(mean(scores)), count=$(mean(counts))")

	# ----------------
	# Visu Tool
	# ----------------

	if visu_graph && occursin("path1d", prob)
		global test_num

		function circleShape(t, s, Δt, Δs)
			θ = LinRange(0, 2*π, 500)
			t .+ Δt*sin.(θ), s .+ Δs*cos.(θ)
		end

		mpc = path1d_mpc()
		println("x_history[1]  : ", x_history[1])
		println("x_history[end]: ", round.(x_history[end]; digits=4))

		s = [x_history[end][i] for i in (1:3:length(x_history[end]))]
		t = collect((1:length(s))) .- 1
		t = t .* mpc.dt

		#println("s=", round.(s; digits=4))
		#println("t=$t")
		#println("obstacles=$(mpc.obstacles)")

		plot(t, s, title="s-t graph", label="final path", xlabel=("t (sec)"), ylabel=("s (m)"), marker=2, legend=:bottomright)
		for (i, obs) in enumerate(mpc.obstacles)
			tt, s = obs
			plot!(circleShape(tt, s, mpc.dt, mpc.dsaf), st=[:shape,], lw=0.5, linecolor=:black, fillalpha=0.2, label="obstacle $i")
		end

		if occursin("interior", prob)
			s = [x_feas[i] for i in (1:3:length(x_feas))]
			plot!(t, s, label="feasibility path", marker=2)
		end

		s = [x0[i] for i in (1:3:length(x0))]
		plot!(t, s, label="initial path", marker=2)

		savefig("plots/$(prob)_st_test$(test_num).png")
		test_num += 1
		if length(constraints_eq_Ax_b(x)) > 0
			A, b = constraints_eq_Ax_b(x)
			#println("equality constraints: $(constraints_eq(x))")
			println("equality constraints: Ax-b=$(A*x-b)")
			println("max   equality constraint violation:$(findmax(abs.(A*x-b))) out of $(length(b))")
		end
		println("max inequality constraint violation:$(findmax(constraints_ineq(x))) out of $(length(constraints_ineq(x)))")
		if mpc.slack_col
			nobs = length(mpc.obstacles)
			println("Collision avoidance constraint: slack_col=", round.(x[end-nobs+1:end], digits=2))
		end
	end

    return x
end
