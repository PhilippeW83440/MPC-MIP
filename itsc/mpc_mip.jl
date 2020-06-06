using Convex
using Mosek, MosekTools # License required www.mosek.com
using LinearAlgebra

using GLPK
using ECOS

#using SCS
#using COSMO
#using CSDP
#using SDPA
#using OSQP
#
#using Gurobi # License required
#using CPLEX # License required

using Plots

testu = false
first_plot = true

solvers = Dict("Mosek" => ()->Mosek.Optimizer(),
               "GLPK" => ()->GLPK.Optimizer(method = GLPK.SIMPLEX), # does not support QP obj
               #"SCS" => ()->SCS.Optimizer(), 
               #"COSMO" => ()->COSMO.Optimizer(),
               #"CSDP" => ()->CSDP.Optimizer(),
               #"SDPA" => ()->SDPA.Optimizer(),
               #"CPLEX" => ()->CPLEX.Optimizer(),
               #"Gurobi" => ()->Gurobi.Optimizer(),
               #"OSQP" => ()->OSQP.Optimizer(),
               "ECOS" => ()->ECOS.Optimizer())

# Use Mosek or ECOS: they are very fast (<20 ms and < 10 ms) and good ...
# ECOS is free but does not support Mixed Integer Programming

pkg = "Mosek"
solver = solvers[pkg]

# For Minnie Simplex tests
testu, pkg = true, "GLPK"


# Hyperparameter(s) of the Collision Avoidance Model
MARGIN_SAF = 1.1

# -----------------------------------
# MPC with Mixed Integer Programming
# -----------------------------------

# We use Mosek to deal with Quadratic cost fct + MIP
# Open source GPLK does not deal with Quadratic cost fct
# Another commercial alternative could be: Gurobi

# We will use Boolean Variables to deal with disjunctive constraints (OR)
# Cf https://optimization.mccormick.northwestern.edu/index.php/Disjunctive_inequalities


mutable struct MpcPath1d
	# MPC planning horizon
	T::Int64 # T horizon: number of time steps
	dt::Float64 # Time Step length

	# Quadratic cost function
	Q::Array{Float64,2}
	R::Array{Float64,2}

	# Feasibility constraints
	smin::Float64 # min pos
	smax::Float64 # max pos
	vmin::Float64 # min speed
	vmax::Float64 # max speed
	umin::Float64 # min acceleration/command
	umax::Float64 # max acceleration/command

	# Safety constraints
	dsaf::Float64 # safety distance wrt other vehicles

	# Dynamic constraints (Constant Acceleration model)
	Ad::Array{Float64,2} # [1 dt; 0 1]
	Bd::Array{Float64,1} # [dt^2/2 dt]

	# # --- NB: note that so far, all constraints are LINEAR ---

	xref::Vector{Float64} # our target position and speed (ref)
	uref::Vector{Float64} # our target control (want to stabilize at u=0)
	xinit::Vector{Float64} # our initial position and speed

	# NB: note that without the obstacles constraints
	# we could deal with linear obj and linear constraints
	# => a simplex algo would be guaranteed to find global optimum (TODO ?)

	# Obstacles: list (1D array) of (t,s) tuples
	# Define when and where an object will cross our path
	# This is based on a prediction: 1 obj may lead to multiple (t,s) tuples
	# to account for uncertainty. NB: could be more efficient to increase dsaf
	obstacles::Array{Tuple{Float64, Float64}, 1}

	nvars_dt::Int64 # Number of variables PER Time Step (2 state vars + 1 ctlr var)

	solv

	MpcPath1d(T=20, dt=0.250,
			  Q=Diagonal([0.2, 10]), R=Diagonal([0.001]), # Normalize <=10 the coefficients
			  smin=0.0, smax=1000.0,
			  vmin=0.0, vmax=25.0, # Between 0 and 25 m.s-1
			  umin=-4.0, umax=2.0, # Between -4 and 2 m.s-2
			  dsaf=10.0,
			  Ad=[1.0 dt; 0 1],
			  Bd=[0.5*dt^2, dt],

			  # NB: xref, uref, xinit, obstacles are expected to be changed/updated everytime we plan (TODO)
			  xref=[100.0, 20.0], # Target pos=100 at v=20 m.s-1
			  uref=[0.0],
			  xinit=[0.0, 20.0], # Start at s=0 with speed v=20 m.s-1
			  obstacles=[(2.0, 40), (3.0, 80)], # In 2 sec a crossing vehicle at s=40 m
			  #obstacles=[(2.0, 35), (3.0, 65)], # In 2 sec a crossing vehicle at s=40 m

			  nvars_dt=3, # x=[s,sd] u=[sdd]
			  solv=Nothing
			  ) = new(T,dt,Q,R,smin,smax,vmin,vmax,umin,umax,dsaf,Ad,Bd,xref,uref,xinit,obstacles,nvars_dt,solv)
end


function path1d_cost(mpc::MpcPath1d, x::Variable, slack_col)
	cost = 0
	# stage cost
	for k in range(1, step=mpc.nvars_dt, length=mpc.T-1)
		xk, uk = x[k:k+1], x[k+2]
		cost += quadform(xk - mpc.xref, mpc.Q)
		cost += quadform(uk - mpc.uref, mpc.R)
	end

	xT, uT = x[end-2:end-1], [x[end]]
	cost += quadform(xT - mpc.xref, mpc.Q)
	cost /= (mpc.nvars_dt * mpc.T) # make the cost num_steps independant

	cost += 1e2*sum(slack_col) # Elastic Model for Collision Avoidance
	return cost
end

function path1d_constraints(mpc::MpcPath1d, x::Variable, p, slack_col, slack_bin)

	# ------------ Constraints that are fixed once for all ----------------

	dt = mpc.dt
	# Time Steps: 2 .. T Dynamics Cosntraints
	for k in range(1+mpc.nvars_dt, step=mpc.nvars_dt, length=mpc.T-1)
		# Dynamic Constraints: Constant Acceleration Model in between 2 Time Steps
		kp = k - mpc.nvars_dt # p for previous
		p.constraints += [x[k] == x[kp] + dt*x[kp+1] + 0.5*dt^2*x[kp+2]]
		p.constraints += [x[k+1] == x[kp+1] + dt*x[kp+2]]
	end

	# --- Inequality constraints ---
	p.constraints += [mpc.umin <= x[3], x[3] <= mpc.umax]

	# Time Steps: 2 .. T
	for k in range(1+mpc.nvars_dt, step=mpc.nvars_dt, length=mpc.T-1)
		p.constraints += [mpc.smin <= x[k], x[k] <= mpc.smax]
		p.constraints += [mpc.vmin <= x[k+1], x[k+1] <= mpc.vmax]
		p.constraints += [mpc.umin <= x[k+2], x[k+2] <= mpc.umax]
	end

	#println("BEFORE: ", length(p.constraints))

	# ------------ Constraints that can be changed on the fly ----------------

	# Order matters !!! comply to a LIFO
	# https://docs.mosek.com/9.2/capi/guidelines-optimizer.html?highlight=lifo

	# --- Equality constraints ---
	p.constraints += [x[1] == mpc.xinit[1]]
	p.constraints += [x[2] == mpc.xinit[2]]

	# Handle obstacle constraint: Elastic Model + Disjunctive Constraints
	M = 1e4
	dsaf = MARGIN_SAF * mpc.dsaf
	for (i, obstacle) in enumerate(mpc.obstacles)
		tcross, scross = obstacle
		tcrossd = floor(Int, tcross/mpc.dt)
		if (tcrossd >= 1) && (tcrossd <= mpc.T)
			# if within MPC horizon
			# BUG FIX: NOT -1 ... Index starts at 1 in Julia ... 
			k = tcrossd * mpc.nvars_dt + 1

			# Elastic Model for Collision Avoidance
			#p.constraints += [x[k] <= scross - dsaf + slack_col[i]]
			p.constraints += [x[k] <= scross - dsaf + slack_col[i] + M * slack_bin[i]]
			p.constraints += [scross + dsaf - slack_col[i] <= x[k] + M * (1 - slack_bin[i]) ]
			p.constraints += [0 <= slack_col[i], slack_col[i] <= dsaf]
		end
	end

end

function path1d_init(mpc::MpcPath1d)
	x0 = zeros(mpc.T * mpc.nvars_dt)
	v = mpc.xinit[2]

	# Time Step: 1
	x0[1:2], x0[3] = mpc.xinit, 0 # init [pos,speed], [accel]

	# Time Steps: 2 .. T
	for k in range(1+mpc.nvars_dt, step=mpc.nvars_dt, length=mpc.T-1)
		x0[k] = x0[k - mpc.nvars_dt] + v * mpc.dt
		x0[k+1] = v
		x0[k+2] = 0
	end
	return x0
end


mutable struct SolverData
	x::Variable
	# Slack Variables for Collision Avoidance
	slack_col::Vector{Variable} # Elastic Model
	slack_bin::Vector{Variable} # Disjunctive Constraints
	p::Problem{Float64}
	num_mpc_obstacles::Int64 # The obstacles visible within MPC time horizon

	function SolverData(mpc)
		s = new()

		s.x = Variable(60)
		s.slack_col, s.slack_bin = [], []
		#for i in 1:length(mpc.obstacles)
		for i in 1:10 # Account (pre allocate) for future obstacles
			push!(s.slack_col, Variable(1, Positive()))
			push!(s.slack_bin, Variable(1, :Bin))
		end

		# GLPK Simplex can be used to find a feasible solution in <= 4ms
		# Could be used then as a starting point for an interior point method
		pkg == "GLPK" ? cost = sum(s.slack_col) : cost = path1d_cost(mpc, s.x, s.slack_col)
		s.p = minimize(cost)
		path1d_constraints(mpc, s.x, s.p, s.slack_col, s.slack_bin)
		s.num_mpc_obstacles = length(mpc.obstacles)

		return s
	end
end

# Dump solution found (text)
function dump_sol(solv::SolverData)
	println("status=", solv.p.status)
	println("cost=", round(solv.p.optval, digits=2))
	println("x=", round.(solv.x.value, digits=2))
	#for (j, col) in enumerate(solv.slack_col)
	for j in 1:solv.num_mpc_obstacles
		println("Collision avoidance constraint $j: slack_col=", round(solv.slack_col[j].value, digits=2))
		println("Collision avoidance        bin $j: slack_bin=", solv.slack_bin[j].value)
	end
end

# Dump solution found (graph)
function dump_graph(x0::Vector{Float64}, solv::SolverData, mpc::MpcPath1d, runtime::Int64)
	function circleShape(t, s, Δt, Δs)
		θ = LinRange(0, 2*π, 500)
		t .+ Δt*sin.(θ), s .+ Δs*cos.(θ)
	end
	
	global first_plot
	x = solv.x.value
	
	s = [x[i] for i in (1:3:length(x))]
	t = collect((1:length(s))) .- 1
	t = t .* mpc.dt
	
	plot!(t, s, label="$(pkg) path ($runtime ms)", marker=2)
	
	if first_plot
		for (i, obs) in enumerate(mpc.obstacles)
			tt, s = obs
			plot!(circleShape(tt, s, mpc.dt, mpc.dsaf), st=[:shape,], lw=0.5, linecolor=:black, fillalpha=0.2, label="obstacle $i")
		end
	
		s = [x0[i] for i in (1:3:length(x0))]
		plot!(t, s, label="initial path", marker=2)
		first_plot = false
	end

end


# ------------- Main for unitary tests --------------------

if testu
	mpc = MpcPath1d()
	solv = SolverData(mpc)
	
	plot(title="s-t graph", xlabel=("t (sec)"), ylabel=("s (m)"), marker=2, legend=:topleft)
	
	#for (pkg, solver) in solvers
	for justone in 1:1
		x0 = path1d_init(mpc)
		runtime = 0
	
		for i in 1:2 # 2 times to get runtime for compiled version
			println("start solver...")
			solv.x.value = x0
			println("x0=$x0")
		
			runtime = @elapsed solve!(solv.p, solver, warmstart=true)
			runtime = convert(Int64, round(1000*runtime))
			println("runtime=$runtime ms")
	
			dump_sol(solv)
		end
		
		dump_graph(x0, solv, mpc, runtime)
	end
	
	#savefig("st_graph_solvers.png")
	savefig("st_graph_$(pkg).png")
end


# Hooks for scn.jl
# -------------------------------------------------------------------------
function mpc_mip()
	mpc = MpcPath1d()
	solv = SolverData(mpc)
	mpc.solv = solv

	# First run (fake one to enable warmstart then)
	x0 = path1d_init(mpc)
	solv.x.value = x0
	solve!(solv.p, solver, warmstart=true)
	dump_sol(solv)

	p = solv.p
	println("#constraints = ", length(p.constraints))

	# 2 initial conditions constraints to delete (new ones needed at each act() call)
	pop!(p.constraints) # p.constraints += [x[1] == mpc.xinit[1]]
	pop!(p.constraints) # p.constraints += [x[2] == mpc.xinit[2]]

	for i in 1:solv.num_mpc_obstacles
		# 4 constraints per obstacle to delete (new ones needed at each act() call)
		#p.constraints += [x[k] <= scross - dsaf + slack_col[i] + M * slack_bin[i]]
		#p.constraints += [scross + dsaf - slack_col[i] <= x[k] + M * (1 - slack_bin[i]) ]
		#p.constraints += [0 <= slack_col[i], slack_col[i] <= dsaf]
		for j in 1:4
			pop!(p.constraints)
		end
	end
	solv.num_mpc_obstacles = 0

	println("#constraints = ", length(p.constraints))
	return mpc
end

function act(mpc::MpcPath1d, state::Array{Float64}, obstacles::Array{Tuple{Float64, Float64}, 1})::Float64

	#println("state=$state")
	#println("obstacles=$obstacles")

	# retrieve solver setup
	solv = mpc.solv
	x, p = solv.x, solv.p
	slack_col, slack_bin = solv.slack_col, solv.slack_bin

	# Update current initial condition pos and speed
	for i in 1:2
		mpc.xinit[i] = state[i]
		p.constraints += [x[i] == mpc.xinit[i]]
	end
	mpc.obstacles = obstacles

	# Update obstacles constraints: Elastic Model + Disjunctive Constraints
	#mpc.obstacles = obstacles
	M = 1e4
	dsaf = MARGIN_SAF * mpc.dsaf
	solv.num_mpc_obstacles = 0
	for (i, obstacle) in enumerate(obstacles)
		tcross, scross = obstacle
		tcrossd = floor(Int, (tcross-state[3])/mpc.dt) # relative time
		if (tcrossd > 0) && (tcrossd < mpc.T)
			solv.num_mpc_obstacles += 1
			# if within MPC horizon
			# BUG FIX: NOT -1 ... Index starts at 1 in Julia ... 
			k = tcrossd * mpc.nvars_dt + 1
			#println("DBG k=$k")

			# Elastic Model for Collision Avoidance
			#p.constraints += [x[k] <= scross - dsaf + slack_col[i]]
			p.constraints += [x[k] <= scross - dsaf + slack_col[i] + M * slack_bin[i]]
			p.constraints += [scross + dsaf - slack_col[i] <= x[k] + M * (1 - slack_bin[i]) ]
			p.constraints += [0 <= slack_col[i], slack_col[i] <= dsaf]
		end
	end

	println("num_mpc_obstacles=$(solv.num_mpc_obstacles)")

	# Update x0
	x0 = path1d_init(mpc)
	x.value = x0

	println("Attempting solve ...")
	solve!(solv.p, solver, warmstart=true)
	println("Finished solve")
	dump_sol(solv)

	println("#constraints = ", length(p.constraints))

	# TODO delete init+obstacles constraints
	# 2 initial conditions constraints to delete (new ones needed at each act() call)
	pop!(p.constraints) # p.constraints += [x[1] == mpc.xinit[1]]
	pop!(p.constraints) # p.constraints += [x[2] == mpc.xinit[2]]

	for i in 1:solv.num_mpc_obstacles
		# 4 constraints per obstacle to delete (new ones needed at each act() call)
		#p.constraints += [x[k] <= scross - dsaf + slack_col[i] + M * slack_bin[i]]
		#p.constraints += [scross + dsaf - slack_col[i] <= x[k] + M * (1 - slack_bin[i]) ]
		#p.constraints += [0 <= slack_col[i], slack_col[i] <= dsaf]
		for j in 1:4
			pop!(p.constraints)
		end
	end
	solv.num_mpc_obstacles = 0

	println("#constraints = ", length(p.constraints))
	println("x = $(x.value)")
	println("---------------------------------> cmd = $(x.value[3])")

	return x.value[3]
end
# -------------------------------------------------------------------------



