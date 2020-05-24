########### "Path1d" Problem Definition ###########

using LinearAlgebra

# via Optimization over station

# ----------------------
# Problem Specification:
# ----------------------

# Variables:
# ----------

# We consider a Time Step of 250 ms
# We plan with an horizon over 5 seconds
# => T=20 time steps

# 20 Time Steps => k in [1, 20]
# state x_k=[s,sd], control u_k=[sdd] for k in [1, 20]

# So we have a problem with 60 variables: in R60 (like secret1 used to be)

# Constraints:
# ------------

# 6 contraints per Time Step
# x_k,min <= x_k <= x_k,max
# u_k,min <= u_k <= u_k,max

# 4 constraints (in inequality form => *2) at every Time Step
# x_k+1 = A * x_k + B * u_k

# 1 constraint for 1st Time Step
# x_0 = x_init

# Obstacles avoidance constraints
# 10 obstacles => 10 abs() (inconsistent) constraints

# So we are dealing with: 20*(6+4) + 1 + 10*2 constraints

# So we have around 250 constraints (twice more than secret1)

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

	MpcPath1d(T=16, dt=0.250,
			  Q=Diagonal([1.0, 50]), R=Diagonal([0.001]),
			  smin=0.0, smax=1000.0,
			  vmin=0.0, vmax=25.0, # Between 0 and 25 m.s-1
			  umin=-4.0, umax=2.0, # Between -4 and 2 m.s-2
			  dsaf=10.0,
			  Ad=[1.0 dt; 0 1],
			  Bd=[0.5*dt^2, dt],

			  # NB: xref, uref, xinit, obstacles are expected to be changed/updated everytime we plan (TODO)
			  xref=[200.0, 20.0], # Target pos=200 at v=20 m.s-1
			  uref=[0.0],
			  xinit=[0.0, 20.0], # Start at s=0 with speed v=20 m.s-1
			  obstacles=[(2.0, 100)], # In 2 sec a crossing vehicle at s=100 m

			  nvars_dt=3 # x=[s,sd] u=[sdd]
			  ) = new(T,dt,Q,R,smin,smax,vmin,vmax,umin,umax,dsaf,Ad,Bd,xref,uref,xinit,obstacles,nvars_dt)
end

global mpc = MpcPath1d()


# We define a quadratic cost function
@counted function path1d(x::Vector)
	cost = 0

	# stage cost
	for k in range(1, step=mpc.nvars_dt, length=mpc.T-1)
		xk, uk = x[k:k+1], [x[k+2]]
		cost += (xk - mpc.xref)' * mpc.Q * (xk - mpc.xref)
		cost += (uk - mpc.uref)' * mpc.R * (uk - mpc.uref)
	end

	xT, uT = x[end-2:end-1], [x[end]]
	cost += (xT - mpc.xref)' * mpc.Q * (xT - mpc.xref)
	cost = 0.5*cost # Just in case someday we want to compute gradient analytically

	println("cost=$cost")
	return cost
end

@counted function path1d_gradient(x::Vector)
	println("UNUSED sofar") # we use finite differences approx in optimize.jl anyways
	exit(1)
end

@counted function path1d_constraints(x::Vector)
	return [-1] # no constraint so far

	# TODO
	constraints = [] # They all must be of the form <= 0 ie Equality constr => 2 ineq

	return [x[1] + x[2]^2 - 1,
			-x[1] - x[2]]
end

function path1d_init()
	return rand(mpc.T * mpc.nvars_dt) * 40.0 # start with garbage values

	# TODO (start with a dynamically feasible traj)

	# define (somewhat) random obstacles
end

