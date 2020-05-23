########### "Path1d" Problem Definition ###########

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


# TODO ...


global T = 20 # Time steps


@counted function path1d(x::Vector)
	println(T)
	return -x[1] * x[2] + 2.0 / (3.0 * sqrt(3.0))
end

@counted function path1d_gradient(x::Vector)
	println("UNUSED sofar") # we use finite differences approx in optimize.jl anyways
	exit(1)
end

@counted function path1d_constraints(x::Vector)
	return [x[1] + x[2]^2 - 1,
			-x[1] - x[2]]
end

function path1d_init()
	return rand(2) * 2.0
end

