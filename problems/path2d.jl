########### "Path2d" Problems Definition ###########

# Via Optimization over time:
# ---------------------------
# We just need a set of waypoints to define a path
# A parametric representation is not required

# Handle longitudinal and lateral safety constraints
# Handle longitudinal and lateral comfort objectives

# TODO


# CurvRoad problem (sample waypoints from a 6th order polynomial)
# positions and headings providing longi and lat directions
function path2d_CurvRoad_gen()
end

@counted function path2d_CurvRoad(x::Vector)
end

@counted function path2d_CurvRoad_gradient(x::Vector)
	println("UNUSED sofar") # we use finite differences approx in optimize.jl anyways
	exit(1)
end

@counted function path2d_CurvRoad_constraints(x::Vector)
end

function path2d_CurvRoad_init()
end

# Overtake problem (sample waypoints from a 6th order polynomial)
# positions and headings providing longi and lat directions
function path2d_Overtake_gen()
end

# Intersec problem (sample waypoints from a 6th order polynomial)
# positions and headings providing longi and lat directions
function path2d_Intersec_gen()
end
