#=
Julia program to find a circumscribed ellipse for given 2d point data set
=#

"""
Module for global parameters and variables
"""
module ParamVar
    struct Parameters
        num_points::Int64  # Number of particles
        x_lim::Float64  # Boundary of square region [-x_lim, x_lim] × [-x_lim, x_lim]
        semimajor_distribution::Float64  # Region of points distribution
        semiminor_distribution::Float64
        angle_distribution::Float64   # Angle of points distribution
    end

    mutable struct Variables
        points_x::Array{Float64}  # x coordinate of point group
        points_y::Array{Float64}  # y coordinate of point group
    end
end


"""
Module for computation
"""
module Compute
using Distributions

    """
    Distribute points in rectangular region & perform rotation
    """
    function distribute_points(param)
        # Distribute points in rectangular region
        x = rand(
            Uniform(-param.semimajor_distribution, param.semimajor_distribution),
            param.num_points)
        y = rand(
            Uniform(-param.semiminor_distribution, param.semiminor_distribution),
            param.num_points)

        # Rotate rectangular region
        x_rot = x .* cos(param.angle_distribution) - y .* sin(param.angle_distribution)
        y_rot = x .* sin(param.angle_distribution) + y .* cos(param.angle_distribution)

        return x_rot, y_rot
    end
end


"""
Module for plot
"""
module Output
using Plots
gr()


    """
    Plot points
    """
    function plot_points(param, var)
        p = scatter(
            var.points_x, var.points_y,
            markercolor = :black,
            aspect_ratio = 1,
            xlims = (-param.x_lim, param.x_lim),
            ylims = (-param.x_lim, param.x_lim),
            axis = nothing,
            size=(640, 640),
            title = "Distributed points"
        )

        savefig(p, "./tmp/distributed_points.png")
    end
end

# ========================================
# Main function
# ========================================

## Declare modules
using Printf
using Plots
gr(
    markerstrokewidth = 0,
    markersize = 10
)
using .ParamVar
using .Compute:
distribute_points
using .Output


# ----------------------------------------
## Declare parameters
# ----------------------------------------
num_points = 10
x_lim = 1.0
semimajor_distribution = 0.5 * x_lim
semiminor_distribution = 0.3 * x_lim
angle_distribution = 0.25 * π

param = ParamVar.Parameters(
    num_points, x_lim,
    semimajor_distribution, semiminor_distribution, angle_distribution
)




# ----------------------------------------
## Declare variables
# ----------------------------------------
points_x = Array{Float64}(undef, param.num_points)
points_y = Array{Float64}(undef, param.num_points)

var = ParamVar.Variables(points_x, points_y)


# ----------------------------------------
## Define 2d point group
# ----------------------------------------
var.points_x, var.points_y = distribute_points(param)

###CHECK###
# Plot distribution
Output.plot_points(param, var)
###CHECK###


# ----------------------------------------
## Find a circumscribed circle of point group
## Define centre & radius
## By iterative method
# ----------------------------------------



# ----------------------------------------
## Find a circumscribed ellipse of point group
## Define
## - centre
## - semimajor axis
## - semininor axis
## - angle of semimajor axis & x axis
# ----------------------------------------
