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

    """
    For point group
    """
    mutable struct Points
        x::Array{Float64}  # x coordinate
        y::Array{Float64}  # y coordinate
    end

    """
    For circumscribed circle
    """
    mutable struct Circle
        centre_x::Float64
        centre_y::Float64
        radius::Float64
    end

    """
    For circumscribed ellipse
    """
    mutable struct Ellipse
        semimajor::Float64
        semiminor::Float64
        angle::Float64
    end
end


"""
Module for computation
"""
module Compute
using Distributions

    """
    Compute distance between two points
    """
    function compute_distance(x1, y1, x2, y2)
        dist_square = (x1-x2)^2 + (y1-y2)^2
        return sqrt(dist_square)
    end


    """
    Compute rotation by counter-clockwise
    """
    function compute_rotation_counterclockwise(x, y, θ)
        x_new = x .* cos(θ) - y .* sin(θ)
        y_new = x .* sin(θ) + y .* cos(θ)

        return x_new, y_new
    end


    """
    Compute rotation by clockwise
    """
    function compute_rotation_clockwise(x, y, θ)
        x_new =  x .* cos(θ) + y .* sin(θ)
        y_new = -x .* sin(θ) + y .* cos(θ)

        return x_new, y_new
    end


    """
    Distribute points in rectangular region & perform rotation
    """
    function distribute_points(param, points)
        # Distribute points in rectangular region
        x = rand(
            Uniform(-param.semimajor_distribution, param.semimajor_distribution),
            param.num_points)
        y = rand(
            Uniform(-param.semiminor_distribution, param.semiminor_distribution),
            param.num_points)

        # Rotate rectangular region
        x_rot, y_rot = compute_rotation_counterclockwise(x, y, param.angle_distribution)

        points.x = x_rot
        points.y = y_rot
    end


    """
    Find centre of circumscribed circle of given point group
    cf. https://www.ipsj.or.jp/07editj/promenade/4309.pdf p.6-7
    """
    function search_circumscribed_circle(param, points, circle)
        # Initial guess = corner of region
        centre_x = -param.x_lim
        centre_y = -param.x_lim

        # Ratio of movement/distance to the most far point
        num_move = 100
        move = 0.5 * param.x_lim
        dist_max = 0.0

        while move >= 1.0e-4
            for itr_move = 1:num_move
                # Find the most far point
                dist_max = 0.0
                object_x = 0.0
                object_y = 0.0
                for itr_point = 1:param.num_points
                    dist = compute_distance(
                        centre_x, centre_y,
                        points.x[itr_point], points.y[itr_point]
                    )
                    if dist > dist_max
                        dist_max = dist
                        object_x = points.x[itr_point]
                        object_y = points.y[itr_point]
                    end
                end

                # Move towards the most far point
                centre_x += move * (object_x - centre_x)
                centre_y += move * (object_y - centre_y)
            end

            # Dicrease move ratio
            move *= 0.5
        end

        circle.centre_x = centre_x
        circle.centre_y = centre_y
        circle.radius = dist_max
    end


    """
    Find circumscribed ellipse of given point group based on circumscribed circle
    """
    function search_circumscribed_ellipse(param, points, circle, ellipse)

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
    function plot_points(param, points)
        p = scatter(
            points.x, points.y,
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


    """
    Define circle shape
    """
    function circle_shape(centre_x, centre_y, radius)
        θ = LinRange(0, 2*π, 100)
        centre_x .+ radius*sin.(θ), centre_y .+ radius*cos.(θ)
    end


    """
    Plot points & circumscribed circle
    """
    function plot_points_circle(param, points, circle)
        p = scatter(
            points.x, points.y,
            markercolor = :black,
            aspect_ratio = 1,
            xlims = (-param.x_lim, param.x_lim),
            ylims = (-param.x_lim, param.x_lim),
            axis = nothing,
            size=(640, 640)
        )
        p! = scatter!(
            [circle.centre_x], [circle.centre_y],
            markercolor = :red
        )
        p! = plot!(
            circle_shape(circle.centre_x, circle.centre_y, circle.radius),
            seriestype = [:shape,],
            lw = 0.5,
            c = :blue,
            linecolor = :black,
            legend = false,
            fillalpha =0.2,
            aspect_ratio = 1,
            title = "circumscribed circle"
        )

        savefig(p!, "./tmp/circumscribed_circle.png")
    end


    """
    Define ellipse shape
    1. Circle of radius=1: x^2+y^2=1
    2. Adjust semimajor/minor axis: x^2 -> x^2/a^2, y^2 -> y^2/b^2
    3. Rotation by angle
    """
    function ellipse_shape(centre_x, centre_y, semimajor, semiminor, angle)
        θ = LinRange(0, 2*π, 100)
        # 1 & 2
        x_tmp = semimajor * (centre_x .+ sin.(θ))
        y_tmp = semiminor * (centre_y .+ cos.(θ))
        # 3
        x_tmp .* cos(angle) - y_tmp .* sin(angle),
        x_tmp .* sin(angle) + y_tmp .* cos(angle)
    end


    """
    Plot points & circumscribed ellipse
    """
    function plot_points_ellipse(param, points, circle, ellipse)
        p = scatter(
            points.x, points.y,
            markercolor = :black,
            aspect_ratio = 1,
            xlims = (-param.x_lim, param.x_lim),
            ylims = (-param.x_lim, param.x_lim),
            axis = nothing,
            size=(640, 640)
        )
        p! = scatter!(
            [circle.centre_x], [circle.centre_y],
            markercolor = :red
        )
        p! = plot!(
            ellipse_shape(
                circle.centre_x, circle.centre_y,  # Centre is same as circle
                ellipse.semimajor, ellipse.semiminor, ellipse.angle),
            seriestype = [:shape,],
            lw = 0.5,
            c = :green,
            linecolor = :black,
            legend = false,
            fillalpha =0.2,
            aspect_ratio = 1,
            title = "circumscribed ellipse"
        )

        savefig(p!, "./tmp/circumscribed_ellipse.png")
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
distribute_points,
search_circumscribed_circle,
search_circumscribed_ellipse
using .Output:
plot_points,
plot_points_circle,
plot_points_ellipse


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
## Define 2d point group
# ----------------------------------------
x = Array{Float64}(undef, param.num_points)
y = Array{Float64}(undef, param.num_points)
points = ParamVar.Points(x, y)

distribute_points(param, points)

###CHECK###
# Plot distribution
# plot_points(param, points)
###CHECK###


# ----------------------------------------
## Find a circumscribed circle of point group
## Define centre & radius
## By iterative method
# ----------------------------------------
centre_x = 0.0
centre_y = 0.0
radius = 0.0
circle = ParamVar.Circle(centre_x, centre_y, radius)

search_circumscribed_circle(param, points, circle)

###CHECK###
# Plot distribution
plot_points_circle(param, points, circle)
###CHECK###


# ----------------------------------------
## Find a circumscribed ellipse of point group
## Define
## - centre
## - semimajor & semininor axis
## - angle of semimajor axis & x axis
# ----------------------------------------
semimajor = 0.0
semiminor = 0.0
angle = 0.0
ellipse = ParamVar.Ellipse(semimajor, semiminor, angle)

search_circumscribed_ellipse(param, points, circle, ellipse)

###CHECK###
# set ellipse for check plot
ellipse.semimajor = param.semimajor_distribution
ellipse.semiminor = param.semiminor_distribution
ellipse.angle = param.angle_distribution
###CHECK###

# Plot points & circumscribed ellipse
plot_points_ellipse(param, points, circle, ellipse)
