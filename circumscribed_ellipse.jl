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
    end

    struct Distribution
        semimajor::Float64  # Region of points distribution
        semiminor::Float64
        angle_x::Float64   # Angle of points distribution
        angle_y::Float64   # Angle of points distribution
        angle_z::Float64   # Angle of points distribution
    end

    """
    For point group
    """
    mutable struct Points
        x::Array{Float64}  # x coordinate
        y::Array{Float64}  # y coordinate
        z::Array{Float64}  # z coordinate
    end

    """
    For circumscribed circle
    """
    mutable struct Circle
        centre_x::Float64
        centre_y::Float64
        centre_z::Float64
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
using Printf
using Distributions

    """
    Compute distance between two points
    """
    function compute_distance(x1, y1, z1, x2, y2, z2)
        dist_square = (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2
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
    function distribute_points(param, dist, points)
        # Distribute points in rectangular region
        # x has semimajor, y & z have semiminor axis
        x = rand(
            Uniform(-dist.semimajor, dist.semimajor),
            param.num_points)
        y = rand(
            Uniform(-dist.semiminor, dist.semiminor),
            param.num_points)
        z = rand(
            Uniform(-dist.semiminor, dist.semiminor),
            param.num_points)

        # Rotate rectangular region
        x_z, y_z = compute_rotation_counterclockwise(x, y, dist.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z, dist.angle_y)  # Along y axis
        y_xz, z_xy = compute_rotation_counterclockwise(y_z, z_y, dist.angle_x)  # Along x axis

        points.x = x_yz
        points.y = y_xz
        points.z = z_xy
    end


    """
    Find centre of circumscribed circle of given point group
    cf. https://www.ipsj.or.jp/07editj/promenade/4309.pdf p.6-7
    """
    function search_circumscribed_circle(param, points, circle)
        # Initial guess = corner of region
        centre_x = -param.x_lim
        centre_y = -param.x_lim
        centre_z = -param.x_lim

        # Ratio of movement/distance to the most far point
        num_move = param.num_points
        move = 0.5 * param.x_lim
        dist_max = 0.0

        while move >= 1.0e-4
            for itr_move = 1:num_move
                # Find the most far point
                dist_max = 0.0
                object_x = 0.0
                object_y = 0.0
                object_z = 0.0
                for itr_point = 1:param.num_points
                    dist = compute_distance(
                        centre_x, centre_y, centre_z,
                        points.x[itr_point], points.y[itr_point], points.z[itr_point]
                    )
                    if dist > dist_max
                        dist_max = dist
                        object_x = points.x[itr_point]
                        object_y = points.y[itr_point]
                        object_z = points.z[itr_point]
                    end
                end

                # Move towards the most far point
                centre_x += move * (object_x - centre_x)
                centre_y += move * (object_y - centre_y)
                centre_z += move * (object_z - centre_z)
            end

            # Dicrease move ratio
            move *= 0.5
        end

        circle.centre_x = centre_x
        circle.centre_y = centre_y
        circle.centre_z = centre_z
        circle.radius = dist_max
    end


    """
    Find circumscribed ellipse of given point group based on circumscribed circle
    1. Find the most distant point
    2. Define semimajor axis length & angle as origin <-> the most distant point
    3. Shift by the centre coordinate
    4. Apply inverse rotation of point group by semimajor axis angle
    5. Adjust semimajor axis: x -> x/a
    -----iteration for semiminor axis length (originally circle radius)
    6. Adjust semiminor axis: y -> y/b
    7. Confirm all points are included in the standard circle of radius=1
        If so, decrease the semiminor axis length by delta
        If not, exit the loop & define the semiminor axis length
    -----end iteration for semiminor axis
    """
    function search_circumscribed_ellipse(param, points, circle, ellipse)
        # 1. Find the most distant point
        dist_max = 0.0
        distant_x = 0.0
        distant_y = 0.0
        for itr_point = 1:param.num_points
            dist = compute_distance(
                circle.centre_x, circle.centre_y,
                points.x[itr_point], points.y[itr_point]
            )
            if dist > dist_max
                dist_max = dist
                distant_x = points.x[itr_point]
                distant_y = points.y[itr_point]
            end
        end

        # 2. Define semimajor axis length & angle
        semimajor_length = compute_distance(
            circle.centre_x, circle.centre_y,
            distant_x, distant_y
        )
        semimajor_angle = atan(distant_y - circle.centre_y, distant_x - circle.centre_x)

        # 3. Shift by the centre coordinate
        x_shift = points.x .- circle.centre_x
        y_shift = points.y .- circle.centre_y

        # 4. Apply inverse rotation of point group by semimajor axis angle
        x_rot, y_rot = compute_rotation_clockwise(x_shift, y_shift, semimajor_angle)

        # 5. Adjust semimajor axis: x -> x/a
        x_std = x_rot / semimajor_length

        # parameter setting for while-loop
        semiminor_length = circle.radius
        flag_all_points_in_circle = true
        delta = 0.01

        while flag_all_points_in_circle

            # 6. Adjust semiminor axis: y -> y/b
            y_std = y_rot / semiminor_length

            # 7. Confirm all points are included in the circle
            for itr_point in 1:param.num_points
                dist = compute_distance(
                    0.0, 0.0,
                    x_std[itr_point], y_std[itr_point]
                )

                # 6-2. If not, set flag to exit the loop
                if dist > 1
                    flag_all_points_in_circle = false
                    break
                end
            end

            # 6-1. If so, decrease the semiminor axis length
            if flag_all_points_in_circle
                semiminor_length -= delta
            end

            ###CHECK###
            s = @sprintf("semimajor %.3f, semiminor %.3f", semimajor_length, semiminor_length)
            println(s)
            ###CHECK###

        end

        ellipse.semimajor = semimajor_length
        ellipse.semiminor = semiminor_length
        ellipse.angle = semimajor_angle
    end
end


"""
Module for plot
"""
module Output
using Printf
using Plots
gr()

    """
    Plot points
    """
    function plot_points(param, points)
        p = plot(
            points.x, points.y, points.z,
            seriestype=:scatter,
            markercolor = :black,
            aspect_ratio = 1,
            xlims = (-param.x_lim, param.x_lim),
            ylims = (-param.x_lim, param.x_lim),
            zlims = (-param.x_lim, param.x_lim),
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
    Define ellipse shape
    1. Circle of radius=1: x^2+y^2=1
    2. Adjust semimajor/minor axis: x -> ax, y -> by
    3. Rotation by angle
    4. Shift the centre of ellipse
    """
    function ellipse_shape(centre_x, centre_y, semimajor, semiminor, angle)
        θ = LinRange(0, 2*π, 100)
        # 2
        x_tmp = semimajor * sin.(θ)
        y_tmp = semiminor * cos.(θ)
        # 3
        x_rot = x_tmp .* cos(angle) - y_tmp .* sin(angle)
        y_rot = x_tmp .* sin(angle) + y_tmp .* cos(angle)
        # 4
        x_rot .+ centre_x, y_rot .+ centre_y
    end


    """
    Define line shape
    """
    function line_shape(param, fixed_x, fixed_y, angle)
        x = -param.x_lim:0.01:param.x_lim
        intercept = fixed_y - fixed_x * tan(angle)
        x, tan(angle) .* x .+ intercept
    end


    """
    Plot points, circumscribed circle & circumscribed ellipse
    """
    function plot_points_circumscribed(param, points, circle, ellipse)
        s = @sprintf("ellipse's semimajor axis %.3f, semiminor axis %.3f & angle %.3f", ellipse.semimajor, ellipse.semiminor, ellipse.angle)
        println(s)

        # Point group
        p = scatter(
            points.x, points.y,
            markercolor = :black,
            aspect_ratio = 1,
            xlims = (-param.x_lim, param.x_lim),
            ylims = (-param.x_lim, param.x_lim),
            axis = nothing,
            size=(640, 640),
            title="Point group & its circumscribe"
        )

        # Centre of circle/ellipse
        p! = scatter!(
            [circle.centre_x], [circle.centre_y],
            markercolor = :red
        )

        # Circumscribed circle
        p! = plot!(
            circle_shape(circle.centre_x, circle.centre_y, circle.radius),
            seriestype = [:shape,],
            lw = 0.5,
            c = :blue,
            linecolor = :black,
            legend = false,
            fillalpha =0.2,
            aspect_ratio = 1
        )

        # Circumscribed ellipse
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
            aspect_ratio = 1
        )

        # Semimajor axis of circumscribed ellipse
        p! = plot!(
            line_shape(
                param,
                circle.centre_x, circle.centre_y, ellipse.angle
            ),
            lw = 2.0,
            linecolor = :red,
            legend = false
        )

        savefig(p!, "./tmp/circumscribed.png")
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
plot_points_circumscribed


# ----------------------------------------
## Declare parameters
# ----------------------------------------
num_points = 10
x_lim = 1.0

param = ParamVar.Parameters(
    num_points, x_lim
    )


# ----------------------------------------
## Define 3d point group
# ----------------------------------------
semimajor = 0.6 * x_lim
semiminor = 0.2 * x_lim
angle_x = 0.2 * π  # Rotation along x axis
angle_y = 0.4 * π  # y axis
angle_z = 0.6 * π  # z axis

dist = ParamVar.Distribution(
    semimajor, semiminor,
    angle_x, angle_y, angle_z
)

x = Array{Float64}(undef, param.num_points)
y = Array{Float64}(undef, param.num_points)
z = Array{Float64}(undef, param.num_points)
points = ParamVar.Points(x, y, z)

distribute_points(param, dist, points)

###CHECK###
# Plot distribution
plot_points(param, points)
###CHECK###


# ----------------------------------------
## Find a circumscribed circle of point group
## Define centre & radius
## By iterative method
# ----------------------------------------
centre_x = 0.0
centre_y = 0.0
centre_z = 0.0
radius = 0.0
circle = ParamVar.Circle(centre_x, centre_y, centre_z, radius)

search_circumscribed_circle(param, points, circle)


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
# # set ellipse for check plot
# ellipse.semimajor = param.semimajor_distribution
# ellipse.semiminor = param.semiminor_distribution
# ellipse.angle = param.angle_distribution
###CHECK###

# Plot points, circumscribed circle & circumscribed ellipse
plot_points_circumscribed(param, points, circle, ellipse)
