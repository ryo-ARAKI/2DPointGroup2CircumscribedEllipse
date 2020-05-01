#=
Julia program to find a circumscribed ellipse for given 3d point data set
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
        angle_y::Float64
        angle_z::Float64
        shift_x::Float64   # Shift of points distribution
        shift_y::Float64
        shift_z::Float64
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
        angle_x::Float64
        angle_y::Float64
        angle_z::Float64
    end
end


"""
Module for rotation computation
"""
module ComputeRotation

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
end


"""
Module for define various shapes
"""
module DefineShape
    using ..ComputeRotation:
        compute_rotation_counterclockwise,
        compute_rotation_clockwise

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
        x_rot, y_rot = compute_rotation_counterclockwise(x_tmp, y_tmp, angle)
        # 4
        x_rot .+ centre_x, y_rot .+ centre_y
    end


    """
    Define sphere shape
    1. Sphere of radius=1: x^2+y^2+z^2=1
    2. Multiply radius
    2. Shift the centre of shpere
    """
    function sphere_shape(circle)
        #1
        x_tmp, y_tmp, z_tmp = compute_sphere()

        # 2
        x_tmp *= circle.radius
        y_tmp *= circle.radius
        z_tmp *= circle.radius

        # 3
        circle.centre_x .+ x_tmp, circle.centre_y .+ y_tmp, circle.centre_z .+ z_tmp
    end


    """
    Compute standard sphere
    Origin = zero, radius=1
    """
    function compute_sphere()
        dim = 30
        θ = range(0, stop=π, length=dim)
        ϕ = range(0, stop=2*π, length=dim)

        x_tmp = sin.(θ) * cos.(ϕ)'
        y_tmp = sin.(θ) * sin.(ϕ)'
        z_tmp = cos.(θ) * ones(dim)'

        return x_tmp, y_tmp, z_tmp
    end


    """
    Define spheroid shape
    1. Sphere of radius=1: x^2+y^2+z^2=1
    2. Adjust semimajor/minor axis: x -> ax, y -> by, z -> bz
    3. Rotation by angle
    4. Shift the centre of ellipse
    """
    function spheroid_shape(circle, ellipse)
        # 1
        x_tmp, y_tmp, z_tmp = compute_sphere()

        # 2
        x_tmp *= ellipse.semimajor
        y_tmp *= ellipse.semiminor
        z_tmp *= ellipse.semiminor

        # 3
        x_z, y_z = compute_rotation_counterclockwise(x_tmp, y_tmp, ellipse.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z_tmp, ellipse.angle_y)  # Along y axis
        y_xz, z_xy = compute_rotation_counterclockwise(y_z, z_y, ellipse.angle_x)  # Along x axis

        # 4
        circle.centre_x .+ x_yz, circle.centre_y .+ y_xz, circle.centre_z .+ z_xy

    end


    """
    Define line shape
    """
    function line_shape(param, fixed_x, fixed_y, angle)
        x = -param.x_lim:0.01:param.x_lim
        intercept = fixed_y - fixed_x * tan(angle)
        x, tan(angle) .* x .+ intercept
    end
end


"""
Module for computation
"""
module Compute
    using ..ComputeRotation:
        compute_rotation_counterclockwise,
        compute_rotation_clockwise
    using ..DefineShape:
        sphere_shape
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
    Distribute points
    1. Distribute points randomly in whole region
    2. Define sphere of radius=1: x^2+y^2+z^2=1
    3. Delete points outside of sphere
    4. Adjust semimajor/minor axis: x -> ax, y -> by, z -> bz
    5. Rotation by angle
    6. Shift the centre of ellipse
    """
    function distribute_points(param, dist, points)
        # 1. Distribute points randomly in whole region
        # -----Prepare 3*num_points points, since points outside of sphere will be deleted
        x = rand(
            Uniform(-param.x_lim, param.x_lim),
            3*param.num_points
        )
        y = rand(
            Uniform(-param.x_lim, param.x_lim),
            3*param.num_points
        )
        z = rand(
            Uniform(-param.x_lim, param.x_lim),
            3*param.num_points
        )

        # 2. Define sphere of radius=1: x^2+y^2+z^2=1
        # 3. Delete points outside of shepre
        for itr_point in 1:3*param.num_points
            distance = compute_distance(
                0.0, 0.0, 0.0,
                x[itr_point], y[itr_point], z[itr_point]
            )
            if distance > param.x_lim
                for itr_overwrite in itr_point:3*param.num_points-1
                    x[itr_overwrite] = x[itr_overwrite+1]  # Overwrite points outside of sphere
                    y[itr_overwrite] = y[itr_overwrite+1]
                    z[itr_overwrite] = z[itr_overwrite+1]

                    x[3*param.num_points] = 0.0
                    y[3*param.num_points] = 0.0
                    z[3*param.num_points] = 0.0
                end
            end
        end

        # 3.5 Pick up num_points points from array
        x_sphere = x[1:param.num_points]
        y_sphere = y[1:param.num_points]
        z_sphere = z[1:param.num_points]

        #=
        They might not have enough number of points in shpere...
        Need to reconsider better implementation
        =#

        # 4. Adjust semimajor/minor axis: x -> ax, y -> by, z -> bz
        x_sphere *= dist.semimajor
        y_sphere *= dist.semiminor
        z_sphere *= dist.semiminor

        # Rotate rectangular region
        x_z, y_z = compute_rotation_counterclockwise(x_sphere, y_sphere, dist.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z_sphere, dist.angle_y)  # Along y axis
        y_xz, z_xy = compute_rotation_counterclockwise(y_z, z_y, dist.angle_x)  # Along x axis

        # Shift rectangular region
        x_shift = x_yz .+ dist.shift_x
        y_shift = y_xz .+ dist.shift_y
        z_shift = z_xy .+ dist.shift_z

        points.x = x_shift
        points.y = y_shift
        points.z = z_shift
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
    1. Shift by the centre coordinate
    2. Find the most distant point
    3. Define semimajor axis length & angle as origin <-> the most distant point
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
        # 1. Shift by the centre coordinate
        x_shift = points.x .- circle.centre_x
        y_shift = points.y .- circle.centre_y
        z_shift = points.z .- circle.centre_z

        # 2. Find the most distant point
        dist_max = 0.0
        distant_x = 0.0
        distant_y = 0.0
        distant_z = 0.0
        for itr_point = 1:param.num_points
            dist = compute_distance(
                0.0, 0.0, 0.0,
                x_shift[itr_point], y_shift[itr_point], z_shift[itr_point]
            )
            if dist > dist_max
                dist_max = dist
                distant_x = x_shift[itr_point]
                distant_y = y_shift[itr_point]
                distant_z = z_shift[itr_point]
            end
        end

        # 3. Define semimajor axis length & angle
        semimajor_length = compute_distance(
            0.0, 0.0, 0.0,
            distant_x, distant_y, distant_z
        )
        semimajor_angle_x = atan(distant_z, distant_y)  # Along x axis
        semimajor_angle_y = atan(distant_z, distant_x)  # Along y axis
        semimajor_angle_z = atan(distant_y, distant_x)  # Along z axis


        # 4. Apply inverse rotation of point group by semimajor axis angle
        y_z, z_y = compute_rotation_clockwise(y_shift, z_shift, semimajor_angle_x)  # Along x axis
        x_z, z = compute_rotation_clockwise(x_shift, z_y, semimajor_angle_y)  # Along y axis
        x, y = compute_rotation_clockwise(x_z, y_z, semimajor_angle_z)  # Along z axis

        # 5. Adjust semimajor axis: x -> x/a
        x_std = x / semimajor_length

        # parameter setting for while-loop
        semiminor_length = circle.radius
        flag_all_points_in_circle = true
        delta = 0.01

        while flag_all_points_in_circle

            # 6. Adjust semiminor axis: y -> y/b & z -> z/b
            y_std = y / semiminor_length
            z_std = z / semiminor_length

            # 7. Confirm all points are included in the circle
            for itr_point in 1:param.num_points
                dist = compute_distance(
                    0.0, 0.0, 0.0,
                    x_std[itr_point], y_std[itr_point], z_std[itr_point]
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

        end

        println("\nResult:")
        println(
            @sprintf "semimajor %.3f, semiminor %.3f, angle %.3f %.3f %.3f" semimajor_length semiminor_length semimajor_angle_x semimajor_angle_y semimajor_angle_z
        )

        ellipse.semimajor = semimajor_length
        ellipse.semiminor = semiminor_length
        ellipse.angle_x = semimajor_angle_x
        ellipse.angle_y = semimajor_angle_y
        ellipse.angle_z = semimajor_angle_z
    end
end


"""
Module for plot
"""
module Output
    using ..DefineShape:
        sphere_shape,
        circle_shape,
        ellipse_shape,
        sphere_shape,
        spheroid_shape,
        compute_sphere,
        line_shape
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
    Plot points, circumscribed circle & circumscribed ellipse
    """
    function plot_points_circumscribed(param, points, circle, ellipse)

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


    """
    Output points, circumscribed circle & circumscribed ellipse in dat file
    """
    function out_points_circumscribed(param, points, circle, ellipse)

        # Points information
        pointsfile = open("./tmp/points.dat","w")
        for itr_point = 1:param.num_points
            write(pointsfile, "$(points.x[itr_point])\t$(points.y[itr_point])\t$(points.z[itr_point])\n")
        end
        close(pointsfile)

        # Sphere information
        sphere = sphere_shape(circle)
        dims_sphere = size.(sphere)
        lens_sphere = length.(sphere)

        spherefile = open("./tmp/sphere.dat","w")
        for itr_point = 1:lens_sphere[1]
            write(spherefile, "$(sphere[1][itr_point])\t$(sphere[2][itr_point])\t$(sphere[3][itr_point])\n")
        end
        close(spherefile)

        # Spheroid information
        spheroid = spheroid_shape(circle, ellipse)
        dims_spheroid = size.(spheroid)
        lens_spheroid = length.(spheroid)

        spheroidfile = open("./tmp/spheroid.dat","w")
        for itr_point = 1:lens_spheroid[1]
            write(spheroidfile, "$(spheroid[1][itr_point])\t$(spheroid[2][itr_point])\t$(spheroid[3][itr_point])\n")
        end
        close(spheroidfile)

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
    plot_points_circumscribed,
    out_points_circumscribed



# ----------------------------------------
## Declare parameters
# ----------------------------------------
num_points = 50
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
angle_y = 0.4 * π
angle_z = 0.6 * π
shift_x = 0.1  # Shift along x axis
shift_y = 0.2
shift_z = 0.3

dist = ParamVar.Distribution(
    semimajor, semiminor,
    angle_x, angle_y, angle_z,
    shift_x, shift_y, shift_z
)

println("Given data:")
println(
    @sprintf "semimajor %.3f, semiminor %.3f, angle %.3f %.3f %.3f" dist.semimajor dist.semiminor dist.angle_x dist.angle_y dist.angle_z
)

x = Array{Float64}(undef, param.num_points)
y = Array{Float64}(undef, param.num_points)
z = Array{Float64}(undef, param.num_points)
points = ParamVar.Points(x, y, z)

distribute_points(param, dist, points)

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
angle_x = 0.0
angle_y = 0.0
angle_z = 0.0
ellipse = ParamVar.Ellipse(
    semimajor, semiminor,
    angle_x, angle_y, angle_z
)

search_circumscribed_ellipse(param, points, circle, ellipse)

###CHECK###
# set ellipse for check plot
# ellipse.semimajor = dist.semimajor
# ellipse.semiminor = dist.semiminor
# ellipse.angle_x = dist.angle_x
# ellipse.angle_y = dist.angle_y
# ellipse.angle_z = dist.angle_z
###CHECK###

"""
# Plot result in 2D
# Plot points, circumscribed circle & circumscribed ellipse
plot_points_circumscribed(param, points, circle, ellipse)
"""

# Output result
out_points_circumscribed(param, points, circle, ellipse)
