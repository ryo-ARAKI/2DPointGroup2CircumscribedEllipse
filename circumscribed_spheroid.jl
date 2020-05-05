#=
Julia program to find a circumscribed spheroid for given 3d point data set
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
        angle_y::Float64
        angle_z::Float64
        shift_x::Float64   # Shift of points distribution
        shift_y::Float64
        shift_z::Float64
    end

    """
    Check parameters
    """
    function ParameterCheck(dist)
        # Check angles
        if dist.angle_y<-π/2 || dist.angle_y>π/2
            throw(DomainError(dist.angle_y, "angle along y axis must be [-π/2, π/2]"))
        end
        if dist.angle_z<-π/2 || dist.angle_z>π/2
            throw(DomainError(dist.angle_z, "angle along z axis must be [-π/2, π/2]"))
        end
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
    For circumscribed sphere
    """
    mutable struct Sphere
        centre_x::Float64
        centre_y::Float64
        centre_z::Float64
        radius::Float64
    end

    """
    For circumscribed spheroid
    """
    mutable struct Spheroid
        semimajor::Float64
        semiminor::Float64
        angle_y::Float64
        angle_z::Float64
    end
end


"""
Module for rotation computation
"""
module ComputeRotation
    using Printf

    """
    Compute rotation by counter-clockwise
    """
    function compute_rotation_counterclockwise(x, y, θ)
        x_new = x .* cos(θ) - y .* sin(θ)
        y_new = x .* sin(θ) + y .* cos(θ)

        return x_new, y_new
    end


    """
    DEBUG
    Check rotation manipulation
    """
    function check_rotation()
        # Define test vector
        x = 0.6; y = 0.0; z = 0.0

        # Perform rotation
        β = π/3.0; γ = π/4.0;
        x_z, y_z = compute_rotation_counterclockwise(x, y, γ)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z, β)  # Along y axis

        semimajor = sqrt(x^2+y^2+z^2)

        # Compute angle of rotation from converted coordinate
        β_ = atan(z_y/semimajor, x_yz/semimajor)  # Along y axis
        γ_ = asin(y_z/semimajor)  # Along z axis

        # Perform inverse rotation
        x_z, z_ = compute_rotation_counterclockwise(x_yz, z_y, -β)  # Along y axis
        x_, y_ = compute_rotation_counterclockwise(x_z, y_z, -γ)  # Along z axis

        println(@sprintf "\nvec = %.3f, %.3f, %.3f (original)" x y z)
        println(@sprintf "vec = %.3f, %.3f, %.3f (rotation)" x_yz y_z z_y)
        println(@sprintf "vec = %.3f, %.3f, %.3f (rotation & inverse rotation)" x_ y_ z_)
        println(@sprintf "angle along y axis: %.3f, z axis: %.3f (parameter)" β γ)
        println(@sprintf "angle along y axis: %.3f, z axis: %.3f (computed)" β_ γ_)
    end
end


"""
Module for define various shapes
"""
module DefineShape
    using ..ComputeRotation:
        compute_rotation_counterclockwise

    """
    Define sphere shape
    """
    function sphere_shape(centre_x, centre_y, radius)
        θ = LinRange(0, 2*π, 100)
        centre_x .+ radius*sin.(θ), centre_y .+ radius*cos.(θ)
    end


    """
    Define spheroid shape
    1. Define sphere of radius=1: x^2+y^2=1
    2. Adjust semimajor/minor axis: x -> ax, y -> by
    3. Rotation by angle
    4. Shift the centre of spheroid
    """
    function spheroid_shape(centre_x, centre_y, semimajor, semiminor, angle)
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
    function sphere_shape(sphere)
        #1
        x_tmp, y_tmp, z_tmp = compute_sphere()

        # 2
        x_tmp *= sphere.radius
        y_tmp *= sphere.radius
        z_tmp *= sphere.radius

        # 3
        sphere.centre_x .+ x_tmp, sphere.centre_y .+ y_tmp, sphere.centre_z .+ z_tmp
    end


    """
    Compute standard sphere
    Origin = zero, radius=1
    """
    function compute_sphere()
        dim = 20
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
    4. Shift the centre of spheroid
    """
    function spheroid_shape(sphere, spheroid)
        # 1
        x_tmp, y_tmp, z_tmp = compute_sphere()

        # 2
        x_tmp *= spheroid.semimajor
        y_tmp *= spheroid.semiminor
        z_tmp *= spheroid.semiminor

        # 3
        x_z, y_z = compute_rotation_counterclockwise(x_tmp, y_tmp, spheroid.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z_tmp, spheroid.angle_y)  # Along y axis

        # 4
        sphere.centre_x .+ y_z, sphere.centre_y .+ x_yz, sphere.centre_z .+ z_y

    end


    """
    Define line shape
    1. Define line x = t, y=0, z=0
    2. Rotate them by predetermined angle
    3. Shift line by centre of spheroid
    """
    function line_shape(param, dist, sphere, spheroid)


        #-----Use sphere & spheroid module information-----
        # 1. Define line x = t, y=0, z=0
        len_line = Int(param.x_lim*100)
        x = range(-param.x_lim, stop=param.x_lim, length=len_line)
        y = range(0.0, stop=0.0, length=len_line)
        z = range(0.0, stop=0.0, length=len_line)

        # 2. Rotate line by angle
        x_z, y_z = compute_rotation_counterclockwise(x, y, spheroid.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z, spheroid.angle_y)  # Along y axis

        # 3. Shift line by centre of spheroid
        x_shift = x_yz .+ sphere.centre_x
        y_shift = y_z .+ sphere.centre_y
        z_shift = z_y .+ sphere.centre_z

        #=
        #-----Use fixed points for debug-----
        # 1. Define points
        x = [-1.0, 1.0]; y = [0.0, 0.0]; z = [0.0, 0.0]

        # 2. Rotate line by angle
        x_z, y_z = compute_rotation_counterclockwise(x, y, dist.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z, dist.angle_y)  # Along y axis

        # 3. Shift line by centre of distribution
        x_shift = x_yz .+ dist.shift_x
        y_shift = y_z .+ dist.shift_y
        z_shift = z_y .+ dist.shift_z
        #-----Use fixed points for debug-----
        =#

        x_shift, y_shift, z_shift
    end
end


"""
Module for computation
"""
module Compute
    using ..ComputeRotation:
        compute_rotation_counterclockwise
    using ..DefineShape:
        sphere_shape,
        compute_sphere
    using Printf
    using Distributions

    """
    DEBUG
    Check computation of angle by coordinates
    """
    function compute_angle(param, dist, points)
        # 1. Shift by the centre coordinate
        x_shift = points.x .- dist.shift_x
        y_shift = points.y .- dist.shift_y
        z_shift = points.z .- dist.shift_z

        # 2. Find the most distant point
        dist_max = 0.0

        distant_x, distant_y, distant_z = [0.0 for _ = 1:3]
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
        semimajor_angle_y = atan(distant_z/semimajor_length, distant_x/semimajor_length)  # Along y axis
        semimajor_angle_z = asin(distant_y/semimajor_length)  # Along z axis

        println(@sprintf "\nComputation from distributed points information")
        println(@sprintf "most distant point x: %.3f, y: %.3f, z: %.3f r: %.3f" distant_x distant_y distant_z semimajor_length)
        println(@sprintf "angle along y axis: %.3f, z axis: %.3f" semimajor_angle_y semimajor_angle_z)
    end

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
    6. Shift the centre of sphere
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

                    x[3*param.num_points], y[3*param.num_points], z[3*param.num_points] = [0.0 for _ = 1:3]
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

        ###CHECK###
        # Use surface points instead of points distributed in the domain
        # x_sphere, y_sphere, z_sphere = compute_sphere()
        ###CHECK###

        # 4. Adjust semimajor/minor axis: x -> ax, y -> by, z -> bz
        x_sphere *= dist.semimajor
        y_sphere *= dist.semiminor
        z_sphere *= dist.semiminor

        # 5. Rotate rectangular region
        x_z, y_z = compute_rotation_counterclockwise(x_sphere, y_sphere, dist.angle_z)  # Along z axis
        x_yz, z_y = compute_rotation_counterclockwise(x_z, z_sphere, dist.angle_y)  # Along y axis

        # 6. Shift rectangular region
        points.x = x_yz .+ dist.shift_x
        points.y = y_z .+ dist.shift_y
        points.z = z_y .+ dist.shift_z

        ###CHECK###
        # Compute angle by spheroid surface information
        # compute_angle(param, dist, points)
        ###CHECK###
    end


    """
    Find centre of circumscribed sphere of given point group
    cf. https://www.ipsj.or.jp/07editj/promenade/4309.pdf p.6-7
    """
    function search_circumscribed_sphere(param, points, sphere)
        # Initial guess = corner of region
        centre_x, centre_y, centre_z = [0.0 for _ = 1:3]

        # Ratio of movement/distance to the most far point
        num_move = param.num_points
        move = 0.5 * param.x_lim
        dist_max = 0.0

        while move >= 1.0e-12
            for itr_move = 1:num_move
                # Find the most far point
                dist_max = 0.0
                object_x, object_y, object_z = [0.0 for _ = 1:3]
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

        println("\nSphere Result:")
        println(@sprintf "radius %.3f" dist_max)
        println(@sprintf "centre x: %.3f, y: %.3f, z: %.3f" centre_x centre_y centre_z)

        sphere.centre_x = centre_x
        sphere.centre_y = centre_y
        sphere.centre_z = centre_z
        sphere.radius = dist_max
    end


    """
    Find the top **% most distant points from zero-centred 3d coordinates
    """
    function search_most_distant_points(param, x, y, z)
        # Set threshold of top **% most distant points
        thresh_distant = 0.10

        # Compute distance for all particles
        distant = zeros(Float64, param.num_points)
        for itr_point = 1:param.num_points
            distant[itr_point] = compute_distance(
                0.0, 0.0, 0.0,
                x[itr_point], y[itr_point], z[itr_point]
            )
        end

        num_distant = Array{Int64}(undef, Int(thresh_distant * param.num_points))
        for itr_distant = 1:Int(thresh_distant * param.num_points)
            # Pick up the most distant point
            index_max = argmax(distant)
            num_distant[itr_distant] = index_max

            # Zeroing index corresponding to the most distant point
            distant[index_max] = 0.0
        end

        return num_distant
    end


    """
    Compute distance as semimajor axis length & its rotation angle
    """
    function compute_semimajor_axis_angle(x, y, z)
        # Compute distance & set it as semimajor axis length
        semimajor_length = compute_distance(
            0.0, 0.0, 0.0,
            x, y, z
        )

        # Compute rotation angle
        semimajor_angle_y = atan(z/semimajor_length, x/semimajor_length)  # Along y axis
        semimajor_angle_z = asin(y/semimajor_length)  # Along z axis

        return semimajor_length, semimajor_angle_y, semimajor_angle_z
    end



    """
    Compute moment along most distant points
    """
    function compute_moment_most_distant_points(param, index_list, x, y, z)

        # Initialisation
        moment_max = 1.0
        index_smallest_moment = 0

        for itr_index in index_list

            # Compute distance as semimajor axis length & its rotation angle
            semimajor_length, semimajor_angle_y, semimajor_angle_z = compute_semimajor_axis_angle(x[itr_index], y[itr_index], z[itr_index])

            # 5. Apply inverse rotation of point group by semimajor axis angle
            x_z, z_ = compute_rotation_counterclockwise(x, z, -semimajor_angle_y)  # Along y axis
            x_, y_ = compute_rotation_counterclockwise(x_z, y, -semimajor_angle_z)  # Along z axis

            # Compute moment (sum of distance to x axis)
            moment = 0.0
            for itr_point = 1:param.num_points
                moment += compute_distance(
                    0.0, 0.0, 0.0,
                    0.0, y_[itr_point], z_[itr_point]
                )
            end
            moment /= param.num_points

            # Update index with minimum moment
            if moment < moment_max
                moment_max = moment
                index_smallest_moment = itr_index
            end
        end

        return index_smallest_moment
    end


    """
    Find circumscribed spheroid of given point group based on circumscribed sphere
    1. Shift by the centre coordinate
    2. Find the top **% most distant points
    3. Compute moment along most distant points
    4. Define semimajor axis length & angle as origin <-> the point among top top **% most distant points with the smallest moment
    5. Apply inverse rotation of point group by semimajor axis angle
    6. Adjust semimajor axis: x -> x/a
    -----iteration for semiminor axis length (originally sphere radius)
    7. Adjust semiminor axis: y -> y/b
    8. Confirm all points are included in the standard sphere of radius=1
        If so, decrease the semiminor axis length by delta
        If not, exit the loop & define the semiminor axis length
    -----end iteration for semiminor axis
    """
    function search_circumscribed_spheroid(param, points, sphere, spheroid)

        # 1. Shift by the centre coordinate
        x_shift = points.x .- sphere.centre_x
        y_shift = points.y .- sphere.centre_y
        z_shift = points.z .- sphere.centre_z


        # 2. Find the top **% most distant points
        num_distant = search_most_distant_points(
            param,
            x_shift, y_shift, z_shift
        )
        println(num_distant)


        # 3. Compute moment along most distant points
        index_smallest_moment = compute_moment_most_distant_points(
            param, num_distant,
            x_shift, y_shift, z_shift
        )


        # 4. Define semimajor axis length & angle
        distant_x = x_shift[index_smallest_moment]
        distant_y = y_shift[index_smallest_moment]
        distant_z = z_shift[index_smallest_moment]
        semimajor_length, semimajor_angle_y, semimajor_angle_z = compute_semimajor_axis_angle(
            distant_x,
            distant_y,
            distant_z
        )
        println(@sprintf "most distant point with minimum moment ind:%i, x: %.3f, y: %.3f, z: %.3f r: %.3f" index_smallest_moment distant_x distant_y distant_z semimajor_length)


        # 5. Apply inverse rotation of point group by semimajor axis angle
        x_z, z = compute_rotation_counterclockwise(x_shift, z_shift, -semimajor_angle_y)  # Along y axis
        x, y = compute_rotation_counterclockwise(x_z, y_shift, -semimajor_angle_z)  # Along z axis


        # 6. Adjust semimajor axis: x -> x/a
        x_std = x / semimajor_length

        # parameter setting for while-loop
        semiminor_length = semimajor_length
        flag_all_points_in_sphere = true
        delta = 0.01

        while flag_all_points_in_sphere


            # 7. Adjust semiminor axis: y -> y/b & z -> z/b
            y_std = y / semiminor_length
            z_std = z / semiminor_length


            # 8. Confirm all points are included in the sphere
            for itr_point in 1:param.num_points
                dist = compute_distance(
                    0.0, 0.0, 0.0,
                    x_std[itr_point], y_std[itr_point], z_std[itr_point]
                )

                # 8-2. If not, set flag to exit the loop
                if dist > 1
                    flag_all_points_in_sphere = false
                    break
                end
            end

            # 8-1. If so, decrease the semiminor axis length
            if flag_all_points_in_sphere
                semiminor_length -= delta
            end

        end

        println("\n Spheroid Result:")
        println(@sprintf "semimajor %.3f, semiminor %.3f" semimajor_length semiminor_length)
        println(@sprintf "angle along y axis: %.3f, z axis: %.3f" semimajor_angle_y semimajor_angle_z)

        spheroid.semimajor = semimajor_length
        spheroid.semiminor = semiminor_length
        spheroid.angle_y = semimajor_angle_y
        spheroid.angle_z = semimajor_angle_z
    end
end


"""
Module for plot
"""
module Output
    using ..DefineShape:
        sphere_shape,
        sphere_shape,
        spheroid_shape,
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
    Plot points, circumscribed sphere & circumscribed spheroid
    """
    function plot_points_circumscribed(param, points, sphere, spheroid)

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

        # Centre of sphere/spheroid
        p! = scatter!(
            [sphere.centre_x], [sphere.centre_y],
            markercolor = :red
        )

        # Circumscribed sphere
        p! = plot!(
            sphere_shape(sphere.centre_x, sphere.centre_y, sphere.radius),
            seriestype = [:shape,],
            lw = 0.5,
            c = :blue,
            linecolor = :black,
            legend = false,
            fillalpha =0.2,
            aspect_ratio = 1
        )

        # Circumscribed spheroid
        p! = plot!(
            spheroid_shape(
                sphere.centre_x, sphere.centre_y,  # Centre is same as sphere
                spheroid.semimajor, spheroid.semiminor, spheroid.angle),
            seriestype = [:shape,],
            lw = 0.5,
            c = :green,
            linecolor = :black,
            legend = false,
            fillalpha =0.2,
            aspect_ratio = 1
        )

        # Semimajor axis of circumscribed spheroid
        p! = plot!(
            line_shape(
                param,
                sphere.centre_x, sphere.centre_y, spheroid.angle
            ),
            lw = 2.0,
            linecolor = :red,
            legend = false
        )

        savefig(p!, "./tmp/circumscribed.png")
    end


    """
    Output points, circumscribed sphere & circumscribed spheroid in dat file
    """
    function out_points_circumscribed(param, dist, points, sphere, spheroid)

        # Points information
        pointsfile = open("./tmp/points.dat","w")
        for itr_point = 1:param.num_points
            write(pointsfile, "$(points.x[itr_point])\t$(points.y[itr_point])\t$(points.z[itr_point])\n")
        end
        close(pointsfile)

        # Sphere information
        sphere_surface = sphere_shape(sphere)
        lens_sphere = length.(sphere_surface)

        spherefile = open("./tmp/sphere.dat","w")
        for itr_point = 1:lens_sphere[1]
            write(
                spherefile,
                "$(sphere_surface[1][itr_point])\t$(sphere_surface[2][itr_point])\t$(sphere_surface[3][itr_point])\n"
            )
        end
        close(spherefile)

        # Spheroid information
        spheroid_surface = spheroid_shape(sphere, spheroid)
        lens_spheroid = length.(spheroid_surface)

        spheroidfile = open("./tmp/spheroid.dat","w")
        for itr_point = 1:lens_spheroid[1]
            write(
                spheroidfile,
                "$(spheroid_surface[1][itr_point])\t$(spheroid_surface[2][itr_point])\t$(spheroid_surface[3][itr_point])\n"
            )
        end
        close(spheroidfile)

        # Semimajor axis information
        semimajor_line = line_shape(param, dist, sphere, spheroid)
        lens_semimajor = length.(semimajor_line)

        semimajor_linefile = open("./tmp/semimajor.dat","w")
        for itr_point = 1:lens_semimajor[1]
            write(
                semimajor_linefile,
                "$(semimajor_line[1][itr_point])\t$(semimajor_line[2][itr_point])\t$(semimajor_line[3][itr_point])\n"
            )
        end
        close(semimajor_linefile)

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
using .ComputeRotation:
    check_rotation
using .Compute:
    distribute_points,
    search_circumscribed_sphere,
    search_circumscribed_spheroid
using .Output:
    plot_points,
    plot_points_circumscribed,
    out_points_circumscribed



# ----------------------------------------
## Declare parameters
# ----------------------------------------
num_points = 400
x_lim = 1.0

param = ParamVar.Parameters(
    num_points, x_lim
    )


# ----------------------------------------
## Define 3d point group
# ----------------------------------------
semimajor = 0.6 * x_lim
semiminor = 0.2 * x_lim
# No rotation for x axis since it does not affect the spheroid shape
angle_y = π/3.0  # Rotation along y axis
angle_z = π/4.0
shift_x = 0.1  # Shift along x axis
shift_y = 0.2
shift_z = 0.3

dist = ParamVar.Distribution(
    semimajor, semiminor,
    angle_y, angle_z,
    shift_x, shift_y, shift_z
)
ParamVar.ParameterCheck(dist)

println("Given data:")
println(@sprintf "semimajor %.3f, semiminor %.3f" dist.semimajor dist.semiminor)
println(@sprintf "angle along y axis: %.3f, z axis: %.3f" dist.angle_y dist.angle_z)
println(@sprintf "centre x: %.3f, y: %.3f, z: %.3f" dist.shift_x dist.shift_y dist.shift_z)

x, y, z = [Array{Float64}(undef, param.num_points) for _ = 1:3]
points = ParamVar.Points(x, y, z)

distribute_points(param, dist, points)

###CHECK###
# Plot distribution
# plot_points(param, points)
# check_rotation()
###CHECK###


# ----------------------------------------
## Find a circumscribed sphere of point group
## Define centre & radius
## By iterative method
# ----------------------------------------
centre_x, centre_y, centre_z = [0.0 for _ = 1:3]
radius = 0.0
sphere = ParamVar.Sphere(centre_x, centre_y, centre_z, radius)

search_circumscribed_sphere(param, points, sphere)


# ----------------------------------------
## Find a circumscribed spheroid of point group
## Define
## - centre
## - semimajor & semininor axis
## - angle of semimajor axis & x axis
# ----------------------------------------
semimajor, semiminor = [0.0 for _ = 1:2]
angle_y, angle_z = [0.0 for _ = 1:2]
spheroid = ParamVar.Spheroid(
    semimajor, semiminor,
    angle_y, angle_z
)

search_circumscribed_spheroid(param, points, sphere, spheroid)

###CHECK###
# set spheroid for check plot
# spheroid.semimajor = dist.semimajor
# spheroid.semiminor = dist.semiminor
# spheroid.angle_y = dist.angle_y
# spheroid.angle_z = dist.angle_z
###CHECK###

"""
# Plot result in 2D
# Plot points, circumscribed circle & circumscribed ellipse
plot_points_circumscribed(param, points, circle, ellipse)
"""

# Output result
out_points_circumscribed(param, dist, points, sphere, spheroid)
