#=
Julia program to find a circumscribed ellipse for given 2d point data set
=#

"""
Module for global parameters and variables
"""
module ParamVar
    struct Parameters
    end

    mutable struct Variables
    end
end


"""
Module for computation
"""
module Compute
end


"""
Module for plot
"""
module Output
end

# ========================================
# Main function
# ========================================

## Declare modules
using Distributions
using Printf
using Plots
gr(
    markerstrokewidth = 0,
    markersize = 10
)
using .ParamVar
using .Compute
using .Output


# ----------------------------------------
## Declare parameters
# ----------------------------------------


# ----------------------------------------
## Define 2d point group
# ----------------------------------------



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
