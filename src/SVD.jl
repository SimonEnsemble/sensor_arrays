#notes
    #MOF_Combinations = n * (n-1) / 2

    #H = [Ha1, Hb1], [Ha2, Hb2]

using CSV, DataFrames # , Gadfly
using LinearAlgebra

#Read in MOFS

    MOFs = CSV.read("../data/henry_constants.csv")

#Perform SVD
function perform_svd(gas1::AbstractString, gas2::AbstractString)
    #Should we create the program for comparing just two MOFS?
    #Or should it be able to compare 3 as well?

    # MOFs[1] indicates the length of the first column, eg. the number of MOFs being screened
    #initialize Henry's Matrix
    H = zeros(2,2)
    sigma = zeros(length(MOFs[1]), length(MOFs[1]), 2)

    for i = 1:length(MOFs[1])

        H[1,1] = MOFs[i, Symbol(gas1 * "(KH_mmol/kgPa)")]
        H[1,2] = MOFs[i, Symbol(gas2 * "(KH_mmol/kgPa)")]


        for j = i+1:length(MOFs[1])

            H[2,1] = MOFs[j, Symbol(gas1 * "(KH_mmol/kgPa)")]
            H[2,2] = MOFs[j, Symbol(gas2 * "(KH_mmol/kgPa)")]

            U, sigma[i, j, :], V = svd(H)
        end
    end

    return sigma
end

function analyze_svd(sigma::Array{Float64, 3})

    #initialize some arrays
    best_indices = argmax(sigma[:, :, 2]) # finds the largest σ₂ value
    worst_indices = argmin(sigma[:, :, 1]) # finds the smallest σ₁ value

    MOF1 = String(MOFs[best_indices[1], 1])
    MOF2 = String(MOFs[best_indices[2], 1])

    println("The most sensitive pair of MOFs is " * MOF1 * " and " * MOF2)

    MOF3 = String(MOFs[worst_indices[1], 1])
    MOF4 = String(MOFs[worst_indices[2], 1])

    println("The least sensitive pair of MOFs is " * MOF3 * " and " * MOF4)

    #error analysis?
    #Delta K/H?

#=
    #create circle array
    n = 1000
    #TODO change to range
    x_lin = linspace(-1,1,n)
    y = zeros(2*n)
    x = zeros(2*n)
    for i = 1:n
        y[2 * i] = sqrt(1 - (x_lin[i] ^ 2))
        y[2 * i - 1] = - sqrt(1 - (x_lin[i] ^ 2))
        x[2 * i] = x_lin[i]
        x[2 * i - 1] = x_lin[i]
    end

    xx = x
    yy = y

    pre_plot = plot(x = xx, y = yy, Geom.point, Guide.xlabel("x"), Guide.ylabel("y"))

    #formats the xy coords appropriately
    xy = hcat(x, y)
    xy = transpose(xy)

    #stretches the xy coords by the sigma from the svd
    highest_sigma = sigma[highest_index[1], highest_index[2]]
    sigma = [highest_sigma[1] 0; 0 highest_sigma[2]]
    println(sigma)
    xy_stretched = sigma*xy

    post_plot = plot(x = xy_stretched[1,:], y = xy_stretched[2,:], Geom.point, Guide.xlabel("x"), Guide.ylabel("y"))

    return post_plot
    =#
end
