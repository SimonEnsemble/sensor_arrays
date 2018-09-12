#notes
    #MOF_Combinations = n * (n-1) / 2

    #H = [Ha1, Hb1], [Ha2, Hb2]

using CSV, DataFrames # , Gadfly
using LinearAlgebra
using PyPlot

#Read in MOFS

df_henry = CSV.read("../data/henry_constants.csv")

#Perform SVD
function perform_svd(gas1::AbstractString, gas2::AbstractString)
    #Should we create the program for comparing just two MOFS?
    #Or should it be able to compare 3 as well?

    # MOFs[1] indicates the length of the first column, eg. the number of MOFs being screened
    #initialize Henry's Matrix
    sigma = zeros(length(df_henry[1]), length(df_henry[1]), 2)

    for i = 1:length(df_henry[1])
        for j = i+1:length(df_henry[1])
            H = make_h_matrix(i, j, gas1, gas2)
            F = svd(H)
            sigma[i, j, :] = F.S
            end
    end
    return sigma
end

function make_h_matrix(mof1::Int, mof2::Int, gas1::AbstractString, gas2::AbstractString)
    H = zeros(2,2)
    H[1, 1] = df_henry[mof1, Symbol(gas1 * "(KH_mmol/kgPa)")] - df_henry[mof1, Symbol("CH4(KH_mmol/kgPa)")]
    H[1, 2] = df_henry[mof1, Symbol(gas2 * "(KH_mmol/kgPa)")] - df_henry[mof1, Symbol("CH4(KH_mmol/kgPa)")]
    H[2, 1] = df_henry[mof2, Symbol(gas1 * "(KH_mmol/kgPa)")] - df_henry[mof2, Symbol("CH4(KH_mmol/kgPa)")]
    H[2, 2] = df_henry[mof2, Symbol(gas2 * "(KH_mmol/kgPa)")] - df_henry[mof2, Symbol("CH4(KH_mmol/kgPa)")]
    return H
end

function analyze_svd(sigma::Array{Float64, 3})

    #initialize some arrays
    best_indices = argmax(sigma[:, :, 2]) # finds the largest σ₂ value
    worst_indices = argmin(sigma[:, :, 1]) # finds the smallest σ₁ value

    MOF1 = String(df_henry[best_indices[1], 1])
    MOF2 = String(df_henry[best_indices[2], 1])

    println("The most sensitive pair of MOFs is " * MOF1 * " and " * MOF2)

    MOF3 = String(df_henry[worst_indices[1], 1])
    MOF4 = String(df_henry[worst_indices[2], 1])

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

function plot_graphs(mof1::Int, mof2::Int, gas1::AbstractString, gas2::AbstractString)
    H = make_h_matrix(mof1, mof2, gas1, gas2)
    H⁻¹ = inv(H)
    F = svd(H⁻¹)
    θ = range(0, stop=2*π, length=500)[1:end-1]
    p = transpose(hcat(cos.(θ), sin.(θ))) # set pt change
    n = H⁻¹ * p # required input

    function plot_vector(x; head_length=0.05, head_width=0.05, color="k", label="", label_dist=0.05)
        x_plot = x - head_length * x / norm(x)
        arrow(0, 0, x_plot[1], x_plot[2], head_width=head_width,
            head_length=head_length, fc=color, ec=color)
        x_label = x + x / norm(x) * label_dist
        text(x_label[1], x_label[2], label)
    end

    figure(figsize=(8, 5))
    subplot(121)
    scatter(y[1, :], y[2, :], c=θ, cmap="hsv")
    xticks([-1, 0, 1])
    yticks([-1, 0, 1])
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    axis("equal")
    xlabel("\$y_{1, sp}^*\$")
    ylabel("\$y_{2, sp}^*\$")
    for k = 1:2
        plot_vector(SVD[:V][:, k], color="k",
            label="\$\\mathbf{v}_$k\$", label_dist=0.1)
    end
    title("Set pt changes")
    tight_layout()

    subplot(122)
    scatter(u[1, :], u[2, :], c=θ, cmap="hsv")
    if incl_valve_limits
        max_valve = 0.15
        plot([-max_valve, -max_valve], [max_valve, -max_valve], color="k", linestyle="--")
        plot([max_valve, max_valve], [max_valve, -max_valve], color="k", linestyle="--")
        plot([max_valve, -max_valve], [max_valve, max_valve], color="k", linestyle="--")
        plot([max_valve, -max_valve], [-max_valve, -max_valve], color="k", linestyle="--")
    end

    xticks([-1, 0, 1])
    yticks([-1, 0, 1])
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    for k = 1:2
        plot_vector(SVD[:U][:, k] * SVD[:S][k], head_length=0.005,
            head_width=0.05/3, color="k", label="\$\\sigma_$k\\mathbf{u}_$k\$", label_dist=0.025)
    end
    axis("equal")
    xlabel("\$u_{1}^*\$")
    ylabel("\$u_{2}^*\$")
    title("Required inputs")
    tight_layout()
    savefig("u_required.png", format="png", dpi=300)
end
