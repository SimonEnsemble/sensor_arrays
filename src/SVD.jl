#notes
    #MOF_Combinations = n * (n-1) / 2

    #H = [Ha1, Hb1], [Ha2, Hb2]

using CSV, DataFrames, Gadfly

#Read in MOFS

    MOFs = CSV.read("C:\\Users\\Caleb\\Documents\\GitHub\\sensor_arrays\\data\\henry_constants.csv")

#Perform SVD
function perform_svd(gas1::Symbol, gas2::Symbol)
    #Should we create the program for comparing just two MOFS?
    #Or should it be able to compare 3 as well?

    #initialize Henry's Matrix
    H = zeros(2,2)
    sigma = Array{Array{Float64, 1}}(length(MOFs[1]), length(MOFs[1]))

    for i = 1:length(MOFs[1])

        gas1 = String(gas1)
        gas2 = String(gas2)

        H[1,1] = MOFs[i, Symbol(gas1*"(KH_mmol/kgPa)")]
        H[1,2] = MOFs[i, Symbol(gas2*"(KH_mmol/kgPa)")]


        for j = 1:length(MOFs[1])
            #Should I erase the entries where the MOFs are the same?
            #or just leave them and do the analysis anyways?
            if i == j
                sigma[i,j] = [0,0]
                continue
            end

            H[2,1] = MOFs[j, Symbol(gas1*"(KH_mmol/kgPa)")]
            H[2,2] = MOFs[j, Symbol(gas2*"(KH_mmol/kgPa)")]

            U, S, V = svd(H)
            sigma[i, j] = S
        end
    end

    return sigma
end

function analyze_svd(sigma)

    #initialize some arrays
    mag = zeros(length(MOFs[1]), length(MOFs[1]))
    highest = 0
    highest_index = [0,0]

    for i = 1:length(MOFs[1])
        for j = 1:length(MOFs[1])
            #TODO Change to max, min
            mag[i, j] = norm(sigma[i, j])
            if mag[i, j] > highest
                highest = mag[i, j]
                highest_index = [i,j]
            end
        end
    end

    MOF1 = String(MOFs[highest_index[1], 1])
    MOF2 = String(MOFs[highest_index[2], 1])

    println("The most sensitive pair of MOFs is " * MOF1 * " and " * MOF2)
    #error analysis?
    #Delta K/H?

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
end
