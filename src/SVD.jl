#notes
    #MOF_Combinations = n * (n-1) / 2

    #H = [Ha1, Hb1], [Ha2, Hb2]

using CSV, DataFrames

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
            #if i == j
            #    sigma[i,j] = [0,0]
            #    continue
            #end

            H[2,1] = MOFs[j, Symbol(gas1*"(KH_mmol/kgPa)")]
            H[2,2] = MOFs[j, Symbol(gas2*"(KH_mmol/kgPa)")]

            #actually do the svd
            #need to store the values constructively
            U, S, V = svd(H)
            sigma[i, j] = S
        end
    end

    return sigma
end

function analyze_svd(sigma)

    mag = zeros(length(MOFs[1]), length(MOFs[1]))
    highest = 0
    highest_index = [0,0]
    for i = 1:length(MOFs[1])
        for j = 1:length(MOFs[1])
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

    #graphing
        #get pyplots to work maybe
        #check out what cory used for his beautiful plot
end
