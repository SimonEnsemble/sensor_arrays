using PorousMaterials
using DataFrames
using JLD2
using HDF5
using CSV
using DelimitedFiles
using Printf

crystals = readdlm("all_MOFs.txt", String) .* ".cif"

data_to_collect = ["henry coefficient [mmol/(g-bar)]", "err henry coefficient [mmol/(g-bar)]", "Qst (kJ/mol)", "elapsed time (min)"]
gases = ["CO2", "C2H6", "CH4"]

df = DataFrame(crystal=crystals[:])
for gas in gases
    for data in data_to_collect
        df[Symbol(gas * "_" * data)] = [-1.0 for i = 1:length(crystals)]
        # allow missing
        df[Symbol(gas * "_" * data)] = convert(Array{Union{Float64, Missing}, 1}, df[Symbol(gas * "_" * data)])
    end
end

insertions_per_volume =
ljforcefield = ".csv"
temperature =

for gas in gases
    for crystal in crystals
#=
        if crystal == "CON.cif"
            @warn "ASK ARNI ABT CON"
            idx_crystal = df[:crystal] .== crystal
            for data in data_to_collect
                df[idx_crystal, Symbol(gas * "_" * data)] = missing
            end
            continue
=#        end
        result_file = "data/henry_sims/" * henry_result_savename(Framework(crystal), Molecule(gas), temperature,
                                           LJForceField(ljforcefield), insertions_per_volume)
        if isfile(result_file)
            @load result_file result

            idx_crystal = df[:crystal] .== crystal
            for data in data_to_collect
                df[idx_crystal, Symbol(gas * "_" * data)] = result[data]
            end
            println("Found " *  gas * " in " * crystal)
        else
            idx_crystal = df[:crystal] .== crystal
            for data in data_to_collect
                df[idx_crystal, Symbol(gas * "_" * data)] = missing
            end
            @warn("did not find " *  gas * " in " * crystal)
        end
    end
end
CSV.write("sensors_crystal_KHs.csv", df)
