using PorousMaterials
using DataFrames
using JLD2
using CSV
using DelimitedFiles
using Printf
using FileIO

#crystals = readdlm("all_MOFs.txt", String) .* ".cif"
#crystals = readdir("data/crystals/")
#crystals = readdir("../new_xtals/")
crystals = ["BEDYEQ_clean_iqeqcharged.cif", "KAMZUV_manual_iqeqcharged.cif", "LICGOW_iqeqcharged.cif", "MERLAZ_clean_iqeqcharged.cif", "NIBBOR_clean.cif", "REWNEO_clean_iqeqcharged.cif",
            "RUBTAK01_clean_h_iqeqcharged.cif", "VOGTIV_clean_h_iqeqcharged.cif", "WIZMAV_iqeqcharged.cif", "ORIWET_iqeqcharged.cif", "OCUNAC_clean5b_iqeqcharged.cif"]
old_crystals = ["BEDYEQ_clean.cif", "KAMZUV_manual.cif", "LICGOW.cif", "MERLAZ_clean.cif", "NIBBOR_clean.cif", "REWNEO_clean.cif",
            "RUBTAK01_clean_h.cif", "VOGTIV_clean_h.cif", "WIZMAV.cif", "ORIWET.cif", "OCUNAC_clean5b.cif"]

data_to_collect = ["henry coefficient [mmol/(g-bar)]", "err henry coefficient [mmol/(g-bar)]", "Qst (kJ/mol)"]
gases = ["CO2", "SO2"]

df = DataFrame(crystal=crystals[:])
for gas in gases
    for data in data_to_collect
        println(gas, "_", data)
        df[Symbol(gas * "_" * data)] = [-1.0 for i = 1:length(crystals)]
        # allow missing
        df[Symbol(gas * "_" * data)] = convert(Array{Union{Float64, Missing}, 1}, df[Symbol(gas * "_" * data)])
    end
end

insertions_per_volume = 250
ljforcefield = "UFF.csv"
temperature = 298.0

for gas in gases
    for crystal in crystals
        result_file = "data/henry_sims/" * henry_result_savename(Framework(crystal), Molecule(gas), temperature,
                                           LJForceField(ljforcefield), insertions_per_volume)
        if isfile(result_file)
            result = load(result_file)
            result = result[collect(keys(result))[1]]
#            @load result_file result

            idx_crystal = df[:crystal] .== crystal
            for data in data_to_collect
                if result[data] < eps()
                    df[idx_crystal, Symbol(gas * "_" * data)] = 0.0
                else
                    df[idx_crystal, Symbol(gas * "_" * data)] = result[data]
                end
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
@printf("Got here\n")
#=
for gas in ["C2H6", "CH4"]#"p_xylene", "o_xylene", "m_xylene", "C2H6", "CH4"]
    for (i, crystal) in enumerate(crystals)
        idx_crystal = df[:crystal] .== crystal
        # If something is missing, lets see if we can grab the info from the non-charged crystal
        #println(df[idx_crystal, Symbol(gas *"_henry coefficient [mmol/(g-bar)]")])
        if ismissing(df[idx_crystal, Symbol(gas *"_henry coefficient [mmol/(g-bar)]")][1])
            result_file = "data/henry_sims/" * henry_result_savename(Framework(old_crystals[i]), Molecule(gas), temperature, LJForceField(ljforcefield), insertions_per_volume) 
            if isfile(result_file)
                result = load(result_file)
                result = result[collect(keys(result))[1]]
                for data in data_to_collect
                    if result[data] < eps()
                        df[idx_crystal, Symbol(gas * "_" * data)] = 0.0
                    else
                        df[idx_crystal, Symbol(gas * "_" * data)] = result[data]
                    end
                end
                println("Found " * gas * " in old " * old_crystals[i])
            else
                println("Didn't find " *gas * "in old " * old_crystals[i])
            end
        end
    end
end
=#
CSV.write("sensors_crystal_KHs_CO2SO2.csv", df)
