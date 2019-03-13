using CSV
using DataFrames
using Printf

f = CSV.read("MOF_KHs.csv")
dropmissing!(f)

ch4_greater = f[Symbol("CH4_henry coefficient [mmol/(g-bar)]")] .> f[Symbol("CO2_henry coefficient [mmol/(g-bar)]")]

mofs = f[ch4_greater, :]

CSV.write("new_MOF_KHs.csv", mofs)
