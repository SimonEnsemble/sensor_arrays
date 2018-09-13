using PorousMaterials
using DataFrames
using CSV
using DelimitedFiles
using Printf

# defines number of insertions of a molecule into a single crystal
insertions_per_volume =
# defines Lennard Jones forcefield
ljforcefield = ".csv"
# defines temperature
temperature =
# defines crystals as the folder which stores all crystals
crystals = readdlm("all_crystals.txt", String)

"""
Write a job submission script to submit for KH simulation of `gas` in `crystals` at `temperature`.
"""
function write_henry_submit_script(molecule::String, crystal::String, temperature::Float64,
                                   ljforcefield::String, insertions_per_volume::Int)
    jobscriptdir = "jobz"
    if ! isdir(jobscriptdir)
        mkdir(jobscriptdir)
    end

    # when this script is run in the cluster, it will print this off as a prompt
    @printf("Writing submission script for Henry Coefficient simulation of %s in %s
            at %f K with %s with %d Widom insertions.\n",
            molecule, crystal, temperature, ljforcefield, insertions_per_volume)

    # build KH_submit.sh script to be run in the cluster
    KH_submit = open("KH_submit.sh", "w")
    @printf(KH_submit,
    """
    #!/bin/bash

    # use current working directory for input and output
    # default is to use the users home directory
    #\$ -cwd

    # name this job
    #\$ -N %s

    #\$ -pe thread 4 # use 4 threads/cores

    # send stdout and stderror to this file
    #\$ -o jobz/%s.o
    #\$ -e jobz/%s.e

    # select queue - if needed; mime5 is SimonEnsemble priority queue but is restrictive.
    ##\$ -q mime5

    # print date and time
    date
    julia -p 4 run_henry.jl %s %s %f %s %d
    """, crystal * "_" * molecule * "_" * string(insertions_per_volume),
        crystal * "_" * molecule, crystal * "_" * molecule, molecule, crystal * ".cif",
        temperature, ljforcefield, insertions_per_volume)
    close(KH_submit)
end

for gas in ["CO2", "CH4"]
    # reads in a crystal from crystals folder
	for crystal in crystals
       	 # submission file filled in by molecule and crystal of interest
       	 write_henry_submit_script(gas, crystal, temperature, ljforcefield, insertions_per_volume)
       	 # runs script to the cluster
       	 run(`qsub KH_submit.sh`)
       	 sleep(1)
	end
end
