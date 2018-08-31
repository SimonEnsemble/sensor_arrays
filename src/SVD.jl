#notes
    #MOF_Combinations = n * (n-1) / 2

    #H = [Ha1, Hb1], [Ha2, Hb2]

using CSV, DataFrames

#Read in MOFS

    MOFs = csv.read(blahblahblah)
    MOFs_henrys = MOFs[blahblahblah]

#Perform SVD
    #Should we create the program for comparing just two MOFS?
    #Or should it be able to compare 3 as well?

    #Somehow choose the gas molecules in question to compare MOFS for
    #like literally have an input for the symbol of the gases
    #like was done in the EOS code


    for i = 1:MOFs_henrys
        for j = 1:MOFs_henrys
            #Create Henry's matrix
            #placeholders for the actual variables
            H = [Hai, Hbi; Haj, Hbj]

            #actually do the svd
            #need to store the values constructively
            U, S, V = svd(H)



#analyze outputs

    #error analysis?
    #Delta K/H?

    #graphing
        #get pyplots to work maybe
        #check out what cory used for his beautiful plot
