import numpy as np
import os

#Script to output non-elastic cross section data ready to be parsed into SDE model
#Cross sections are evaluated as described in https://iopscience.iop.org/article/10.1088/1361-6560/ae5586 DOI 10.1088/1361-6560/ae5586

from endf_parserpy import EndfParser
parser = EndfParser()
File_name = "Materials/atoms.txt"
element_data = []

#Reads in atom data, should be in the format Element, Z, A, extra empty lines will throw error
with open(File_name) as f:
    for line in f:
        word, num1, num2 = line.split()
        element_data.append((word, int(num1), float(num2)))
for data in element_data:
    #Element info + Raw ENDF data file locations
    Element = data[0]
    CrossSectionFile = "Splines/ENDF_data/" + Element + "_cs_ne.txt"
    EnergyAngleFile = "Splines/ENDF_data/" + Element + "_enang_ne.txt"
    if os.path.exists(CrossSectionFile) and os.path.exists(EnergyAngleFile):
        endf_dict = parser.parsefile(EnergyAngleFile)
        densities_out = {}
        #key item pairs for incidient energy
        for key, item in endf_dict[6][5]["subsection"][2]["E"].items():
            goinginen = item
            denstemp = 0
            denstemp2 = []
            energytemp = []
            r = []
            #key item pairs for energy,density,rvalue at given incident energy
            for key2 in endf_dict[6][5]["subsection"][2]["b"][key]:
                #density value
                denstemp += endf_dict[6][5]["subsection"][2]["b"][key][key2][0]
                #output energy value
                energytemp.append(endf_dict[6][5]["subsection"][2]["Ep"][key][key2])
                denstemp2.append(denstemp)
                #r value 
                r.append(endf_dict[6][5]["subsection"][2]["b"][key][key2][1])
            for j in range(len(denstemp2)):
                denstemp2[j] *= 1 / denstemp
            denstemp2[-1] = 1.0
            densities_out[item] = [energytemp, denstemp2, r]
        max_length = 0
        for key in densities_out.keys():
            #This should not flag
            if len(densities_out[key][0]) != len(densities_out[key][1]) or len(densities_out[key][0]) != len(densities_out[key][2]):
                print("Error in CDF length")
                print(key, len(densities_out[key][0]), len(densities_out[key][1]))
            max_length = max(max_length, len(densities_out[key][0]))
        for key in densities_out.keys():
            #Include if all output rows are to be equal length. 
            while len(densities_out[key][0]) < max_length and False:
                densities_out[key][0].append(densities_out[key][0][-1])
                densities_out[key][1].append(densities_out[key][1][-1])
                densities_out[key][2].append(densities_out[key][2][-1])
        endf_dict = parser.parsefile(CrossSectionFile)
        avaenergies = []
        avaenergies2 = []
        for key in densities_out.keys():
            avaenergies.append(key)
            avaenergies2.append(key / 1e6)
        f = open("Splines/" + Element + "_ne_rate.txt", "w")
        #cross section energy
        f.write(" ".join([str(z / 1e6) for z in endf_dict[3][5]["xstable"]["E"]]) + "\n")
        #cross section r value
        tmp = " ".join([str(max(z, 0)) for z in endf_dict[3][5]["xstable"]["xs"]])
        f.write(tmp + "\n")
        f.close()
        f = open("Splines/" + Element + "_ne_energyangle_cdf.txt", "w")
        f.write(" ".join([str(z) for z in avaenergies2]) + "\n")
        for keys in densities_out.keys():
            tmp = " ".join([str(z / 1e6) for z in densities_out[keys][0]])
            f.write(tmp + "\n")
            tmp = " ".join([str(max(z, 0)) for z in densities_out[keys][1]])
            f.write(tmp + "\n")
            tmp = " ".join([str(z) for z in densities_out[keys][2]])
            f.write(tmp + "\n")
        f.close()
