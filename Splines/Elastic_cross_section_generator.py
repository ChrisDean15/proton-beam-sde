import math
import numpy as np
import scipy
import scipy.interpolate
import os
from endf_parserpy import EndfParser
#Script to output Elastic cross section data ready to be parsed into SDE model
#Cross sections are evaluated as described in https://iopscience.iop.org/article/10.1088/1361-6560/ae5586 DOI 10.1088/1361-6560/ae5586
parser = EndfParser()
File_name = "Materials/atoms.txt"
element_data = []
#Cuts off cross section data at this forward scattering angle in RAD. Can be kept small as the variable parameter rutherford_cutoff in the main code
# can be used to modify this cutoff variably. 
cut_off = 0.0001
#Spherical brownian motion scattering density diverges from Moliere scattering density at 2.5x standard deviation
#Exclude cross section data above Gaussian_SD_cut_off's of this standard deviation. Should be at least 2.5x
#but can be fitted parameter
Gaussian_SD_cut_off = 2.5
#Reads in atom data, should be in the format Element, Z, A, extra empty lines will throw error
with open(File_name) as f:
    for line in f:
        word, num1, num2 = line.split()
        element_data.append((word, int(num1), float(num2)))
for data in element_data:
    #Element info + Raw ENDF data file locations
    Element = data[0]
    CrossSectionFile = "Splines/ENDF_data/" + Element + "_cs_el.txt"
    AngleFile = "Splines/ENDF_data/" + Element + "_ang_el.txt"
    #Checks if ENDF files exist before processing
    if os.path.exists(CrossSectionFile) and os.path.exists(AngleFile):
        A_targ = data[2]
        A_inc = 1
        Z_targ = data[1]
        Z_inc = 1



        # SD of multiple scattering contribution, individual scattering should start at 2.5SD. Taken with track length step 0.5mm (first constant in chi_c_sq).
        def multiple_scattering_sd(e, z,a):
            mpcsq = 938.346  # mass of proton * speed of light squared, MeV
            pv = (2 * mpcsq + e) * e / (mpcsq + e)
            betasq = (2 * mpcsq + e) * e / ((mpcsq + e) ** 2)
            p = pv / math.sqrt(betasq) #momentum in MeV / c.
            chi_c_sq = 0.05*z*(z+1)*0.157/(a*pv*pv)
            chi_a_sq_vec = 2.007e-5 * (z**(2 / 3)) * (1 + 3.34 * (z / (137 * math.sqrt(betasq)))**2) / (p * p)
            omega = chi_c_sq / chi_a_sq_vec
            F = 0.98
            v = omega / (2 * (1 - F))
            ret = math.sqrt(
                chi_c_sq * (((1 + v) * math.log(1 + v) / v) - 1) / (1 + F * F)
            )
            return ret
        # Rutherford total CrossSection in Barns from -1 to y in [-1,1)
        def RutherfordCrossSecCalc(A, a, Z, z, E, y):
            if y < -1 or y >= 1:
                print("y out of bounds for Rutherford Cross Section")
                return -1
            A_ratio = A / a
            #constant (1e-2)2.08 is evaluated as alpha^2 hbar^2c^2 1e16 as in Section 6.2.6 of the ENDF data manual
            out = (
                ((1e-2) * 2.08 * ((Z * z) * (1 + A_ratio)) ** 2) / ((2 * A_ratio * E) ** 2)
            ) * (1 / (1 - y) - 1 / 2)
            return out


        def CM_to_Lab_Frame(ang0, E, A, a):
            ang = math.acos(ang0)
            mp = a * 938.346  # mass of proton * c^2, MeV
            mn = A * 938.346  # mass of colliding nucleus * c^2, MeV
            p = math.sqrt((E) * (E + 2 * mp))
            u = p / (E + mp + mn) 
            g = 1 / math.sqrt(1 - u * u)
            e = E + mp
            v_ratio = u * (e - u * p) / (p - u * e)
            if np.fabs((g * (math.cos(ang) + v_ratio))) == 0:
                out = math.pi / 2
            else:
                out = math.atan(math.sin(ang) / (g * (math.cos(ang) + v_ratio)))
            if out < 0:
                out += math.pi
            return out


        def DiscreteDensityConstructor(endf_dict):
            densities = {}
            #incident energies
            incident_energies = endf_dict[6][2]["subsection"][1]["E"].items()
            for key, value in incident_energies:
                #outgoing angle densities for given incident energy
                AngleProb = endf_dict[6][2]["subsection"][1]["A"][key].items()
                temp1 = []  # store angles
                temp2 = []  # store P
                for key2, value2 in AngleProb:
                    if key2 % 2 == 1:
                        temp1.append(value2)
                    else:
                        temp2.append(value2)
                temp1[-1] = 1  # final angle is slightly away from one, converting to 1 for ease
                densities[value] = [temp1, temp2]
            return densities


        def splineConstructorAngle(k, densities):
            spline_PDF = {}
            for key in densities.keys():
                splt = scipy.interpolate.splrep(densities[key][0], densities[key][1], k=k)
                spline_PDF[key] = splt
            return spline_PDF


        def Total_CDF_constructor(cross_sec_data, spline_PDF, A, a, Z, z, cutoff, Gaussian_SD_cutoff):
            Tail_spacing_1 = []
            Tail_spacing_2 = []
            for i in range(39): # Setup interpolation points to match growth of rutherford crossection near boundary
                Tail_spacing_1.append((11 + 10 * i - 1) / (11 + 10 * i))
            for i in range(9000): #Sparser interpolation points if needed to minimise data size
                Tail_spacing_2.append((401 + 400 * i - 1) / (401 + 400 * i))
            Tail_spacing_1 = np.array(Tail_spacing_1)
            Tail_spacing_2 = np.array(Tail_spacing_2)
            Tail_spacing = np.concatenate((Tail_spacing_1, Tail_spacing_2))
            xx = np.concatenate((np.linspace(-1, 0.90, 100), Tail_spacing))
            Total_CDF_dict = {}
            for key in spline_PDF.keys():
                if key not in cross_sec_data:
                    print(
                        "Energy value for "
                        + str(Element) + " at " + str(1e-6*key)
                        + " MeV in angle data does not exist in cross section data, skipping"
                    )
                    continue
                yy = []
                for x in xx:
                    if CM_to_Lab_Frame(x, key * 1e-6, A, a) < max(Gaussian_SD_cutoff * multiple_scattering_sd(key * 1e-6, Z, A),cutoff):
                        break
                    yy.append(CM_to_Lab_Frame(x, key * 1e-6, A, a))
                total_CDF_key = []
                for x in xx:
                    if CM_to_Lab_Frame(x, key * 1e-6, A, a) < max(Gaussian_SD_cutoff * multiple_scattering_sd(key *1e-6, Z, A),cutoff):
                        break
                    temp_cdf_value = 0
                    temp_cdf_value += cross_sec_data[key] * scipy.interpolate.splint(
                        -1, x, spline_PDF[key]
                    )
                    #integral density of rutherford cross section (6.9) Section 6.2.6 of ENDF manual
                    temp_cdf_value += RutherfordCrossSecCalc(A, a, Z, z, key * 1e-6, x)
                    if temp_cdf_value < 0:
                        temp_cdf_value = 0
                    total_CDF_key.append(2 * math.pi * temp_cdf_value)
                if len(yy) == 0 or total_CDF_key[-1] == 0:
                    Total_CDF_dict[key] = [[0],[1e-10]] #scattering is zero, set small scattering angle with low rate
                else:
                    Total_CDF_dict[key] = [yy, total_CDF_key]
            max_length = 0
            for key in Total_CDF_dict.keys():
                max_length = max(max_length, len(Total_CDF_dict[key][0]))
                # This should not flag
                if len(Total_CDF_dict[key][0]) != len(Total_CDF_dict[key][1]):
                    print("Error in CDF length")
            for key in Total_CDF_dict.keys():
                #Include if all output rows are to be equal length.
                while len(Total_CDF_dict[key][0]) < max_length and False:
                    Total_CDF_dict[key][0].append(Total_CDF_dict[key][0][-1])
                    Total_CDF_dict[key][1].append(Total_CDF_dict[key][1][-1])
            return Total_CDF_dict


        def ErrorCheck(CDF_Data):
            for key in CDF_Data.keys():
                for i in range(1, len(CDF_Data[key][0])):
                    # This should not flag
                    if CDF_Data[key][0][i] > CDF_Data[key][0][i - 1]:
                        print("Error in CM to Lab Conv")
                        print(data,key, i, CDF_Data[key][0][i], CDF_Data[key][0][i - 1])
                for i in range(1, len(CDF_Data[key][1])):
                    # This may flag due to numerical error and is not a concern, numerical integration step of CDF is negative, take previous value
                    if CDF_Data[key][1][i] < CDF_Data[key][1][i - 1]:
                    #    print(
                    #        "Error in CDF combination"
                    #    )  
                    #    print(key, i, CDF_Data[key][1][i], CDF_Data[key][1][i - 1])
                        CDF_Data[key][1][i] = CDF_Data[key][1][i - 1]
            return CDF_Data
        endf_dictCS = parser.parsefile(CrossSectionFile)
        true_incident_energy = []
        #incident energy value
        incident_energy = endf_dictCS[3][2]["xstable"]["E"]
        #rate value
        cross_section = endf_dictCS[3][2]["xstable"]["xs"]
        cross_sec_data = {}
        for i in range(len(incident_energy)):
            cross_sec_data[incident_energy[i]] = cross_section[i]

        endf_dict = parser.parsefile(AngleFile)

        density = DiscreteDensityConstructor(endf_dict)
        splines = splineConstructorAngle(3, density)
        CDF_Data = Total_CDF_constructor(cross_sec_data, splines, A_targ, A_inc, Z_targ, Z_inc, cut_off, Gaussian_SD_cut_off)
        CDF_Data = ErrorCheck(CDF_Data)
        incident_energyout = []
        for key in CDF_Data.keys():
            incident_energyout.append(key)
        f = open("Splines/" + Element+"_el_ruth_cross_sec.txt", "w")
        tmp = " ".join([str(z * 1e-6) for z in incident_energyout])
        f.write(tmp + "\n")
        for key in CDF_Data.keys():
            tmp = " ".join([str(z) for z in CDF_Data[key][0]])
            f.write(tmp + "\n")
            tmp = " ".join([str(z) for z in CDF_Data[key][1]])
            f.write(tmp + "\n")
        f.close()
