import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate

from endf_parserpy import EndfParser
Element ="carbon"
CrossSectionFile ="Splines/"+Element+"_cs_el.txt"
AngleFile ="Splines/"+Element+"_ang_el.txt"
A_targ = 12
A_inc = 1
Z_targ = 6
Z_inc = 1
#Rutherford total CrossSection in Barns from -1 to y in [-1,1)
def RutherfordCrossSecCalc(A,a,Z,z,E,y):
    if y<-1 or y>=1:
        print("y out of bounds for Rutherford Cross Section")
        return -1
    A_ratio = A/a
    out = (((1e-2)*2.08*((Z*z)*(1+A_ratio))**2)/((2*A_ratio*E)**2))*(1/(1-y)-1/2)
    return out
def CM_to_Lab_Frame(ang0, E, A, a):
    ang = math.acos(ang0)
    mp = a*938.346  # mass of proton * c^2, MeV
    mn = A*938.346  # mass of colliding nucleus * c^2, MeV
    p = math.sqrt((E) * (E + 2 * mp))
    u = p / (E + mp + mn) #typo should be 2mp
    g = 1 / math.sqrt(1 - u * u)
    e = E + mp
    v_ratio = u * (e - u * p) / (p - u * e)
    if np.fabs((g * (math.cos(ang) + v_ratio))) == 0:
        out = math.pi/2
    else:
        out = math.atan(math.sin(ang) / (g * (math.cos(ang) + v_ratio)))
    if out < 0:
        out += math.pi
    return out


def DiscreteDensityConstructor(endf_dict):
    densities = {}
    incident_energies = endf_dict[6][2]["subsection"][1]["E"].items()
    for key, value in incident_energies:
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
def Total_CDF_constructor(cross_sec_data,spline_PDF,A,a,Z,z):
    Tail_spacing = []
    for i in range(10000):
        Tail_spacing.append((11+5*i-1)/(11+5*i))
    Tail_spacing = np.array(Tail_spacing)
    xx = np.concatenate((np.linspace(-1,0.9,100),Tail_spacing))
    Total_CDF_dict = {}
    for key in spline_PDF.keys():
        if key not in cross_sec_data:
            print("Energy value "+str(key)+" in angle data does not exist in cross section data, skipping")
            continue
        yy = []
        for x in xx:
            yy.append(CM_to_Lab_Frame(x,key*1e-6,A,a))
        total_CDF_key = []
        for x in xx:
            temp_cdf_value = 0
            temp_cdf_value +=cross_sec_data[key]*scipy.interpolate.splint(-1, x, spline_PDF[key])
            temp_cdf_value += RutherfordCrossSecCalc(A,a,Z,z,key*1e-6,x)
            total_CDF_key.append(2*math.pi*temp_cdf_value)
        Total_CDF_dict[key]=[yy,total_CDF_key]
    return Total_CDF_dict
def ErrorCheck(CDF_Data):
    for key in CDF_Data.keys():
        for i in range(1,len(CDF_Data[key][0])):
            if CDF_Data[key][0][i] > CDF_Data[key][0][i-1]:
                print("Error in CM to Lab Conv")
                print(key,i,CDF_Data[key][0][i],CDF_Data[key][0][i-1])
        for i in range(1,len(CDF_Data[key][1])):
            if CDF_Data[key][1][i] <= CDF_Data[key][1][i-1]:
                print("Error in CDF combination") #This may flag due to numerical error, numerical integration step of CDF is negative, take previous value 
                print(key,i,CDF_Data[key][1][i],CDF_Data[key][1][i-1])
                CDF_Data[key][1][i] = CDF_Data[key][1][i-1]
    return CDF_Data
def Threshold_Truncator(CDF_Data,Threshold):
    CDF_Data0 = CDF_Data
    for key in CDF_Data0.keys():
        i = -1
        while CDF_Data0[key][0][i] <= Threshold:
            i -=1
        Distance = (Threshold-CDF_Data0[key][0][(i+1)])/(CDF_Data0[key][0][i]-CDF_Data0[key][0][(i+1)])  
        New_CDF_Value = (1-Distance)*CDF_Data0[key][1][(i+1)]+Distance*CDF_Data0[key][1][i]
        if i <-2:
            Temp_yy = CDF_Data0[key][0][:(i+2)]
            Temp_CDF = CDF_Data0[key][1][:(i+2)]
            CDF_Data0[key]=[Temp_yy,Temp_CDF]
        CDF_Data0[key][0][-1] = Threshold
        CDF_Data0[key][1][-1] = New_CDF_Value
    return CDF_Data0
def Renormalise_CDF(CDF_Data):
    CDF_Data0 = CDF_Data
    IncomingEnergy = []
    PoissonRate = []
    for key in CDF_Data0.keys():
        IncomingEnergy.append(key*1e-6)
        Renormalise_temp = CDF_Data0[key][1][-1]
        PoissonRate.append(Renormalise_temp)
        CDF_Data0[key][0].reverse()
        CDF_Data0[key][1].reverse()
        for i in range(len(CDF_Data0[key][1])):
            CDF_Data0[key][1][i] = 1-(CDF_Data0[key][1][i]/Renormalise_temp)
        CDF_Data0[key][1][0]=0
    return ([IncomingEnergy,PoissonRate],CDF_Data0)



parser = EndfParser()
endf_dictCS = parser.parsefile(CrossSectionFile)
true_incident_energy = []

incident_energy = endf_dictCS[3][2]["xstable"]["E"]
cross_section = endf_dictCS[3][2]["xstable"]["xs"]
cross_sec_data = {}
for i in range(len(incident_energy)):
    cross_sec_data[incident_energy[i]] = cross_section[i]

endf_dict = parser.parsefile(AngleFile)

density = DiscreteDensityConstructor(endf_dict)
for i in density.keys():
    true_incident_energy.append(i)
splines = splineConstructorAngle(3, density)
CDF_Data = Total_CDF_constructor(cross_sec_data,splines,A_targ,A_inc,Z_targ,Z_inc)
CDF_Data = ErrorCheck(CDF_Data)
incident_energyout = []
for key in CDF_Data.keys():
    incident_energyout.append(key)
f = open("Splines/"+Element+"_el_ruth_cross_sec.txt", "w")
tmp = " ".join([str(z*1e-6) for z in true_incident_energy])
f.write(tmp + "\n")
for key in CDF_Data.keys():
    tmp = " ".join([str(z) for z in CDF_Data[key][0]])
    f.write(tmp + "\n")
    tmp = " ".join([str(z) for z in CDF_Data[key][1]])
    f.write(tmp + "\n")
f.close()
