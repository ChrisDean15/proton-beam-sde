import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate

from endf_parserpy import EndfParser

Threshold = 0.1


def Nuclear_trig_constants(n, eta, sign):  # recursively evaluate antiderivative
    if n == 0:
        return [0, 1]
    else:
        sin_const = (
            sign
            * (
                (n - 1) * Nuclear_trig_constants(n - 1, eta, sign)[0]
                + eta * Nuclear_trig_constants(n - 1, eta, sign)[1]
            )
            / ((n - 1) ** 2 + eta**2)
        )
        cos_const = (
            sign
            * (
                (n - 1) * Nuclear_trig_constants(n - 1, eta, sign)[1]
                - eta * Nuclear_trig_constants(n - 1, eta, sign)[0]
            )
            / ((n - 1) ** 2 + eta**2)
        )
        return [sin_const, cos_const]


def Trig_Integral_Eval(
    deg, cosval, sinval, y, eta, sign
):  # evaluation of antiderivative at a given value y
    coeff = Nuclear_trig_constants(deg, eta, sign)
    return ((1 + (y * sign)) ** (deg - 1)) * (coeff[0] * sinval + coeff[1] * cosval)


def CDF_Poly_Eval(deg, eta, c, y, sign):  # evaluation of integration by parts form
    cosval = np.cos(eta * np.log((1 + (y * sign)) / 2) + c)
    sinval = np.sin(eta * np.log((1 + (y * sign)) / 2) + c)
    temp = 0
    for i in range(deg + 1):
        temp += (
            Trig_Integral_Eval(i + 1, cosval, sinval, y, eta, sign)
            * (y ** (deg - i))
            * math.factorial(deg)
            * ((-1) ** (i))
            / math.factorial(deg - i)
        )
    temp += (
        ((-1) ** (deg + 1))
        * math.factorial(deg)
        * Trig_Integral_Eval(deg + 1, cosval, sinval, 0, eta, sign)
    )
    return temp


def Legendre_Coeff_Cal(deg):  # outputs coefficients of legendre polynomial as a list
    if deg == 0:
        return [1]
    elif deg == 1:
        return [0, 1]
    else:
        temp1 = Legendre_Coeff_Cal(deg - 1)
        temp1.insert(0, 0)
        for i in range(len(temp1)):
            temp1[i] = (2 * (deg - 1) + 1) * temp1[i] / deg
        temp2 = Legendre_Coeff_Cal(deg - 2)
        temp2.append(0)
        temp2.append(0)
        for i in range(len(temp2)):
            temp2[i] = (1 - deg) * temp2[i] / deg
        temp3 = []
        for i in range(len(temp2)):
            temp3.append(temp1[i] + temp2[i])
        return temp3


def CDF_Legendre_Eval(
    deg, eta, c, y, sign
):  # evaluation of integration by parts for Legendre polynomial
    Legendre_temp = Legendre_Coeff_Cal(deg)
    out = 0
    for i in range(deg):
        out += Legendre_temp[i] * CDF_Poly_Eval(deg, eta, c, y, sign)
    return out


def Legendre_int(deg, y):  # evaluation of Legendre polynomial integral
    Legendre_temp = Legendre_Coeff_Cal(deg)
    Total = 0
    for i in range(len(Legendre_temp)):
        Total += (y ** (i + 1)) / (i + 1)
    return Total


def Complex_Polar(re, im):
    r = np.sqrt(re**2 + im**2)
    theta = np.arctan2(im, re)
    return [r, theta]


def eta_eval(Z, z, E, m):
    sqrtterm = 0.0496 * m / (E * 2)
    out = Z * z * np.sqrt(sqrtterm)
    return out


def Total_Rate_calc(
    Z, z, E, A, a, y, two_s, a_coef_re, a_coef_im, b_coef
):  # total rate between 0 and y
    A_ratio = A / a
    eta_val = eta_eval(Z, z, E, a)
    Rutherford_Coef = ((1e-2) * 4.16 * ((Z * z) * (1 + A_ratio)) ** 2) / (
        (2 * A_ratio * E) ** 2
    )
    Rate = 0
    Rate += Rutherford_Coef * (
        (y / (1 - y**2))
        + ((((-1) ** two_s) / (two_s + 1)))
        * math.sin(2 * eta_val * math.log((1 + y) / (1 - y)))
        / (2 * eta_val)
    )
    Nuclear_term_1 = 0
    for i in range(len(a_coef_im)):
        polar_temp = Complex_Polar(a_coef_re[i], a_coef_im[i])
        int_temp = 0
        int_temp += ((-1) ** i) * CDF_Legendre_Eval(i, eta_val, polar_temp[1], y, 1)
        int_temp += CDF_Legendre_Eval(i, eta_val, polar_temp[1], y, -1)
        Nuclear_term_1 += polar_temp[0] * (2 * i + 1) * int_temp / 2
    Nuclear_term_1 *= 2 * eta_val
    Rate -= Nuclear_term_1
    Nuclear_term_2 = 0
    for i in range(len(a_coef_im)):
        Nuclear_term_2 += (4 * i + 1) * b_coef[i] * Legendre_int(2 * i, y) / 2
    Rate += Nuclear_term_2
    Rate *= 2 * np.pi
    return Rate


def Threshold_Truncator(
    ang_data, Rate_data, Threshold
):  # truncates CDF at both ends to avoid explosion in rutherford crosseciton
    i = -1
    while ang_data[i] <= Threshold:
        i -= 1
    if i == -1:
        print("error end point does not reach threshold")
    Distance = (Threshold - ang_data[(i + 1)]) / (ang_data[i] - ang_data[(i + 1)])
    New_CDF_Value = (1 - Distance) * Rate_data[(i + 1)] + Distance * Rate_data[i]
    if i < -2:
        ang_data = ang_data[: (i + 2)]
        Rate_data = Rate_data[: (i + 2)]
    ang_data[-1] = Threshold
    Rate_data[-1] = New_CDF_Value
    return ang_data, Rate_data


def ErrorCheck(Data, sign):  # checks for monotonicity in CDF and CM to Lab frame
    if sign == 1:
        for i in range(1, len(Data)):
            if Data[i] < Data[i - 1]:
                print("Error in CDF monotonicity")
    if sign == -1:
        for i in range(1, len(Data)):
            if Data[i] >= Data[i - 1]:
                print("Error in CM to Lab monotonicity")
                print(Data[i], Data[i - 1])
    if sign != -1 and sign != 1:
        print("sign value error")
    return 0


def CM_to_Lab_Frame(ang0, E, A, a):
    ang = math.acos(ang0)
    mp = a * 938.346  # mass of proton * c^2, MeV
    mn = A * 938.346  # mass of colliding nucleus * c^2, MeV
    p = math.sqrt((E) * (E + 2 * mp))
    u = p / (E + mp + mn)  # typo should be 2mp
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


parser = EndfParser()
endf_dict = parser.parsefile("Hydrogen_elastic_data.txt")
density = {}
angle_discretize = np.linspace(0, 0.9, 100)
angle_discretize_end = []
for i in range(
    600
):  # Setup interpolation points to match growth of rutherford crossection near boundary
    angle_discretize_end.append(np.sqrt((5.5 + 5 * i - 1) / (5.5 + 5 * i)))
angle_discretize_pos = np.concatenate((angle_discretize, angle_discretize_end))
Energy_discrete = []
Rate_discrete = []
for i in endf_dict[6][2]["subsection"][1]["E"].keys():
    print(i)
    Total_rate = []
    Real_coef = []
    Imaginary_coef = []
    b_coef = []
    Energy = endf_dict[6][2]["subsection"][1]["E"][i] * 1e-6
    if Energy >= 1:
        for r in range(8, 22):
            if r % 2 == 0:
                Real_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
            else:
                Imaginary_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
        for r in range(1, 8):
            b_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
        for ang in angle_discretize_pos:
            Total_rate.append(
                Total_Rate_calc(
                    1, 1, Energy, 1, 1, ang, 1, Real_coef, Imaginary_coef, b_coef
                )
            )
        #  CM_frame_angle, Total_rate = Threshold_Truncator(CM_frame_angle, Total_rate,Threshold)
        ErrorCheck(Total_rate, 1)
        Energy_discrete.append(Energy)
        max_rate = Total_rate[-1]
        Rate_discrete.append(max_rate)
        Total_rate_True = []
        angle_discrete_True = []
        Bool_temp = True
        for i in range(len(Total_rate)):
            if True:
                angle_discrete_True.insert(
                    0, CM_to_Lab_Frame(-angle_discretize_pos[i], Energy, 1, 1)
                )
                Total_rate_True.insert(0, -Total_rate[i])
            angle_discrete_True.append(
                CM_to_Lab_Frame(angle_discretize_pos[i], Energy, 1, 1)
            )
            Total_rate_True.append(Total_rate[i])
        max_rate = -Total_rate_True[0]
        for i in range(len(Total_rate_True)):
            Total_rate_True[i] += max_rate
        density[Energy] = [angle_discrete_True, Total_rate_True]
f = open("hydrogen_el_ruth_cross_sec.txt", "w")
tmp = " ".join([str(z) for z in Energy_discrete])
f.write(tmp + "\n")
bool = True
for key in density.keys():
    tmp = " ".join([str(z) for z in density[key][0]])
    f.write(tmp + "\n")
    tmp = " ".join([str(z) for z in density[key][1]])
    f.write(tmp + "\n")
    bool = False
f.close()
