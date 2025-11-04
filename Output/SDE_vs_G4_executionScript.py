"""
Calculations and plots for the SDE article
Author: Maria L. Perez-Lara
"""

from SDE_vs_G4 import (retrieve_sde_output, retrieve_g4_output, plot_lateral_profiles,
                       plot_slice, pymedphys_gamma, compare_bragg_peaks, include_new_material,
                       range_comparison, plot_multiple_bragg_peaks)
import numpy as np


"""
Suffix must be selected to compare appropriate files. The options are:
- "": For homogeneous water phantom. 100 and 150 MeV.
- "_bone+h2o": 2 cm bone, water. 100 and 150 MeV.
- "_air+bone+h2o": 5 cm air, 2 cm bone, water. 100 MeV beam only.
"""

N_primaries = "1E6"
suffix = ""
energy1 = 100 # MeV
energy2 = 150 # MeV
comparePhysics = True # Set to true if data from BERT and EMY are available to plot comparison between phys lists
figure_type = ".png"  # Change this to {figure_type} for the article figures
voxel_volume = 0.001 # cm3

# Adjust parameters according to case:
if suffix == "_bone+h2o":
    density_matrix = include_new_material(np.ones((200, 200, 200)), tmin=0, tmax=20, rho=1.45) # add bone
    cuts_100 = [10, 40, 65]
    cuts_150 = [10, 50, 145]
elif suffix == "_air+bone+h2o":
    density_matrix = include_new_material(np.ones((200, 200, 200)), tmin=0, tmax=50, rho=0.001225) # add air
    density_matrix = include_new_material(density_matrix, tmin=50, tmax=70, rho=1.45)  # add bone
    cuts_100 = [25, 60, 115]
    cuts_150 = [0, 0, 0, 0]
else: # Homogeneous case with water
    density_matrix = np.ones((200, 200, 200))
    cuts_100 = [30, 50, 75]
    cuts_150 = [50, 100, 150]

g4_dose1 = (retrieve_g4_output(f"G4_TotalEneDepMap_{N_primaries}_{energy1}MeV{suffix}.csv")
            / (float(N_primaries)*density_matrix * voxel_volume))

sde_dose1 = (retrieve_sde_output(f"SDE_{N_primaries}_{energy1}MeV{suffix}.txt")
             / (float(N_primaries)*density_matrix * voxel_volume))

# Create images for configs with two available beam energies:
if suffix == "" or suffix == "_bone+h2o":
    g4_dose2 = (retrieve_g4_output(f"G4_TotalEneDepMap_{N_primaries}_{energy2}MeV{suffix}.csv")
                / (float(N_primaries) * density_matrix * voxel_volume))
    sde_dose2 = (retrieve_sde_output(f"SDE_{N_primaries}_{energy2}MeV{suffix}.txt")
                 / (float(N_primaries) * density_matrix * voxel_volume))

    # Plots for Energy1:
    f1 = plot_slice(g4_dose1, sde_dose1, 'z', minval=-5, xmin=0, xmax=9, ymin=-5, ymax=5, zmin=-5, zmax=5)
    f2 = plot_lateral_profiles(g4_dose1, sde_dose1, xmin=0, xmax=10, ymin=-5, ymax=5, zmin=-5, zmax=5, lowlim=-1,
                                   uplim=1, cuts=cuts_100)
    print("Gamma values for 100 MeV: ")
    _ = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=10, is_local=False, xmin=0, xmax=9,
                           ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)
    _ = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=10, is_local=True, xmin=0, xmax=9,
                           ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)
    f3 = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=1, is_local=True, xmin=0, xmax=9,
                            ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)
    _ = pymedphys_gamma(g4_dose1, sde_dose1, dta=1, dd=2, th_percent=1, is_local=True, xmin=0, xmax=9,
                           ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)

    f1.savefig(f"../Output/2DSlice_z_{energy1}MeV{suffix}{figure_type}")
    f2.savefig(f"../Output/lateral_profiles_{energy1}MeV{suffix}{figure_type}")
    f3.savefig(f"../Output/gammaTest_{energy1}MeV{suffix}{figure_type}")

    # Plots for Energy2:
    f4 = plot_slice(g4_dose2, sde_dose2, 'z', minval=-5, xmin=0, xmax=18, ymin=-10, ymax=10, zmin=-10, zmax=10)
    f5 = plot_lateral_profiles(g4_dose2, sde_dose2, xmin=0, xmax=18, ymin=-10, ymax=10, zmin=-10, zmax=10, lowlim=-1, uplim=1,
                            cuts=cuts_150)
    print("Gamma values for 150 MeV: ")
    _ = pymedphys_gamma(g4_dose2, sde_dose2, dta=2, dd=3, th_percent=10, is_local=False, xmin=0, xmax=18,
    ymin=-5, ymax=5, zmin=-5, zmax=5)
    _ = pymedphys_gamma(g4_dose2, sde_dose2, dta=2, dd=3, th_percent=10, is_local=True, xmin=0, xmax=18,
    ymin=-5, ymax=5, zmin=-5, zmax=5)
    f6 = pymedphys_gamma(g4_dose2, sde_dose2, dta=2, dd=3, th_percent=1, is_local=True, xmin=0, xmax=18,
                           ymin=-5, ymax=5, zmin=-5, zmax=5)
    _ = pymedphys_gamma(g4_dose2, sde_dose2, dta=1, dd=2, th_percent=1, is_local=True, xmin=0, xmax=18,
    ymin=-5, ymax=5, zmin=-5, zmax=5)

    f4.savefig(f"../Output/2DSlice_z_{energy2}MeV{suffix}{figure_type}")
    f5.savefig(f"../Output/lateral_profiles_{energy2}MeV{suffix}{figure_type}")
    f6.savefig(f"../Output/gammaTest_{energy2}MeV{suffix}{figure_type}")

    # # Comparative plots in 1D:
    f7, f8 = compare_bragg_peaks(sde_dose1, g4_dose1, energy1, sde_dose2, g4_dose2, energy2)
    f7.savefig(f"../Output/1DProj_{energy1}_vs_{energy2}MeV{suffix}{figure_type}")
    f8.savefig(f"../Output/1DSlice_{energy1}_vs_{energy2}MeV{suffix}{figure_type}")

    # Proton range numbers:
    r90_1 = range_comparison(g4_dose1, sde_dose1, 0.9)
    r50_1 = range_comparison(g4_dose1, sde_dose1, 0.5)
    print(f"For {energy1} MeV, proton range agreement is within {r90_1*10:.2f} mm (R90) and {r50_1*10:.2f} mm (R50)")
    r90_2 = range_comparison(g4_dose2, sde_dose2, 0.9)
    r50_2 = range_comparison(g4_dose2, sde_dose2, 0.5)
    print(f"For {energy2} MeV, proton range agreement is within {r90_2*10:.2f} mm (R90) and {r50_2*10:.2f} mm (R50)")

# Comparison with other physics lists, only available for homogeneous phantom configs:
if comparePhysics and suffix == "":
    g4_dose_emy1 = (retrieve_g4_output(f"G4_TotalEneDepMap_{N_primaries}_{energy1}MeV_emy.csv")
                / (float(N_primaries)*density_matrix * voxel_volume))
    g4_dose_bert1 = (retrieve_g4_output(f"G4_TotalEneDepMap_{N_primaries}_{energy1}MeV_bert.csv")
                / (float(N_primaries)*density_matrix * voxel_volume))
    g4_dose_emy2 = (retrieve_g4_output(f"G4_TotalEneDepMap_{N_primaries}_{energy2}MeV_emy.csv")
                / (float(N_primaries)*density_matrix * voxel_volume))
    g4_dose_bert2 = (retrieve_g4_output(f"G4_TotalEneDepMap_{N_primaries}_{energy2}MeV_bert.csv")
                / (float(N_primaries)*density_matrix * voxel_volume))
    f9 = plot_multiple_bragg_peaks(g4_dose1, g4_dose_emy1, g4_dose_bert1, sde_dose1, how='proj', xmax=8,
                         names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"], ref_name='QGSP_BIC_EMZ', maxdif=10)
    f10 = plot_multiple_bragg_peaks(g4_dose1, g4_dose_emy1, g4_dose_bert1, sde_dose1, how='slice', xmax=8,
                          names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"], ref_name='QGSP_BIC_EMZ')

    f11 = plot_multiple_bragg_peaks(g4_dose2, g4_dose_emy2, g4_dose_bert2, sde_dose2, how='proj', xmax=17,
                          names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"], ref_name='QGSP_BIC_EMZ', maxdif=25)
    f12 = plot_multiple_bragg_peaks(g4_dose2, g4_dose_emy2, g4_dose_bert2, sde_dose2, how='slice', xmax=17,
                          names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"], ref_name='QGSP_BIC_EMZ', maxdif=30)

    f9.savefig(f"../Output/1DProj_multiComparison_{energy1}MeV{figure_type}")
    f10.savefig(f"../Output/1DSlice_multiComparison_{energy1}MeV{figure_type}")
    f11.savefig(f"../Output/1DProj_multiComparison_{energy2}MeV{figure_type}")
    f12.savefig(f"../Output/1DSlice_multiComparison_{energy2}MeV{figure_type}")

# # 3 media config, only 100 MeV beam available
if suffix == "_air+bone+h2o":
    f1 = plot_slice(g4_dose1, sde_dose1, 'z', minval=-5, xmin=0, xmax=14, ymin=-7, ymax=7, zmin=-7, zmax=7)
    f1.savefig(f"../Output/2DSlice_z_{energy1}MeV{suffix}{figure_type}")
    f2, f3 = compare_bragg_peaks(sde_dose1, g4_dose1, energy1, max_diff1=20, max_diff2=30)
    f2.savefig(f"../Output/1DProj_{energy1}MeV{suffix}{figure_type}")
    f3.savefig(f"../Output/1DSlice_{energy1}MeV{suffix}{figure_type}")
    f4 = plot_lateral_profiles(g4_dose1, sde_dose1, xmin=0, xmax=18, ymin=-7, ymax=7, zmin=-7, zmax=7, lowlim=-1, uplim=1,
                            cuts=cuts_100)
    f4.savefig(f"../Output/lateral_profiles_{energy1}MeV{suffix}{figure_type}")
    f5 = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=1, is_local=True, xmin=0, xmax=18,
                          ymin=-7, ymax=7, zmin=-7, zmax=7)
    f5.savefig(f"../Output/gammaTest_{energy1}MeV{suffix}{figure_type}")

    # Proton range numbers:
    r90_1 = range_comparison(g4_dose1, sde_dose1, 0.9)
    r50_1 = range_comparison(g4_dose1, sde_dose1, 0.5)
    print(f"For {energy1} MeV, proton range agreement is within {r90_1*10:.2f} mm (R90) and {r50_1*10:.2f} mm (R50)")
