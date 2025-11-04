"""
Calculations and plots for the SDE article
This script does not include Geant4 data, it is used for SDE data visualisation purposes.
Author: Maria L. Perez-Lara
"""

from SDE_vs_G4 import (
    retrieve_sde_output,
    plot_lateral_profiles,
    plot_slice,
    compare_bragg_peaks,
    include_new_material,
    plot_bragg_peaks_SDE,
)
import numpy as np


"""
Suffix must be selected to compare appropriate files. The options are:
- "": For homogeneous water phantom. 100 and 150 MeV.
- "_bone+h2o": 2 cm bone, water. 100 and 150 MeV.
- "_air+bone+h2o": 5 cm air, 2 cm bone, water. 100 MeV beam only.
"""

N_primaries = "1E6"
suffix = ""
energy1 = 100  # MeV
energy2 = 150  # MeV
voxel_volume = 0.001  # cm3
figure_type = ".png"  # Change this to .eps for the article figures

# Adjust parameters according to case:
if suffix == "_bone+h2o":
    density_matrix = include_new_material(
        np.ones((200, 200, 200)), tmin=0, tmax=20, rho=1.45
    )  # add bone
    cuts_100 = [10, 40, 65]
    cuts_150 = [10, 50, 145]
elif suffix == "_air+bone+h2o":
    density_matrix = include_new_material(
        np.ones((200, 200, 200)), tmin=0, tmax=50, rho=0.001225
    )  # add air
    density_matrix = include_new_material(
        density_matrix, tmin=50, tmax=70, rho=1.45
    )  # add bone
    cuts_100 = [25, 60, 115]
    cuts_150 = [0, 0, 0, 0]
else:  # Homogeneous case with water
    density_matrix = np.ones((200, 200, 200))
    cuts_100 = [30, 50, 75]
    cuts_150 = [50, 100, 150]

sde_dose1 = retrieve_sde_output(f"SDE_{N_primaries}_{energy1}MeV{suffix}.txt") / (
    float(N_primaries) * density_matrix * voxel_volume
)

# Create images for configs with two available beam energies:
if suffix == "" or suffix == "_bone+h2o":
    sde_dose2 = retrieve_sde_output(f"SDE_{N_primaries}_{energy2}MeV{suffix}.txt") / (
        float(N_primaries) * density_matrix * voxel_volume
    )

    # Plots for Energy1:
    f1 = plot_slice(
        None,
        sde_dose1,
        "z",
        minval=-5,
        xmin=0,
        xmax=9,
        ymin=-5,
        ymax=5,
        zmin=-5,
        zmax=5,
    )
    f2 = plot_lateral_profiles(
        None,
        sde_dose1,
        xmin=0,
        xmax=10,
        ymin=-5,
        ymax=5,
        zmin=-5,
        zmax=5,
        lowlim=-1,
        uplim=1,
        cuts=cuts_100,
    )

    f1.savefig(f"../Output/SDE_2DSlice_z_{energy1}MeV{suffix}{figure_type}")
    f2.savefig(f"../Output/SDE_lateral_profiles_{energy1}MeV{suffix}{figure_type}")

    # Plots for Energy2:
    f4 = plot_slice(
        None,
        sde_dose2,
        "z",
        minval=-5,
        xmin=0,
        xmax=18,
        ymin=-10,
        ymax=10,
        zmin=-10,
        zmax=10,
    )
    f5 = plot_lateral_profiles(
        None,
        sde_dose2,
        xmin=0,
        xmax=18,
        ymin=-10,
        ymax=10,
        zmin=-10,
        zmax=10,
        lowlim=-1,
        uplim=1,
        cuts=cuts_150,
    )

    f4.savefig(f"../Output/SDE_2DSlice_z_{energy2}MeV{suffix}{figure_type}")
    f5.savefig(f"../Output/SDE_lateral_profiles_{energy2}MeV{suffix}{figure_type}")

    # Comparative plots in 1D:
    f7, f8 = plot_bragg_peaks_SDE(sde_dose1, energy1, sde_dose2, energy2)
    f7.savefig(f"../Output/SDE_1DProj_{energy1}_vs_{energy2}MeV{suffix}{figure_type}")
    f8.savefig(f"../Output/SDE_1DSlice_{energy1}_vs_{energy2}MeV{suffix}{figure_type}")

# 3 media config, only 100 MeV beam available
if suffix == "_air+bone+h2o":
    f1 = plot_slice(
        None,
        sde_dose1,
        "z",
        minval=-5,
        xmin=0,
        xmax=14,
        ymin=-7,
        ymax=7,
        zmin=-7,
        zmax=7,
    )
    f1.savefig(f"../Output/SDE_2DSlice_z_{energy1}MeV{suffix}{figure_type}")
    f2, f3 = plot_bragg_peaks_SDE(sde_dose1, energy1)
    f2.savefig(f"../Output/SDE_1DProj_{energy1}MeV{suffix}{figure_type}")
    f3.savefig(f"../Output/SDE_1DSlice_{energy1}MeV{suffix}{figure_type}")
    f4 = plot_lateral_profiles(
        None,
        sde_dose1,
        xmin=0,
        xmax=18,
        ymin=-7,
        ymax=7,
        zmin=-7,
        zmax=7,
        lowlim=-1,
        uplim=1,
        cuts=cuts_100,
    )
    f4.savefig(f"../Output/SDE_lateral_profiles_{energy1}MeV{suffix}{figure_type}")
