# Proton-beam-SDE
This repository houses a simulator for an SDE describing a the path of a proton in proton beam therapy.
It is under active development and **not fit for standalone use**.

# Installation

In Unix-like environments, simply call `make` in the project root. You will need:

- the g++ compiler (tested on version 13.2.0),
- the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (tested on version 2.7.1),
- the [libconfig](https://hyperrealm.github.io/libconfig/) package (tested on version 1.5).

# Usage

Simulation and output parameters are specified in the config file `sim.cfg`. In particular, you can
modify:
- the number of protons to simulate,
- the radius of the nozzle from which protons are emitted,
- the Gaussian distribution of the initial position of the protons, conditioned to lie in the nozzle,
- the Gaussian distribution of the initial energy of the protons,
- the change points between materials,
- the index of the material in each gap between change points, based on the ordering in Materials/materials.txt,
- the simulator time-step,
- the energy at which a proton is absorbed and its path ends,
- the side length of output voxelation,
- the path to which a sparse 3d grid of dose deposition is stored.

After compilation, run the simulation by calling `./simulate sim.cfg` at the project root.
You can also replace the `sim.cfg` argument with the path to any other appropriate config file.
**Warning**: The default number of protons in `sim.cfg` is set to one million, and voxel sizes and time
steps are also set to a moderate resolution. As a result, the simulation can take 1-2 minutes.

# Post-processing

## Python scripts

Default output paths are specified into the `Output` folder, which contains the Python scripts for output
visualisation. To use it you will need:

- [Python 3](https://www.python.org/) (tested on version 3.12.3), as well as the following Python modules:
- numpy (tested on version 1.26.4),
- matplotlib (tested on version 3.10.0),
- itk (tested on version 5.4.4),
- pymedphys (tested on version 0.40.0),
- scipy (tested on version 1.14.1)

All the necessary functions are written in `SDE_vs_G4.py`, whereas the execution scripts to produce the relevant figures 
is written in `SDE_vs_G4_executionScript.py` (if the Geant4 files are available) or `SDE_plots.py` (for exclusively
plotting SDE output). The scripts include integrated and central-axis depth-dose 1D curves, 2D dose distributions, 
lateral profiles, gamma analysis plots, pass rate calculations and proton range calculations. Simply modify the variables 
defined at the start of the script according to what you need to plot. More detailed information is found in the comments
of each script.

## Geant4 output generation

Geant4 is used as reference to assess the accuracy of the SDE. The source code to generate the outputs that are used for 
comparison are found in the `Geant4-dose-and-secondaries` folder, which has its own README file with a detailed explanation
of the model and how to run the simulations.

# Nuclear cross sections

The `Splines` folder contains splines fitted to nuclear cross section data for hydrogen, oxygen, and nitrogen.
Eventually these will get replaced with a Python script to generate the same cross sections.

# Materials

The interface for specifying materials is via the text files in the `Materials` folder.

- The `atoms.txt` file lists all the atoms needed for all desired materials. Each row specifies an atom name as
well as its atomic number and mass.
- The `materials.txt` lists names of all desired materials. Each material name must be accompanied by a corresponding
`name.txt` file.
- Each `name.txt` file lists three pieces of information. The first two rows contain the material density in `g / cm^3`
and mean excitation energy in `eV`. Rows from the third specify the chemical composition of the material. Each of these
rows consists of an atom ID and the corresponding fraction by mass of that atom in the material. The numerical IDs refer
to the rows of `atoms.txt`, counting from zero.

For example, if the rows of `atoms.txt` are:

hydrogen 1 1\
nitrogen 7 14\
oxygen 8 16

then the stoichiometric rows of water are specified as

0 0.111\
2 0.889

because water is approximately 1/9 hydrogen and 8/9 oxygen by mass.
