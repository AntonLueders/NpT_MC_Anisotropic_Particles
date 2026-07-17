# Introduction

This repository contains the NpT MC simulation software used to obtain the data in the manuscript

[1] Lüders, A., Franosch, T., and Schilling, R. *How the number of vertices determines the thermodynamics of fluids of hard polygons*.

The different folders contain individual programs for the different hard convex particles. All information and methods employed are described in detail in the supporting information of [1]. While the programs can be used to repeat the corresponding simulations, we strongly encourage a full independent reproduction using the information in the supplemental information.

# Table of contents:
 - [ Requirements and how to build the executables](#Requirements)
 - [How to use the MC programs](#Use)
 - [How to interprete the generated data](#Data)
 - [Disclaimer](#Disclaimer)
 - [Cite this software](#Cite)

 <a id="Requirements"></a>
# Requirements and how to build the executables

The program is written in *C* (C99 standard) and it can successfully be built with gcc. The simulations presented in [1] were performed using binaries build with *gcc version 4.8.5* and *gcc version 7.5.0*. To compile BDHM, the library GSL must be linked. 

The source files of the different programs can be found in the correspondings folders in this repository. An example makefiles are given. As usual, it can be used to compile the program by typing 

``
make
``

in the console while being in the folder that contains the source files.

 <a id="Use"></a>
# How to use the MC programs

The particular MC programs must be executed in the console using the ".out" file. Note that the parameters of the simulation must be stated in a file called **inputfile** which is in the same folder as the binary. Note that this input file is "static" meaning that the parameters must always be stated at the same position. If inputfile is not present in the folder where the MC program is executed, the program exits with an error. To generate a new example input file, use "Setup" as additional command when excecuting the binary in the console. The input file has always the form of the following example:

```
N	2000
p	5.000000
T	1.000000
sweeps	10000000
rate	10000
a       0.886
seed	1
delta	0.40000
delta_r	1.00000
delta_V	1.100000
r_v	5.7500
m	4
lambda	1.0
```

Here, particular input lines can be vary between the different shapes (the example shows the input for polygons with a finite number of vertices).

The particular lines of such an input file and the meaning of the numbers are summarized in the following list:

- **N:** Number of particles
- **p:** Pressure
- **T:** Temperature (the program is only tested for the case where this value is set to unity)
- **sweeps:** Number of MC sweeps (i.e., steps per particle)
- **rate:** The particle configuration is saved every "rate" sweeps.
- **a:** Length scale of the particles. For polygons, this is the distance between the center and the vertices for the regular reference body before deformation. For ellipses, this is the larger semi-axis.
- **seed:** Random seed.
- **delta:** Displacement for translational random moves.
- **delta_r:** Displacements for rotational random moves.
- **delta_V:** Increment corresponding to the MC step of the system area.
- **r_v:** Verlet-radius of the used Verlet lists.
- **m:** Number of vertices.
- **lambda:** Scaling factor used to deform the regular polygons.

Note that the word and number separation in the input file should be done with the **tab key** (the program is only tested for this case).  

 <a id="Data"></a>
# How to interprete the generated data

All simulation programs in this repository output a single data file:

- *Position_For_Ovito.xyz:* This file stores the particle configurations. Every "rate" MC sweeps, the center coordinates are stored in this file. The last column of the data file contains the current total area of the simulation box, which is needed to calculate the average particle density.

To visualize the configurations of the system, it is suggested to use external software that can interpret the xyz position files. For this, *VMD* [Humphrey et al. J. Mol. Graph. 1996] or *OVITO* [https://www.ovito.org/] can be used.

 <a id="Disclaimer"></a>
# Disclaimer

Note that the presented code is **not** an ongoing software project with active support. Instead, it is part of the supplemental information of [1] and reflects the state of the software used to obtain the corresponding numerical results. For any questions regarding [1] and the corresponding data, please contact the corresponding authors.

 <a id="Cite"></a>
# Cite this software

If you use this simulation software or parts of it, please cite

[1] Lüders, A., Franosch, T., and Schilling, R. *How the number of vertices determines the thermodynamics of fluids of hard polygons*.
