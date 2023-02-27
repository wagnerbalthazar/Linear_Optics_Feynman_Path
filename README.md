# Linear_Optics_Feynman_Path

In this repository we find a c++ code to simulate Linear Optical Boson Sampling using Feynman's formalism of Quantum Dynamics. 
A pdf file is available for a naive understanding of the method and with a complete description of how to use the code. Read the pdf file 'Linear Optics Feynman Path'. The paper with a complete description of the method will be available soon.
The codes available use different functions. All of them are in a same file to make things readable and user-friendly. The user can only make changes in the 'main' function to run the code.

Contact: wagner.balthazar@inl.int; wagner.balthazar@ifrj.edu.br
-------------------------------------------------------------------------------


**Prerequisites:
** This code was written using the following versions:
* c++ 20;
* The code runs in other early versions of c++, but the Coroutine library must be installed since Coroutine was added as standard only in c++20.

*** Brief description: ** *
The modules supplied allow us to execute two different tasks (as described in our pre-print []):
** Main routine **
* Probability Feynman (“probabilty_feynman.cpp”) – This package contains a source file with all functions to calculate the probability of a Fock-state boson sampling experiment for chosen input/output states, and a given interferometer made out of a finite number of locally connected beam-splitter layers.
** Additional routines **
* Probability Ryser (“probability_Ryser.cpp”) - This code runs to calculate the probability of boson sampling for one amplitude using Ryser's formula to calculate the permanent.
* Probability Glynn (“probability_Glynn.cpp”) - This code runs to calculate the probability of boson sampling for one amplitude using Glynn's formula to calculate the permanent. 
* Example 1 (”Example_1_correctness.cpp”) - The code is used to verify the correctness of the code by comparing Feynman, Glynn, and Ryser.
* Example 2 (”Example_2_increase_modes_photons.cpp”) - The code is used to verify for different depths the runtime code when the number of modes is increasing always
with 1 photon per mode. 
* Example 3 (”Example_3_increase_depth.cpp”) - The code is used to study the runtime and memory for different layers with a constant number of modes.
* Example 4 (”Example_4_Higher_input.cpp”) - is used to study the runtime and memory increasing the photon density, i.e. the number of photons per mode.
* Example 5 (”Example_5_read_csv_bs_parameters.cpp”) - The code is used to call the beam splitters parameters from a matrix in a CSV file. This aims to make the code easier to use.

**User instructions**:
To use the code you need to read the pdf file ‘Linear-Optical Feynman Paths (LOFP): C++ code for strong simulation of Fock-state boson sampling’. In section II we present the codes available. In section III we have presented the physics background. In section IV, there is a complete explanation of Running the code. In section V, we have presented the benchmark. Finally, the conclusions about the package.
We run the code in Visual Studio, an integrated development environment (IDE). It is very easy to run the code in Visual Studio if you never worked with c++. You only need to create a project. Change the language for c++ 20. Open the source file available and run the code. Another option is to use a compiler and text editor of your preference. 
  
**Citation**
message: "If you use this software, please cite it as below."
authors: Wagner F Balthazar; Ernesto F Galvão.
  
orcid: https://orcid.org/0000-0002-0135-1668; https://orcid.org/0000-0003-2582-9918
  
title: Linear-Optical Feynman Paths (LOFP): code for strong simulation of Fock-state boson sampling

version: 1.0.0

doi: 10.5281/zenodo.7681675

date-released: 2023-02-27

url: "https://github.com/wagnerbalthazar/Linear_Optics_Feynman_Path/"
     
url: "https://zenodo.org/record/7681675#.Y_zeRXbP2Uk"

We acknowledge the financial support of H2020-FETOPEN Grant PHOQUSING (GA no.: 899544).

**Copyright:**
Copyright (C) 2022  W.F.B
Creative Commons Attribution 4.0 International
