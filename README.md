# MLPrecon
Machine-Learned Preconditioners for Linear Solvers in Geophysical Flows


Step 1) Compiling:

Compile module "sheusp_functions.f90":
gfortran -O2 -c sheusp_functions.f90

Compile shallow-water model versions that use different preconditioners:
No Preconditioner:                sheusp_impl_NOprecon_1M10_res05.f90
Implicit Richardson:              sheusp_impl_ImpRichardprecon_1M10_res05.f90
ML Precon; perturbed run:         sheusp_impl_MLprecon_1M10_perturbed_res05.f90
ML Precon; initial period:        sheusp_impl_MLprecon_1M10_initial_15days_res05.f90
ML Precon reduced; perturbed run: sheusp_impl_MLprecon_reduced_1M10_perturbed_res05.f90

Example using no preconditioner:

gfortran -O2 -c sheusp_impl_NOprecon_1M10_res05.f90

gfortran -O2 -o sheusp_impl_NOprecon_1M10_res05.exe sheusp_functions.o sheusp_impl_NOprecon_1M10_res05.o


Step 2) Run Simulation

Create required directories and run the respective shallow-water model version
The weights required by the machine-learned preconditioners are provided in directory "weights".


Step 3) Create Convergence plot

Convergence plots were done using:
plot_residual_Convergence_MaxMinMedian.py
Directory names and simulation periods are specified in the script header.
scritp was run with python3.6
