# Atmospheric Cluster Dynamics Code

Atmospheric Cluster Dynamics Code, also known as ACDC, is a tool for studying the very first steps of atmospheric new particle formation from gas-phase molecules. ACDC simulates the dynamics of molecular cluster and small nanoparticle populations by generating and solving the cluster birth-death equations for given ambient conditions, yielding the time evolution of cluster concentrations and formation rates.

ACDC can be applied to

* Simulate cluster formation from different chemical compounds using quantum chemical input data
* Study the details of cluster growth processes by, for example, tracking the growth pathways in simulations
* Evaluate experimental methods to interpret measured molecular cluster data by generating synthetic test data
* Include the initial particle formation dynamics in an aerosol dynamics model by direct coupling
* Generate nanoparticle formation rate data to be used in larger-scale models

This repository contains the cluster equation generator, which is the core of ACDC, as well as ready-made templates for easy conducting of standard simulations. Detailed information on all available features can be found in the technical manual, and the quick guide provides the essentials to consider when applying the model on a given set of molecular clusters.

Please note that ACDC is a kinetic model, and does not contain assumptions on the thermodynamic properties of the simulated clusters. Instead, the thermodynamic data is given as input to ACDC. Recent data for clusters of different chemical compositions can be found, for example, in the [Atmospheric Cluster Database (ACDB)](https://github.com/elmjonas/ACDB).

## Citation

If you use the codes provided in this repository in any study, please reference the following papers:

* Simulations of molecular cluster sets by e.g. the standard Matlab tools: Olenius et al.: Free energy barrier in the growth of sulfuric acid–ammonia and sulfuric acid–dimethylamine clusters, J. Chem. Phys. 139, 084312 (2013), https://doi.org/10.1063/1.4819024
* Coupling of ACDC to an aerosol dynamics model by the Fortran routines: Roldin et al.: The role of highly oxygenated organic molecules in the Boreal aerosol-cloud-climate system, Nat. Commun. 10, 4370 (2019), https://doi.org/10.1038/s41467-019-12338-8

## License

This project is licensed under the terms of the GNU General Public License v3.0, as provided with this repository.

Note that the repository contains a few individual functions or routines by other authors. For each such function, the original license or acknowledgement is found in the same folder where the function source code is located and must not be separated from the source code.

## References

The following freely available functions and routines are used in some of the standard ACDC set-ups:

### Matlab functions available in the Matlab File Exchange (https://se.mathworks.com/matlabcentral/fileexchange/)

arrow.m<br/>
Erik Johnson (2020). arrow (https://www.mathworks.com/matlabcentral/fileexchange/278-arrow), MATLAB Central File Exchange. Retrieved October 29, 2020.

axesLabelsAlign3D.m<br/>
Matthew Arthington (2020). Align axes labels in 3D plot (https://www.mathworks.com/matlabcentral/fileexchange/27450-align-axes-labels-in-3d-plot), MATLAB Central File Exchange. Retrieved October 29, 2020.

makeColorMap.m<br/>
Doug Hull (2020). MakeColorMap (https://www.mathworks.com/matlabcentral/fileexchange/17552-makecolormap), MATLAB Central File Exchange. Retrieved October 29, 2020.

### Fortran routines

dvode.f<br/>
The VODE solver for ordinary differential equation systems, Lawrence Livermore National Laboratory (https://computing.llnl.gov/casc/odepack/).<br/>
Brown, P. N., Byrne, G. D., and Hindmarsh, A. C.: VODE, a variable-coefficient ODE solver, SIAM J. Sci. Stat. Comput. 10, 1038-1051 (1989), https://doi.org/10.1137/0910062.

## Prerequisites

In order to use ACDC (on either Windows or Unix), the following software must be installed:

* Perl
* Either Matlab or Fortran

Note that some features, such as the growth pathway tracking, are available only in the Matlab version.

## Authors

Code maintainer: **Tinja Olenius** (tinja.olenius@alumni.helsinki.fi)
