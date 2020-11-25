# Atmospheric Cluster Dynamics Code

Atmospheric Cluster Dynamics Code, better known as ACDC, is a tool for studying the very first steps of atmospheric new particle formation from gas-phase molecules. ACDC simulates the dynamics of molecular cluster and small nanoparticle populations by generating and solving the cluster birth-death equations for given ambient conditions, yielding the time evolution of cluster concentrations and formation rates.

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
* Coupling of ACDC to an aeorsol dynamics models by the Fortran routines: Roldin et al.: The role of highly oxygenated organic molecules in the Boreal aerosol-cloud-climate system, Nat. Commun. 10, 4370 (2019), https://doi.org/10.1038/s41467-019-12338-8

## Prerequisites

In order to use ACDC (on either Windows or Unix), the following software must be installed:

* Perl
* Either Matlab or Fortran

Note that some features, such as the growth pathway tracking, are available only in the Matlab version.

## Authors

Code maintainer: **Tinja Olenius** (tinja.olenius@alumni.helsinki.fi)

## License

This project is licensed under the terms of the GNU General Public License v3.0, as provided with this repository.

Note that the repository contains a few individual functions by other authors. For each such function, the license is found in the same folder where the function source code is located and must not be separated from the source code.

## Acknowledgments

* The authors of the freely available functions and routines used in some of the standard set-ups are acknowledged for their work
* Sincere thanks to everyone who has participated in the development of the model by direct contributions, feedback or suggestions
