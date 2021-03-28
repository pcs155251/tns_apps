# TNS algorithms
## Introduction
This is a collection of TNS algorithms, including:

* computation of environment of 2D tensor network with translational invariance:
    * Corner transfer matrix (CTM) for iPEPS and 3-PESS.
    * CTM with dimension reduction technique for 3-PESS.
    * Channel environment.
* optimization method: 
    * Simple update for iPEPS and PESS. 
    * Full update for iPEPS.
    * Variational update for iPEPS (under development).
* auxiliary TNS functions
 
Where iPEPS is defined on square lattice and 3-PESS is defined on kagome lattice.

## Usage
* Requirement
    * uni10 v2.1.0: https://gitlab.com/uni10/uni10
* Executable: Below directories include main function and corresponding input parameters which can be compiled and execute.
    * channel: Compute channel environment for 2D transverse Ising model.
    * ipeps: Simple update and full update with CTM environment for transverse Ising model.
    * pess: Simple update and CTM, CTM with dimension reduction scheme for Heisenberg model with DM interations.
The dependency is handled by Makefile in each directory. 
To build the executable, just type `make`.
* Input Files: 
In the directories with executable, there is corresponding \*.rc file where input parameters are set.
For the meaning of each parameters, please refer to class `paraIpeps` and `paraPess`.
* Output Files:
    * The data will stored in data/ automatically
    * The tensors will stored in OutputT/ automatically
* Auxiliary TNS functions:
    * tools/ some useful functions.
    * nsy_ham/ hamiltonians.
    * nsy_op/ operators.
* One can create its own main function based on the class and method which have been implemented.
* For more detail documentation, please see html/index.html.
 
## Referece
* TNS:
    * Annals of Physics 349 (2014) 117-158
* iPEPS:
    * PRL 101, 250602 (2008)
    * Phys. Rev. B 92, 035142 (2015)
* PESS:
    * Phys. Rev. X 4, 011025 (2014)
    * Phys. Rev. B 93, 075154 (2016)
* CTM:
    * J. Phys. Soc. Jpn. 66 (1997) 3040
    * Phys. Rev. B 84, 041108(R) (2011)
* Dimension reduction technique:
    * Phys. Rev. B 96, 045128 (2017)
* Channel environment:
    * Phys. Rev. B 94, 155123 (2016)
 
## Author
Chih-Yuan Lee 
