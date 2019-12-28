<p align="center">
  <img width="460" height="300" src="https://github.com/df398/flocky/blob/master/flocky-logo.jpg">
</p>

Copyright (c) 2019 David Furman

flocky is an open-source platform for development of high-quality ReaxFF reactive force fields for reactive molecular dynamics simulations. flocky is based on the RiPSOGM swarm intelligence algorithm for meta-heuristic global optimization.

Features:

* Interfaced to the standalone Fortran implementation of tapered ReaxFF
* MPI support for an asynchronous parallel performance
* Flexible parameter space exploration: supports use of bounds or percentage change from current values
* Implicit multi-objective fitness function by the use of relative weights
* Training set supports both finite-size and periodic systems
* On the fly detection of over-fitting during optimization
* On the fly Bayesian errors analysis for uncertainty quantification


For theoretical background behind tapered ReaxFF, RiPSOGM and flocky:

David Furman and David J. Wales,
[Transforming the Accuracy and Numerical Stability of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.9b02810),
J. Phys. Chem. Lett. 2019, 10 (22)

David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput. 2018, 14 (6)

Click [here](https://df398.github.io/flocky/) for documentation and usage examples.


## Installation
Requirements:
* Any recent version of GCC
* OpenMPI (tested with 1.10.2 and 4.0.1) or MPICH2 for parallel flocky
* [Boost C++](https://www.boost.org/) (tested with 1.70.0)
* [OptimLib](https://www.kthohr.com/optimlib.html#installation-method-1-shared-library) and its prerequisites: [Armadillo](http://arma.sourceforge.net/download.html)

Install procedure:

(1) Install OpenMPI or MPICH for parallel execution (for a serial version, skip to step 2):

On Ubuntu to install mpich type:
```bash
sudo apt-get install mpich
```

For OpenMPI type:
```bash
sudo apt-get install openmpi-bin
```

(2) Install Boost C++ required libraries and header files:
```bash
cd path/to/boost_1_61_0
./bootstrap.sh --with-libraries=system,filesystem --prefix=path/to/installation/prefix
./b2 install link=static
```

(3) Set your BOOST_LIB_PATH and BOOST_INC_PATH environment variables:
```bash
export BOOST_LIB_PATH="path/to/installation/prefix/lib" 
export BOOST_INC_PATH="path/to/installation/prefix/include"
```

(4) Install OptimLib as a shared library:
```bash
export ARMA_INCLUDE_PATH=/path/to/armadillo
./configure -i "/path/to/install-dir" -p -m "-lblas"
make
make install
```
For full OpenMP support, consider using OpenBLAS (i.e. replace "-lblas" with "-lopenblas")

(5) Finally, compile flocky:
```bash
make -f Makefile.mpi
```
for a parallel version, or:
```bash
make -f Makefile.serial
```
for a serial version.

## Usage
To run flocky in parallel mode, use the corresponding MPI command. Note: N should be less than or equal to the number of swarm agents.
```bash
mpiexec -np N flocky_mpi
```
or, launch a serial version using:
```bash
./flocky_serial
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to contribute.

## License
flocky is distributed under the terms of the [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/) license.

## Citations
If you find flocky useful, please consider citing:

David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput., 2018, 14 (6)
