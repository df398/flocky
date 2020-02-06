<p align="center">
  <img src="https://github.com/df398/flocky/blob/master/flocky-logo.jpg">
</p>

flocky is an open-source platform for development of high-quality ReaxFF reactive force fields. It is based on the RiPSOGM swarm intelligence algorithm for global optimization.

Features:

* Interfaced to the Tapered ReaxFF fortran code
* 2-level MPI parallelisation over swarm members and over training set
* Flexible parameter space exploration: supports use of bounds or percentage change from current values
* Implicit multi-objective fitness function uses relative weights
* Powerful local minimisation algorithms
* Regularization methods (L1 or L2) and on-the-fly detection of overfitting during optimization
* On-the-fly Bayesian errors analysis for uncertainty quantification
* Training set supports both finite-size and periodic systems


For theoretical background behind Tapered ReaxFF, RiPSOGM and flocky, please consult:

David Furman and David J. Wales,
[Transforming the Accuracy and Numerical Stability of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.9b02810),
J. Phys. Chem. Lett. 2019, 10 (22)

David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput. 2018, 14 (6)

Click [here](https://df398.github.io/flocky/) for additional documentation, usage examples and a user forum.


## Installation
Requirements:

* A recent C++ compiler that supports the C++11 standard
* A fortran compiler (for Tapered ReaxFF)
* [OpenMPI](https://www.open-mpi.org/) or [MPICH2](https://www.mpich.org/) for parallel performance (highly recommended)
* [Boost C++](https://www.boost.org/)
* [NLopt](https://github.com/stevengj/nlopt)

Install procedure:

(0) Install Tapered ReaxFF:

```bash
cd Tapered-ReaxFF-for-flocky
make -f Makefile.my
cd ..
```

(1) Install OpenMPI or MPICH2 for parallel performance (to build a serial version, skip to step 2):

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
export BOOST_LIB_PATH=path/to/installation/prefix/lib
export BOOST_INC_PATH=path/to/installation/prefix/include
```

(4) Install NLopt as a shared library:
```bash
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=path/to/preferred/installdir
make
make install
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:path/to/preferred/installdir/lib
export NLOPT_PATH=path/to/preferred/installdir
```

(5) Finally, compile flocky with the appropriate makefile:
```bash
make -f Makefile.mpi
```
or:
```bash
make -f Makefile.serial
```

## Usage
To run flocky in parallel mode, use the corresponding MPI command.
```bash
mpiexec -np N flocky_mpi
```
or, launch a serial version using:
```bash
flocky_serial
```
> note: **tapreaxff** has to be present in the working directory. To run **flocky** on HPC cluster machines, refer to the [documentation](https://df398.github.io/flocky/).

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to contribute.

## License
flocky is distributed under the terms of the [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/) license.
Copyright © 2019 David Furman, University of Cambridge.

## Citations
If you used flocky for your research, please cite the relevant publications:

David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput., 2018, 14 (6)

David Furman and David J. Wales,
[Transforming the Accuracy and Numerical Stability of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.9b02810),
J. Phys. Chem. Lett. 2019, 10 (22)


