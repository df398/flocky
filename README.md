# RiPSOGM

RiPSOGM is a software package for a simple and efficient parametrization of ReaxFF reactive force fields.

It is based on the Rotation Invariant Particle Swarm Optimization algorithm and is augmented with Gaussian mutation operators to enhance its exploration abilities.
This provides an efficient way of generating new, high-quality ReaxFF reactive force fields for new systems of interest.

Features:

* Interfaced to the "gold standard" Fortran implementation of ReaxFF developed by Prof. Adri van Duin
* MPI support for an asynchronous parallel performance; swarm members are distributed across processors
* Flexible parameter space exploration: supports use of bounds or percentage change from current values
* Implicit multi-objective fitness function by the use of relative weights
* Training set supports both finite-size and periodic systems


For theoretical background behind RiPSOGM, please see:

David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput., 2018, 14 (6)

A user's manual can be downloaded from: [link](http://insertlink)


## Installation
Requirements:

* OpenMPI or MPICH2
* Boost C++

Install procedure:

(1) Install OpenMPI or MPICH for parallel execution:

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

(4) Finally, compile RiPSOGM:
```bash
make 
```

## Usage
To run RiPSOGM, use the corresponding MPI command. Note: N should be less than or equal to the number of swarm agents.
```bash
mpiexec -np N ripsogm_public
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to contribute.

## License
RiPSOGM is distributed under the terms of the GNU GPLv3 license.
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## Citations
If you find RiPSOGM useful, please consider citing:
David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput., 2018, 14 (6)
