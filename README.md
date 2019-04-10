# RiPSOGM

RiPSOGM is a software package (written in C++11) for simple and efficient parametrization of ReaxFF reactive force fields.

It is based on the Rotation Invariant Particle Swarm Optimization algorithm and is augmented with Gaussian mutation operators to enhance its exploration abilities.
This provides an efficient way of generating new, high-quality ReaxFF reactive force fields for new systems of interest.

Features:

* Interfaced to the "gold standard" Fortran implementation of ReaxFF developed by Prof. Adri van Duin
* MPI support for a asynchronous parallel performance; swarm members are distributed across processors
* Flexible parameter space exploration: supports use of bounds or percentage change from current values for parameters
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

(2) Install Boost C++:
* Follow the instructions in [install_boost](https://www.boost.org/doc/libs/1_61_0/more/getting_started/unix-variants.html). 
Note that RiPSOGM requires the compilation of boost::filesystem and boost:system static libraries in addition to the header files in boost.

(3) Adjust the PATHS for Boost in the provided Makefile

(3) To compile RiPSOGM:
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
RiPSOGM is distributed under ther terms of the GNU GPLv3 license.
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
