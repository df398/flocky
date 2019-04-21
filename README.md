![flocky](https://github.com/df398/flocky/blob/master/FLOCKY-SMALL.png)


Copyright (c) 2019 David Furman

flocky is an open-source software for efficient development of ReaxFF reactive force fields.

It is based on the recently developed RiPSOGM algorithm to generate high-quality ReaxFF reactive force fields for new systems of interest.

Features:

* Interfaced to the "gold standard" Fortran implementation of ReaxFF developed by Prof. Adri van Duin
* MPI support for an asynchronous parallel performance
* Flexible parameter space exploration: supports use of bounds or percentage change from current values
* Implicit multi-objective fitness function by the use of relative weights
* Training set supports both finite-size and periodic systems
* Automated detection of over-fitting during optimization (soon)
* Bayesian error analysis for uncertainty quantification (soon)


For theoretical background behind RiPSOGM and flocky:

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

(4) Finally, compile flocky:
```bash
make 
```

## Usage
To run flocky, use the corresponding MPI command. Note: N should be less than or equal to the number of swarm agents.
```bash
mpiexec -np N ripsogm
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
