# RiPSOGM

RiPSOGM is a software package (written in C++11) for efficient parametrization of ReaxFF reactive force fields.

It is based on a improved version of the Rotation Invariant Particle Swarm Optimization algorithm with less risks of being trapped in local minima solutions due to the use of (Gaussian) mutation operators.

Features:

* Interfaced to the "gold standard" stand-alone implementation of ReaxFF
* MPI support for a parallel performance
* Flexible parameter space exploration with on-the-fly validity tests
* Implicit multi-objective fitness function 
* Training set supports both finite-size and periodic systems


For theoretical background behind RiPSOGM, please see:


David Furman, Benny Carmeli, Yehuda Zeiri and Ronnie Kosloff,
[Enhanced Particle Swarm Optimization Algorithm: Efficient Training of ReaxFF Reactive Force Fields](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01272),
J. Chem. Theory Comput., 2018, 14 (6)

A user's manual can be downloaded from: [link](http://insertlink)


## Installation
Requirements:

* OpenMPI / MPICH2
* Boost C++ 

To compile RiPSOGM, adjust the provided makefile accordingly to your system configuration (point to path of boost library), and execute the following:
```bash
make 
```


## Usage
To use RiPSOGM,  execute it with the corresponding MPI task shceduler (mpiexec,mpirun,etc.) where N is the number of processors to run on. For efficiency, N should be less than or equal to the number of swarm agents.
```bash
mpiexec -np N ./ripsogm_public
```


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
