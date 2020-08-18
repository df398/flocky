/*  ---------------------------------------------------------------------- *
    flocky v1.0 Copyright (C) 2018-2020 David Furman, PhD.
    df398@cam.ac.uk, University of Cambridge, UK.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
* ------------------------------------------------------------------------ */

#pragma once
#ifndef PAR_H
#define PAR_H
#include <stdio.h>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_numeric.hpp>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/process.hpp>
#include <boost/process/async_system.hpp>
#include <nlopt.hpp>
#include <map>
# define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
# undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <regex>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <include/pcg_random.hpp>

using namespace std;

extern int ierr, core, numcores;
extern int mycore_swarmcore;
extern int size_swarmcores;
extern int mycore_reaxffcore;
extern int size_reaxffcores;

extern int ptrainset;
extern vector <int> swarmcores;
extern int swarmcores_size;
extern vector <int> reaxffcores;
extern int reaxffcores_size;
extern int totmembers;
extern int dim;
extern int NumP;
extern double c1;
extern double c2;
extern double inertiamax;
extern double inertiamin;
extern bool chang;
extern bool verbose;
extern bool lg_yn;
extern bool contff;
extern bool fixcharges;
extern int regular;
extern double dropvalue;
extern double hlambda;
extern bool ofit;
extern bool uq;
extern bool ensembleave;
extern int preppath;
extern bool gbfitfound;
extern bool printonce;
extern bool firstovfit;
extern int freq;
extern double levyscale;
extern int parid_gbfit;
extern int faili;
extern int ninf;
extern int iter;
extern int maxiters;
extern int cycle;
extern int maxcycles;
extern bool perc_yn;
extern double perc;
extern int localmin;
extern int lm_iter_max;
extern int numffbags;
extern double lm_err_tol;
extern bool lm_vals_bound;
extern double econv;
extern double frqconv;
extern int count_lmin_success;
extern int lmin_success;
const long double pi = 3.14159265358979323846;
const long double econst = 2.71828182845904523536;

extern int    numel;
extern int    max_line_atompar;
extern int    numbty;
extern int    max_line_bondpar;
extern int    numodty;
extern int    max_line_offdpar;
extern int    numaty;
extern int    max_line_angles;
extern int    numtoty;
extern int    max_line_tors;
extern int    numhbty;
extern int    max_line_hbs;
extern vector <vector <int >> inversep; // stores inverse params

class Par { // declaration of a particle
  public:

    Par(); // constructor for class Par of position and velocity components
  ~Par() {

  };

  struct {
   int cycle;
   int iter;
   int parid;
  } state;

  void read_bounds(); // read in min/max bounds of params and set min/max domain. Output dim = number of lines.
  void check_bounds_contff(); // check if parameters are outside bounds and reset them back into domain
  double get_min_dim();
  void read_ffield(); // read ffield file into matrix
  void write_trainset();
  void write_geo();
  void write_ffield_lg(const vector <double> &active_params, int cycle, int iter, int parid);
  void write_ffield(const vector <double> &active_params, int cycle, int iter, int parid);
  void write_aveffield_lg();
  void write_aveffield();

  double get_vel(int n);
  double get_pos(int k);
  vector < double > get_pos_vec();
  double get_bpos(int u);
  vector < double > get_normdir(); // generates a vector (point) uniform on a hypersphere 

  void update_vel(double inertiaf, double CF, vector < double > globpos, double time);
  void update_pos();
  void set_pos(vector < double > pos_best_particle);
  void set_posdim(int i, double posx);
  void set_vel(vector < double > vel_best_particle);
  void set_veldim(int i, double velx);
  void update_bpos();
  void update_pos_levy(vector < double > globpos, double time, double inertiaf);
  double get_levy_McCul(double time, double maxtime); // generate Levy deviate via McCullouch Algorithm

  double get_fitness(); // get particle fitness
  void set_fitness(double fit); // set particle fitness
  double eval_fitness(const vector <double> &active_params, void *mydata); // evalulate fitness
  void iterate();
  double minf;
  vector <double> x;
  vector <double> eval_numgrad(const vector <double> &active_params, void *mydata); // evaluate numerical gradients of fitness
  double get_bfit();
  void set_bfit(double bfit);

  vector < vector < string >> ffieldmat; // ffield inside a vector of vectors (matrix) (for each particle)
  vector < vector < vector <string> > > ffensemble; // ensemble of gbest ffields
  vector < vector < string >> averageffield; // same as above but for ensemble-averaged force field
  vector < int > ffline; // line in ffield file
  vector < int > ffcol; // colum in line in ffield file
  vector < double > mindomain; // min function domain boundary
  vector < double > maxdomain; // max function domain boundary
  vector < double > minpar; // min boundary at initialization step
  vector < double > maxpar; // max boundary at initialization step
  //vector <double> minpos;							// min position value during optimization 
  //vector <double> maxpos;							// max position value during optimization
  int fails; // number of trials particle did not improve personal fitness

  double get_reg(); // calculate regularization
  double reg;
  void dropout(double dropprobability);
  vector <double> numgrad;


  vector < double > pos; // vector of position components of a particle
  vector < double > vel; // vector of velocity components of a particle
  vector < double > bpos; // particle's best own position
  vector <int> dropped_dimns; // dimensions of position to be dropped out during local minimization

  double bfitness; // particle previous best fitness
  double fitness; // particle current fitness
  //double pfitness; // particle temporary fitness when parallelizing training set
};

class Swarm { // declaration of a Swarm
  public:

    Swarm();
  ~Swarm() {

  };

  Par & GetPar(int ParID);
  void get_userinp();
  //void MPI_Split();
  //void MPI_Free();
  void read_icharg_control();
  void AddPar(Par & newPar);
  void MakeEnsembleFF(Swarm & newSwarm);
  void MakeFilesfromPath(Swarm & newSwarm);
  void Populate(Swarm & newSwarm, int iter);
  void Propagate(Swarm & newSwarm, int iter);
  void printdisp(Swarm & newSwarm, int timestep, int iter, int fr);
  void printopt(Swarm & newSwarm, int timestep, int iter, int fr);
  void printpos(Swarm & newSwarm, int timestep, int iter, int fr);
  void printUQFF(Swarm & newSwarm, int timestep, int iter, int fr);
  void printUQQoI(Swarm & newSwarm, int timestep, int iter, int fr);
  void printvel(Swarm & newSwarm, int timestep, int iter, int fr);
  void printdeg(Swarm & newSwarm, int timestep, int iter, int fr);
  void write_ffield_gbest(int core, int cycle, int iter, int par);
  void detovfit(Swarm & newSwram, int cpuid_gbfit, int cycle, int iter, int parID);
  vector < double > get_gbpos();
  vector <double> get_com(Swarm newSwarm);
  double get_disp(Swarm newSwarm); 
  void update_gbpos(Par & newPar);
  double get_gbfit();
  void set_gbfit(double bfit);
  vector < double > gbpos;
  double gbfit;
  int cpuid_gbfit; // CPU with the global best fitness
  int parid_gbfit; // ID of the global best particle
  int get_worse(Swarm newSwarm);
  int get_best(Swarm newSwarm);

  private:

    vector < Par > AllParticles;

  // global best fitness
};

#
#endif
