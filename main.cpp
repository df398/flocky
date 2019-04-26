/*  ---------------------------------------------------------------------- *
    flocky v1.0 Copyright (C) 2019 David Furman, PhD. 
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

#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include "par.h"

#ifndef WITH_MPI
#include <chrono>
using sec = chrono::seconds;
#endif

using namespace std;

int main(int argc, char * argv[]) {
#ifdef WITH_MPI
  ierr = MPI_Init( & argc, & argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, & core);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, & numcores);

  if (core == 0) {
    cout << "*************************************************************************\n";
    cout << "*                                                                       *\n";
    cout << "*        flocky v1.0 Copyright (C) 2019, GNU L-GPL 3.0                  *\n";
    cout << "*              David Furman, PhD, df398@cam.ac.uk                       *\n";
    cout << "*                 University of Cambridge, UK.                          *\n";
    cout << "*                                                                       *\n";
    cout << "*        Publications using flocky should cite:                         *\n";
    cout << "*        Furman, David; Carmeli, Benny; Zeiri, Yehuda; Kosloff,         *\n";
    cout << "*        Ronnie, J. Chem. Theory Comput., 2018, 14 (6)                  *\n";
    cout << "*                                                                       *\n";
    cout << "*                       https://www.furmanlab.com                       *\n";
    cout << "*                                                                       *\n";
    cout << "*************************************************************************\n";

    // check for necessary files
    std::ifstream fin1("ffield");
    if (fin1.fail()) {
      cout << "Unable to open ReaxFF 'ffield' file on CPU " << core << ". Aborting! \n";
      exit(EXIT_FAILURE);
    }

    std::ifstream fin2("control");
    if (fin2.fail()) {
      cout << "Unable to open 'control' file on CPU " << core << ". Aborting! \n";
      exit(EXIT_FAILURE);
    }

    std::ifstream fin3("geo");
    if (fin3.fail()) {
      cout << "Unable to open geo file on CPU " << core << ". Aborting! \n";
      exit(EXIT_FAILURE);
    }

    std::ifstream fin4("params.mod");
    if (fin4.fail()) {
      cout << "Unable to open 'params.mod' file on CPU " << core << ". Aborting! \n";
      exit(EXIT_FAILURE);
    }

    // get optimization settings from user
    cout << "\n";
    cout << "Is ffield file of the ReaxFF_lg type? [1=true/0=false]" << endl;
    cin >> lg_yn;

    if (lg_yn == true) {
      cout << "\n";
      cout << "Proceeding with LG-augmented ReaxFF ffield. " << endl;
      cout << "Important: make sure your params.mod file is adjusted accordingly!\n" << endl;
    } else {
      cout << "\n";
      cout << "Proceeding with standard ReaxFF ffield." << endl;
      cout << "Important: make sure your params.mod file is adjusted accordingly!\n" << endl;
    };

    cout << "Take initial positions of first search agent from current ffield file? [1=true/0=false]" << endl;
    cout << "Note: Initial positions of all other agents will be random." << endl;
    cin >> contff;

    cout << "Use percent change in parameter bounds? [1=true/0=false: use bounds from params.mod]" << endl;
    cin >> perc_yn;
    if (perc_yn == true){
      cout << "Enter %change [0.0 - 1.0]: " << endl;
      cin >> perc;
    };

    while (NumP <= 0) {
      cout << "Enter swarm size to distribute over " << numcores << " processors:" << endl;
      cin >> NumP;
      NumP = int(floor(NumP / numcores)); // each CPU gets its own allocation of swarm members
    };
    cout << "Change flocky optimization defaults? [1=true/0=false]" << endl;
    cin >> chang;

    if (chang == true) {
      cout << "New c1 value:" << endl;
      cin >> c1;
      cout << "New c2 value:" << endl;
      cin >> c2;
      cout << "New w1 value:" << endl;
      cin >> inertiamax;
      cout << "New w2 value (w2 >= w1):" << endl;
      cin >> inertiamin;
      cout << "New fail_i value (int > 0): " << endl;
      cin >> faili;
      cout << "New gamma value (> 0): " << endl;
      cin >> levyscale;
      cout << "\n";
    } else {
      cout << "Proceeding with default optimization settings:" << endl;
      cout <<"c1 = c2 = 2.0; w1 = 0.9, w2 = 0.4, fail_i = 1, gamma = 0.01" << endl;
    };

    //cout << "Enter levy scale:" << endl;
    //cin >> levyscale;
    cout << "\n";
    cout << "Enter max number of iterations:" << endl;
    cin >> maxiters;
    cout << "Enter output frequency:" << endl;
    cin >> freq;
    cout << "Enter number of optimization cycles:" << endl;
    cin >> maxcycles;
  };
  // prepare directories for each CPU process
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core);
  // get present working directory full path
  boost::filesystem::path pwd(boost::filesystem::current_path());

  boost::filesystem::copy_file(pwd.string() + "/ffield", pwd.string() + "/CPU." + str_core + "/ffield",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/control", pwd.string() + "/CPU." + str_core + "/control",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/geo", pwd.string() + "/CPU." + str_core + "/geo",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/CPU." + str_core + "/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/params.mod", pwd.string() + "/CPU." + str_core + "/params.mod",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/trainset.in", pwd.string() + "/CPU." + str_core + "/trainset.in",
    boost::filesystem::copy_option::overwrite_if_exists);
  // check if reaxff runs with fixed charges (= require charges file)
  boost::filesystem::ifstream chargeStream("charges");
  if (!chargeStream.fail()) {
    boost::filesystem::copy_file(pwd.string() + "/charges", pwd.string() + "/CPU." + str_core + "/charges",
      boost::filesystem::copy_option::overwrite_if_exists);
   fixcharges = true;
  };
  chargeStream.close();

  boost::filesystem::ofstream iopt_file(pwd.string() + "/CPU." + str_core + "/fort.20");
  iopt_file << "0";
  iopt_file.close();
  ofstream outfile35(pwd.string() + "/CPU." + str_core + "/fort.35");
  outfile35 << "23434.1" << endl;
  outfile35.close();

  // broadcast between all processes from process 0
  // MPI_Bcast must be visible to all processes!
  MPI_Bcast( & lg_yn, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & chang, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & contff, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & NumP, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & c1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & c2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & inertiamax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & inertiamin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & faili, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & levyscale, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & maxiters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & freq, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & maxcycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & perc_yn, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & perc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & fixcharges, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

  for (cycle = 0; cycle < maxcycles; cycle++) {
    string cyclecount = std::to_string(cycle);
    double starttime, endtime, diff, sumdiff, avgdiff;
    starttime = MPI_Wtime();
    // bkup original ffield before eval_fitness generates one for current particle fitness evaluation
    boost::filesystem::copy_file("ffield", "ffield.initial." + cyclecount, 
    boost::filesystem::copy_option::overwrite_if_exists);

    Swarm MySwarm;
    MySwarm.Populate(MySwarm, cycle);
    MySwarm.Propagate(MySwarm, cycle);
    // Each core should time its own operation, then send its data to core 0. Then, core 0 receives
    // timings from all cores, adds them up and divides by the sizecpus. Then prints this to the log
    // cpu file in the cpu folder. This will be the average CPU time of all CPUS.
    endtime = MPI_Wtime();
    diff = endtime - starttime;
    MPI_Reduce( & diff, & sumdiff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    avgdiff = sumdiff / numcores;

    if (core == 0) {
      cout << "Total CPU time: " << avgdiff << " seconds" << endl;
    };
    MySwarm = Swarm();

  };
  ierr = MPI_Finalize();
  return 0;
};
#endif

#ifndef WITH_MPI
    cout << "*************************************************************************\n";
    cout << "*                                                                       *\n";
    cout << "*        flocky v1.0 Copyright (C) 2019, GNU L-GPL 3.0                  *\n";
    cout << "*              David Furman, PhD, df398@cam.ac.uk                       *\n";
    cout << "*                 University of Cambridge, UK.                          *\n";
    cout << "*                                                                       *\n";
    cout << "*        Publications using flocky should cite:                         *\n";
    cout << "*        Furman, David; Carmeli, Benny; Zeiri, Yehuda; Kosloff,         *\n";
    cout << "*        Ronnie, J. Chem. Theory Comput., 2018, 14 (6)                  *\n";
    cout << "*                                                                       *\n";
    cout << "*                       https://www.furmanlab.com                       *\n";
    cout << "*                                                                       *\n";
    cout << "*************************************************************************\n";

    // check for necessary files
    std::ifstream fin1("ffield");
    if (fin1.fail()) {
      cout << "Unable to open ReaxFF 'ffield' file. Aborting! \n";
      exit(EXIT_FAILURE);
    }

    std::ifstream fin2("control");
    if (fin2.fail()) {
      cout << "Unable to open 'control' file. Aborting! \n";
      exit(EXIT_FAILURE);
    }

    std::ifstream fin3("geo");
    if (fin3.fail()) {
      cout << "Unable to open geo file on CPU. Aborting! \n";
      exit(EXIT_FAILURE);
    }

    std::ifstream fin4("params.mod");
    if (fin4.fail()) {
      cout << "Unable to open 'params.mod' file. Aborting! \n";
      exit(EXIT_FAILURE);
    }
    // get optimization settings from user
    cout << "\n";
    cout << "Is ffield file of the ReaxFF_lg type? [1=true/0=false]" << endl;
    cin >> lg_yn;

    if (lg_yn == true) {
      cout << "\n";
      cout << "Proceeding with LG-augmented ReaxFF ffield. " << endl;
      cout << "Important: make sure your params.mod file is adjusted accordingly!\n" << endl;
    } else {
      cout << "\n";
      cout << "Proceeding with standard ReaxFF ffield." << endl;
      cout << "Important: make sure your params.mod file is adjusted accordingly!\n" << endl;
    };

    cout << "Take initial positions of first search agent from current ffield file? [1=true/0=false]" << endl;
    cout << "Note: Initial positions of all other agents will be random." << endl;
    cin >> contff;

    cout << "Use percent change in parameter bounds? [1=true/0=false: use bounds from params.mod]" << endl;
    cin >> perc_yn;
    if (perc_yn == true){
      cout << "Enter %change [0.0 - 1.0]: " << endl;
      cin >> perc;
    };
    while (NumP <= 0) {
      cout << "Enter swarm size (recommended 4-20 for serial job):" << endl;
      cin >> NumP;
    };
    cout << "Change flocky optimization defaults? [1=true/0=false]" << endl;
    cin >> chang;

    if (chang == true) {
      cout << "New c1 value:" << endl;
      cin >> c1;
      cout << "New c2 value:" << endl;
      cin >> c2;
      cout << "New w1 value:" << endl;
      cin >> inertiamax;
      cout << "New w2 value (w2 >= w1):" << endl;
      cin >> inertiamin;
      cout << "New fail_i value (int > 0): " << endl;
      cin >> faili;
      cout << "New gamma value (> 0): " << endl;
      cin >> levyscale;
      cout << "\n";
    } else {
      cout << "Proceeding with default optimization settings:" << endl;
      cout <<"c1 = c2 = 2.0; w1 = 0.9, w2 = 0.4, fail_i = 1, gamma = 0.01" << endl;
    };

    //cout << "Enter levy scale:" << endl;
    //cin >> levyscale;
    cout << "\n";
    cout << "Enter max number of iterations:" << endl;
    cin >> maxiters;
    cout << "Enter output frequency:" << endl;
    cin >> freq;
    cout << "Enter number of optimization cycles:" << endl;
    cin >> maxcycles;
  
  // check if reaxff runs with fixed charges (= require charges file)
  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::ifstream chargeStream("charges");
  if (!chargeStream.fail()) {
    boost::filesystem::copy_file(pwd.string() + "/charges", pwd.string() + "/charges",
      boost::filesystem::copy_option::overwrite_if_exists);
   fixcharges = true;
  };
  chargeStream.close();

  boost::filesystem::ofstream iopt_file(pwd.string() + "/fort.20");
  iopt_file << "0";
  iopt_file.close();
  ofstream outfile35(pwd.string() + "/fort.35");
  outfile35 << "23434.1" << endl;
  outfile35.close();
  for (cycle = 0; cycle < maxcycles; cycle++) {
    string cyclecount = std::to_string(cycle);

    auto starttime = chrono::steady_clock::now();

    // bkup original ffield before eval_fitness generates one for current particle fitness evaluation
    boost::filesystem::copy_file("ffield", "ffield.initial." + cyclecount,
    boost::filesystem::copy_option::overwrite_if_exists);

    Swarm MySwarm;
    MySwarm.Populate(MySwarm, cycle);
    MySwarm.Propagate(MySwarm, cycle);
    
    auto endtime = chrono::steady_clock::now();
    auto diff = endtime - starttime;
    cout << " " << endl;
    cout << "Total CPU time: " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;
    MySwarm = Swarm();

  };
  return 0;
};
#endif
