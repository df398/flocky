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
  boost::filesystem::ofstream log("log.flocky", ofstream::trunc);
#ifdef WITH_MPI
  ierr = MPI_Init( & argc, & argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, & core);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, & numcores);

  if (core == 0) {
    log << "*********************************************************************\n";
    log << "*                __   _                  _                          *\n";
    log << "*               / _| | |                | |                         *\n";
    log << "*              | |_  | |   ___     ___  | | __  _   _               *\n";
    log << "*              |  _| | |  / _ \\   / __| | |/ / | | | |              *\n";
    log << "*              | |   | | | (_) | | (__  |   <  | |_| |              *\n";
    log << "*              |_|   |_|  \\___/   \\___| |_|\\_\\  \\__, |              *\n";
    log << "*                                                __/ |              *\n";
    log << "*                                               |___/               *\n";
    log << "*                                                                   *\n";
    log << "*                                                                   *\n";
    log << "*            f l o c k y  - training reactive force fields          *\n";
    log << "*                            Version 1.0                            *\n";
    log << "*                                                                   *\n";
    log << "*                David Furman, PhD ; df398@cam.ac.uk                *\n";
    log << "*                       University of Cambridge                     *\n";
    log << "*                  Copyright (C) 2018-2019 GPL-3.0                  *\n";
    log << "*                                                                   *\n";
    log << "*       Publications using flocky are requested to cite:            *\n";
    log << "*       Furman, David; Carmeli, Benny; Zeiri, Yehuda;               *\n";
    log << "*       Kosloff, Ronnie, J. Chem. Theory Comput., 2018, 14 (6)      *\n";
    log << "*                                                                   *\n";
    log << "*                    https://www.furmanlab.com                      *\n";
    log << "*                                                                   *\n";
    log << "*********************************************************************\n";
    log.close();

    // check for necessary files
    std::ifstream fin1("ffield");
    if (fin1.fail()) {
      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      log << "Error: unable to open ReaxFF 'ffield' file on CPU " << core << ". Aborting! \n";
      log.close();
      exit(EXIT_FAILURE);
    };

    std::ifstream fin2("control");
    if (fin2.fail()) {
      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      log << "Error: unable to open 'control' file on CPU " << core << ". Aborting! \n";
      log.close();
      exit(EXIT_FAILURE);
    };

    std::ifstream fin3("geo");
    if (fin3.fail()) {
      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      log << "Error: unable to open geo file on CPU " << core << ". Aborting! \n";
      log.close();
      exit(EXIT_FAILURE);
    };

    std::ifstream fin4("params.mod");
    if (fin4.fail()) {
      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      log << "Error: unable to open 'params.mod' file on CPU " << core << ". Aborting! \n";
      log.close();
      exit(EXIT_FAILURE);
    };
  }; // end opening message and checks on core 0

  for (cycle = 0; cycle < maxcycles; cycle++) {
    string cyclecount = std::to_string(cycle);

    // Start collecting timing info
    // ------------------------------------------------------------------------//
    double starttime, endtime, diff, sumdiff, avgdiff;
    starttime = MPI_Wtime();
    // ------------------------------------------------------------------------//

    // bkup original ffield
    boost::filesystem::copy_file("ffield", "ffield.initial." + cyclecount, 
      boost::filesystem::copy_option::overwrite_if_exists);
    
    Swarm MySwarm;
    MySwarm.get_userinp();
    MySwarm.Populate(MySwarm, cycle);
    MySwarm.Propagate(MySwarm, cycle);

    // End collecting timing info
    // ------------------------------------------------------------------------//
    endtime = MPI_Wtime();
    diff = endtime - starttime;   
    MPI_Reduce( & diff, & sumdiff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    avgdiff = sumdiff / numcores;
    if (core == 0) {
      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      log << "Total CPU time: " << boost::format("%6.1f") %avgdiff << " seconds" << endl;
      log.close();
    };
    // ------------------------------------------------------------------------//

    MySwarm = Swarm();

  }; // end cycles

  ierr = MPI_Finalize();
  return 0;
};
#endif

#ifndef WITH_MPI
    log << "*********************************************************************\n";
    log << "*                __   _                  _                          *\n";
    log << "*               / _| | |                | |                         *\n";
    log << "*              | |_  | |   ___     ___  | | __  _   _               *\n";
    log << "*              |  _| | |  / _ \\   / __| | |/ / | | | |              *\n";
    log << "*              | |   | | | (_) | | (__  |   <  | |_| |              *\n";
    log << "*              |_|   |_|  \\___/   \\___| |_|\\_\\  \\__, |              *\n";
    log << "*                                                __/ |              *\n";
    log << "*                                               |___/               *\n";
    log << "*                                                                   *\n";
    log << "*                                                                   *\n";
    log << "*            f l o c k y  - training reactive force fields          *\n";
    log << "*                            Version 1.0                            *\n";
    log << "*                                                                   *\n";
    log << "*                David Furman, PhD ; df398@cam.ac.uk                *\n";
    log << "*                       University of Cambridge                     *\n";
    log << "*                  Copyright (C) 2018-2019 GPL-3.0                  *\n";
    log << "*                                                                   *\n";
    log << "*       Publications using flocky are requested to cite:            *\n";
    log << "*       Furman, David; Carmeli, Benny; Zeiri, Yehuda;               *\n";
    log << "*       Kosloff, Ronnie, J. Chem. Theory Comput., 2018, 14 (6)      *\n";
    log << "*                                                                   *\n";
    log << "*                    https://www.furmanlab.com                      *\n";
    log << "*                                                                   *\n";
    log << "*********************************************************************\n";
    log.close();


  for (cycle = 0; cycle < maxcycles; cycle++) {
    string cyclecount = std::to_string(cycle);

    auto starttime = chrono::steady_clock::now();

    // bkup original ffield before eval_fitness generates one for current particle fitness evaluation
    boost::filesystem::copy_file("ffield", "ffield.initial." + cyclecount,
    boost::filesystem::copy_option::overwrite_if_exists);

    Swarm MySwarm;
    MySwarm.get_userinp();
    MySwarm.Populate(MySwarm, cycle);
    MySwarm.Propagate(MySwarm, cycle);
    
    auto endtime = chrono::steady_clock::now();
    auto diff = endtime - starttime;

    boost::filesystem::ofstream log("log.flocky", ofstream::app);
    log << " " << endl;
    log << "Total CPU time: " << chrono::duration_cast<sec>(diff).count() << " seconds" << endl;
    MySwarm = Swarm();
    log.close();
  };
  return 0;
};
#endif
