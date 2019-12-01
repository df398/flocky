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

#include "par.h"

// initialize non-determinstic seed
std::random_device seed;
// Mersene Twister: Good quality random number generator
std::mt19937 generator(seed() + core);
std::uniform_real_distribution < double > dist1(0.0, 1.0);

int ierr, core, numcores = 0;
int funceval = 0;
double curfit;
// tag if particle position in dimension lies inside [mindomain, maxdomain]
//bool inside = true;
int dim = 0;
int NumP = 0;
int freq = 1;
int cycle = 0;
int maxcycles = 1;
int iter = 0;
int maxiters = 0;
int parid_gbfit = 0;
bool fixcharges = false;
bool lg_yn = false;
bool contff = false;
bool chang = false;
bool perc_yn = true;
bool ofit = false;
bool uq = false;
bool gbfitfound = false;
bool firstovfit = true;
double ovfitness = 0.0;
double currovfitness = 0.0;
double initial_disp = 0.0;
double perc = 0.2;
double c1 = 2.0;
double c2 = 2.0;
double inertiafac = 0.9;
double inertiamax = 0.9;
double inertiamin = 0.4;
double levyscale = 1.0;
double confac = 1.0;
int faili = 1;

// --------------------------------------------------------------- //
//             General functions definitions                       //
// --------------------------------------------------------------- //

// define int of modulo
int mod(int x, int m) {
  return (x % m + m) % m;
}


// return L2-norm (magnitude) of a vector. uses a range-based loop.
double l2_norm(vector < double > const & u) {
  double accum = 0.;
  for (double x: u) {
    accum += x * x;
  }
  return sqrt(accum);
};

// --------------- End general functions definitions -------------- //


Par::Par() {

  // read ffield file into matrix ffieldmat. split by each entry.
  read_ffield();

  // set min/max domains from params.mod file and set dim = numlines in params.mod file
  read_bounds();

  for (int m = 0; m < dim; m++) {
    std::uniform_real_distribution < double > dist2(mindomain.at(m), maxdomain.at(m));
    //	std::normal_distribution <double> distnorm(0.0,20.0);
    double x = dist2(generator);
    double v = 0.5 * (dist2(generator) - x);
    // initialize particle's position vector
    pos.push_back(x);
    // initialize particle's velocity vector
    vel.push_back(v);
    // initialize particle's best own position vector
    bpos = pos;
  };
};

void Par::read_ffield() {
  // read ffield file into ffieldmat matrix split by entries
  // following: https://stackoverflow.com/questions/10521658/reading-files-columns-into-array
  ffieldmat.clear();
#ifdef WITH_MPI
  string str_core = std::to_string(core);
  std::ifstream fin("CPU." + str_core + "/ffield");
#endif
#ifndef WITH_MPI
  std::ifstream fin("ffield"); 
#endif
  if (fin.fail()) {
#ifdef WITH_MPI
    cout << "Error: unable to open 'ffield' file on CPU" << core << " \n";
#endif
#ifndef WITH_MPI
    cout << "Unable to open 'ffield' file \n"; 
#endif
    fin.close();
    exit(EXIT_FAILURE);
  } else {
    std::string line;
    int numlines = 0;

    // for each line
    while (std::getline(fin, line)) {
      numlines++;
      // create a new row
      std::vector < string > lineData;
      string val;
      std::istringstream lineStream(line);
      // for each value in line
      while (lineStream >> val) {
        // add to the current row
        lineData.push_back(val);
      };
      // add row to ffieldmat matrix
      ffieldmat.push_back(lineData);

    };
    fin.close();
  };
};

void Par::write_ffield(int cycle, int iter, int par) {
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(par);
  // update particle ffield
  int index = 0;
  // for each line value in ffline and corresponding column value
  // in ffcol cp corresponding pos value into ffieldmat. Range-based for
  // statement --> Following: https://msdn.microsoft.com/en-us/library/jj203382.aspx
  for (int line: ffline) {
    // buffer for sprintf
    char buffer[50];
    double n;
    n = sprintf(buffer, "%9.4f", pos.at(index));
    ffieldmat[line][ffcol.at(index)] = buffer;
    index = index + 1;
  };
  string str_core = std::to_string(core);
  // write updated ffield to ffield.tmp.cycle.iter.parid file
  ofstream output_file;
  // current ffield file stream
  ifstream ffield_file;
  string comment;
#ifdef WITH_MPI
  output_file.open("CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID, ios::out);
  ffield_file.open("CPU." + str_core + "/ffield", ios:: in );
#endif
#ifndef WITH_MPI
  output_file.open("ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID, ios::out);
  ffield_file.open("ffield", ios:: in );
#endif
  // read all ffield file into vector of lines to use it next to write comments/header lines
  vector < string > ffield_lines;
  while (getline(ffield_file, comment)) {
    ffield_lines.push_back(comment);
  };
  ffield_file.close();
  /* -------------------------------------
   * WRITE HEADER LINE FOR FFIELD
   * -------------------------------------
   */
  for (int m = 0; m < 2; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* -------------------------------------
   * WRITE GENERAL PARAMS SECTION
   * -------------------------------------
   */

  for (int m = 2; m < 41; m++) {
    boost::format f("%10.4f %s");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };
  /* ------------------------------------
   * WRITE ATOM PARAMS HEADERS
   * ------------------------------------
   */

  for (int m = 41; m < 45; m++) {
    output_file << ffield_lines.at(m) << endl;
  };
  
  using boost::is_any_of;
  // the line that contains the number of elements
  string numel_line = ffield_lines.at(41);
  // vector to store the words after split
  vector < string > results;
  boost::trim(numel_line);
  boost::split(results, numel_line, is_any_of("\t "));
  int numel = stoi(results.at(0));
  int max_line_atompar = 4 * numel + 45;
  /* -------------------------------------
   * WRITE ATOM PARAMS SECTION
   * -------------------------------------
   */

  int m = 45;
  while (m < max_line_atompar) {
    boost::format f("% s%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions( f.exceptions() &
    ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
    for (std::vector<std::string>::iterator it=ffieldmat.at(m).begin();it!=ffieldmat.at(m).end();++it){
            f = f % (*it);
    };
    output_file << f << endl;
    for (int i = 1; i < 4; i++){
      m = m + 1;
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions( f.exceptions() &
      ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
      for (std::vector<std::string>::iterator it=ffieldmat.at(m).begin();it!=ffieldmat.at(m).end();++it){
              f = f % (*it);
      };
      output_file << f << endl;
    };
    m = m + 1;
  };

  /* ------------------------------------
   * WRITE BONDS PARAMS HEADERS
   * ------------------------------------
   */

  for (int m = max_line_atompar; m < max_line_atompar + 2; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  using boost::is_any_of;
  // the line the contains the number of bond types
  string numbty_line = ffield_lines.at(max_line_atompar);
  // vector to store the words after split
  vector < string > results_bonds;
  boost::trim(numbty_line);
  boost::split(results_bonds, numbty_line, is_any_of("\t "));
  int numbty = stoi(results_bonds.at(0));
  int max_line_bondpar = 2 * numbty + max_line_atompar;

  /* ------------------------------------
   * WRITE BONDS PARAMS SECTION
   * ------------------------------------
   */

  for (int m = max_line_atompar + 2; m < max_line_bondpar + 2; m++) {
    if (mod(m, 2) != 0) {
      boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    } else {
      boost::format f("      %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  /* ------------------------------------
   * WRITE OFF-DIAG HEADERS
   * ------------------------------------
   */

  for (int m = max_line_bondpar + 2; m < max_line_bondpar + 3; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE OFF-DIAG PARAMS SECTION
   * ------------------------------------
   */

  using boost::is_any_of;
  // the line that contains the number of off-diag types
  string numodty_line = ffield_lines.at(max_line_bondpar + 2);
  // vector to store the words after split
  vector < string > results_offdiag;
  boost::trim(numodty_line);
  boost::split(results_offdiag, numodty_line, is_any_of("\t "));
  int numodty = stoi(results_offdiag.at(0));
  int max_line_offdpar = max_line_bondpar + 3 + numodty;

  for (int m = max_line_bondpar + 3; m < max_line_offdpar; m++) {
    // the last entry is 10.4f because dispersion coeff. can get to 4 digits long.
    // So, to prevent it sticking to the left column
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };
  /* ------------------------------------
   * WRITE ANGLE HEADERS
   * ------------------------------------
   */

  for (int m = max_line_offdpar; m < max_line_offdpar + 1; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE ANGLE PARAMS SECTION
   * ------------------------------------
   */

  using boost::is_any_of;
  // the line that contains the number of angle types
  string numaty_line = ffield_lines.at(max_line_offdpar);
  // vector to store the words after split
  vector < string > results_angle;
  boost::trim(numaty_line);
  boost::split(results_angle, numaty_line, is_any_of("\t "));
  int numaty = stoi(results_angle.at(0));
  int max_line_angles = max_line_offdpar + 1 + numaty;

  for (int m = max_line_offdpar + 1; m < max_line_angles; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };
  /* ------------------------------------
   * WRITE TORSION HEADERS
   * ------------------------------------
   */

  for (int m = max_line_angles; m < max_line_angles + 1; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE TORSION PARAMS SECTION
   * ------------------------------------
   */

  using boost::is_any_of;
  // the line that contains the number of torsion types
  string numtoty_line = ffield_lines.at(max_line_angles);
  // vector to store the words after split
  vector < string > results_tors;
  boost::trim(numtoty_line);
  boost::split(results_tors, numtoty_line, is_any_of("\t "));
  int numtoty = stoi(results_tors.at(0));
  int max_line_tors = max_line_angles + 1 + numtoty;

  for (int m = max_line_angles + 1; m < max_line_tors; m++) {
    boost::format f("  %i  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };
  /* ------------------------------------
   * WRITE H-BOND HEADERS
   * ------------------------------------
   */

  for (int m = max_line_tors; m < max_line_tors + 1; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE H-BOND PARAMS SECTION
   * ------------------------------------
   */

  using boost::is_any_of;
  // the line that contains the number of Hbond types
  string numhbty_line = ffield_lines.at(max_line_tors);
  // vector to store the words after split
  vector < string > results_hb;
  boost::trim(numhbty_line);
  boost::split(results_hb, numhbty_line, is_any_of("\t "));
  int numhbty = stoi(results_hb.at(0));
  int max_line_hbs = max_line_tors + 1 + numhbty;

  for (int m = max_line_tors + 1; m < max_line_hbs; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  output_file.close();

  // replace ffield file with the new ffield (ffield.tmp.cycle.iter.parid)
  boost::filesystem::path pwd(boost::filesystem::current_path());
#ifdef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    pwd.string() + "/CPU." + str_core + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
#endif
#ifndef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    pwd.string() + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
#endif
};





void Par::write_ffield_lg(int cycle, int iter, int par) {
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(par);

  // update particle ffield
  int index = 0;
  // for each line value in ffline and corresponding column value in ffcol 
  // cp corresponding pos value into ffieldmat. Range-based for statement 
  // following: https://msdn.microsoft.com/en-us/library/jj203382.aspx
  for (int line: ffline) {
    // buffer for sprintf
    char buffer[50];
    double n;
    n = sprintf(buffer, "%9.4f", pos.at(index));
    ffieldmat[line][ffcol.at(index)] = buffer;
    index = index + 1;

  };
  string str_core = std::to_string(core);
  // write updated ffield to ffield file (ffield.tmp.cycle.iter.parid)
  ofstream output_file;
  // current ffield file stream
  ifstream ffield_file;
  string comment;

  output_file.open("CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID, ios::out);
  ffield_file.open("CPU." + str_core + "/ffield", ios:: in );
  // read all ffield file into vector of lines to use it next 
  // to write comments/header lines
  vector < string > ffield_lines;
  while (getline(ffield_file, comment)) {
    ffield_lines.push_back(comment);
  };
  ffield_file.close();

  /* -------------------------------------
   * WRITE HEADER LINE FOR FFIELD
   * -------------------------------------
   */

  for (int m = 0; m < 2; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* -------------------------------------
   * WRITE GENERAL PARAMS SECTION
   * -------------------------------------
   */

  for (int m = 2; m < 41; m++) {
    boost::format f("%10.4f %s");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  /* ------------------------------------
   * WRITE ATOM PARAMS HEADERS
   * ------------------------------------
   */

  for (int m = 41; m < 45; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  using boost::is_any_of;
  // the line that contains the number of elements
  string numel_line = ffield_lines.at(41);
  // vector to store the words after split
  vector < string > results;
  boost::trim(numel_line);
  boost::split(results, numel_line, is_any_of("\t "));
  int numel = stoi(results.at(0));
  int max_line_atompar = 5 * numel + 45;

  /* -------------------------------------
   * WRITE ATOM PARAMS SECTION
   * -------------------------------------
   */

  int m = 45;
  while (m < max_line_atompar) {
    boost::format f("% s%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions( f.exceptions() &
    ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
    for (std::vector<std::string>::iterator it=ffieldmat.at(m).begin();it!=ffieldmat.at(m).end();++it){
            f = f % (*it);
    };
    output_file << f << endl;
    
    for (int i = 1; i < 4; i++){
      m = m + 1;
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions( f.exceptions() &
      ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
      for (std::vector<std::string>::iterator it=ffieldmat.at(m).begin();it!=ffieldmat.at(m).end();++it){
              f = f % (*it);
      };
      output_file << f << endl;
    };
    
    m = m +1;
    boost::format f2("   %9.4f%9.4f");
    f2.exceptions( f2.exceptions() &
    ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
    for (std::vector<std::string>::iterator it=ffieldmat.at(m).begin();it!=ffieldmat.at(m).end();++it){
            f2 = f2 % (*it);
    };
    output_file << f2 << endl;


    m = m + 1;
  };

  /* ------------------------------------
   * WRITE BONDS PARAMS HEADERS
   * ------------------------------------
   */

  for (int m = max_line_atompar; m < max_line_atompar + 2; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* 
    06.03.19 Check how many elements are there in the ffield then limit the write of bond params
    to only those lines that belong to number of bond types. This is to avoid writing lines
    that belong to the next section (off-dia) in the bond section.
  */
  using boost::is_any_of;
  // the line the contains the number of bond types
  string numbty_line = ffield_lines.at(max_line_atompar);
  // vector to store the words after split
  vector < string > results_bonds;
  boost::trim(numbty_line);
  boost::split(results_bonds, numbty_line, is_any_of("\t "));
  int numbty = stoi(results_bonds.at(0));
  int max_line_bondpar = 2 * numbty + max_line_atompar;

  /* ------------------------------------
   * WRITE BONDS PARAMS SECTION
   * ------------------------------------
   */

  for (int m = max_line_atompar + 2; m < max_line_bondpar + 2; m++) {
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };

    output_file << f << endl;
    m = m + 1;
    boost::format f2("      %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f2.exceptions(f2.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f2 = f2 % ( * it);
    };

    output_file << f2 << endl;
  };

  /* ------------------------------------
   * WRITE OFF-DIAG HEADERS
   * ------------------------------------
   */

  for (int m = max_line_bondpar + 2; m < max_line_bondpar + 3; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE OFF-DIAG PARAMS SECTION
   * ------------------------------------
   */

  /* 
    06.03.19 Check how many elements are there in the ffield then limit the write of bond params
    to only those lines that belong to number of bond types. This is to avoid writing lines
    that belong to the next section (off-dia) in the bond section.
  */
  using boost::is_any_of;
  // the line that contains the number of off-diag types
  string numodty_line = ffield_lines.at(max_line_bondpar + 2);
  // vector to store the words after split
  vector < string > results_offdiag; 
  boost::trim(numodty_line);
  boost::split(results_offdiag, numodty_line, is_any_of("\t "));
  int numodty = stoi(results_offdiag.at(0));
  int max_line_offdpar = max_line_bondpar + 3 + numodty;

  for (int m = max_line_bondpar + 3; m < max_line_offdpar; m++) {
    // the last entry is 10.4f because dispersion coeff. can get to 4 digits long.
    // So, to prevent it sticking to the left column
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%10.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  /* ------------------------------------
   * WRITE ANGLE HEADERS
   * ------------------------------------
   */

  for (int m = max_line_offdpar; m < max_line_offdpar + 1; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE ANGLE PARAMS SECTION
   * ------------------------------------
   */

  /* 
    06.03.19 Check how many angles are there in the ffield then limit the write of angle params
    to only those lines that belong to number of angle types. This is to avoid writing lines
    that belong to the next section (torsions) in the angles section.
  */
  using boost::is_any_of;
  // the line that contains the number of angle types
  string numaty_line = ffield_lines.at(max_line_offdpar);
  // vector to store the words after split
  vector < string > results_angle;
  boost::trim(numaty_line);
  boost::split(results_angle, numaty_line, is_any_of("\t "));
  int numaty = stoi(results_angle.at(0));
  int max_line_angles = max_line_offdpar + 1 + numaty;

  for (int m = max_line_offdpar + 1; m < max_line_angles; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  /* ------------------------------------
   * WRITE TORSION HEADERS
   * ------------------------------------
   */

  for (int m = max_line_angles; m < max_line_angles + 1; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE TORSION PARAMS SECTION
   * ------------------------------------
   */

  /* 
    06.03.19 Check how many torsion types are there in the ffield then limit the write of torsion params
    to only those lines that belong to number of torsion types. This is to avoid writing lines
    that belong to the next section (h-bonds) in the torsions section.
  */
  using boost::is_any_of;
  // the line that contains the number of torsion types
  string numtoty_line = ffield_lines.at(max_line_angles);
  // vector to store the words after split
  vector < string > results_tors;
  boost::trim(numtoty_line);
  boost::split(results_tors, numtoty_line, is_any_of("\t "));
  int numtoty = stoi(results_tors.at(0));
  int max_line_tors = max_line_angles + 1 + numtoty;

  for (int m = max_line_angles + 1; m < max_line_tors; m++) {
    boost::format f("  %i  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  /* ------------------------------------
   * WRITE H-BOND HEADERS
   * ------------------------------------
   */

  for (int m = max_line_tors; m < max_line_tors + 1; m++) {
    output_file << ffield_lines.at(m) << endl;
  };

  /* ------------------------------------
   * WRITE H-BOND PARAMS SECTION
   * ------------------------------------
   */

  /* 
    06.03.19 Check how many Hbond types are there in the ffield then limit the write of Hbond params
    to only those lines that belong to number of Hbond types. This is to avoid accessing nonexistent cells.
  */
  using boost::is_any_of;
  // the line that contains the number of Hbond types
  string numhbty_line = ffield_lines.at(max_line_tors);
  // vector to store the words after split
  vector < string > results_hb;
  boost::trim(numhbty_line);
  boost::split(results_hb, numhbty_line, is_any_of("\t "));
  int numhbty = stoi(results_hb.at(0));
  int max_line_hbs = max_line_tors + 1 + numhbty;

  for (int m = max_line_tors + 1; m < max_line_hbs; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  output_file.close();

  // replace ffield file with the new ffield (ffield.tmp.cycle.iter.parid)
  boost::filesystem::path pwd(boost::filesystem::current_path());
#ifdef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    pwd.string() + "/CPU." + str_core + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp.*");
#endif
#ifndef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    pwd.string() + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  //boost::filesystem::remove(pwd.string() + "/ffield.tmp.*");
#endif

};

void Par::read_bounds() {
  // Reading params columns into allData array split by entries
  // following: https://stackoverflow.com/questions/10521658/reading-files-columns-into-array

  std::vector < std::vector < double > > allData;
  string str_core = std::to_string(core);
  std::ifstream fin("params.mod");
  if (fin.fail()) {
    cout << "Unable to open parameters file 'params.mod' on CPU " << core << ". \n";
    fin.close();
    exit(EXIT_FAILURE);
  } else {
    std::string line;
    int numlines = 0;

    // for each line
    while (std::getline(fin, line)) {
      numlines++;
      // create a new row
      std::vector < double > lineData;
      double val;
      std::istringstream lineStream(line);
      // for each value in line
      while (lineStream >> val) {
        // add to the current row
        lineData.push_back(val);
      };
      // add row to allData
      allData.push_back(lineData);
    };
    fin.close();

    dim = numlines;

    if (perc_yn == false){
      for (int i=0; i<dim; i++){
        // read bounds from modified params file
        mindomain.push_back(allData[i][5]);
        maxdomain.push_back(allData[i][6]);
        // read line number of parameter from modified params file
        ffline.push_back(allData[i][0]-1);
        // read column number of parameter from modified params file
        ffcol.push_back(allData[i][1]-1);
      };
    } else {
      for (int i=0; i<dim; i++){
        // read line number of parameter from modified params file
        ffline.push_back(allData[i][0]-1);
        // read column number of parameter from modified params file
        ffcol.push_back(allData[i][1]-1);
        mindomain.push_back(stod(ffieldmat.at(ffline.at(i)).at(ffcol.at(i))));
        maxdomain.push_back(stod(ffieldmat.at(ffline.at(i)).at(ffcol.at(i))));
        mindomain.at(i) = (1.0 - perc)*mindomain.at(i);
        maxdomain.at(i) + (1.0 + perc)*maxdomain.at(i);
        };

    };


  };
};

double Par::get_min_dim() {
  // get length of smallest domain dimension 
  double min_dim = 1.0E99;
  for (int i = 0; i < dim; i++) {
    if (min_dim > (maxdomain.at(i) - mindomain.at(i))) {
      min_dim = maxdomain.at(i) - mindomain.at(i);
    };
  };

  return min_dim;
};

double Par::get_vel(int n) {
  return vel.at(n);
};

double Par::get_pos(int k) {
  return pos.at(k);
};

vector < double > Par::get_pos_vec() {
  return pos;
};

double Par::get_bpos(int u) {
  return bpos.at(u);
};

vector < double > Par::get_normdir() {

  // generate d points uniformly distributed on a d-dimensional sphere (i.e. random direction)
  // following http://mathworld.wolfram.com/HyperspherePointPicking.html
  std::normal_distribution < double > normdist2(0.0, 1.0);
  vector < double > dir;
  double length = 0.0;
  for (int n = 0; n < dim; n++) {
    double d = normdist2(generator);
    dir.push_back(d);
    length = length + pow(dir.at(n), 2);
  };
  length = sqrt(length);

  for (int n = 0; n < dim; n++) {
    // uniformly positioned vector on a hypersphere
    dir.at(n) = dir.at(n) / length;
  };

  return dir;
};

double Par::get_levy_McCul(double iter, double maxiters) {
  // Generation of Levy symmetric distribution (beta=0) using
  // McCulloch's Algorithm
  // When alpha approaches 2, the distribution becomes Gaussian with mean tau
  // and variance 2 * c ^ 2, and beta has no effect.

  //std::uniform_real_distribution <double> dist1(0.0,1.0)
  std::normal_distribution < double > distnorm(0.0, 1.0);
  // generate variable from exponential distribution
  //double w = -log(dist1(generator));
  // generate variable from uniform distribution
  //double phi = (dist1(generator) - 0.5)*pi;
  //double alpha_init = 1.0;
  //double alpha_fin = 2.0;
  //double alpha = 2.0; // alpha_init + iter * (alpha_fin - alpha_init) / maxiters; //min(2.0,exp(log(2.0)/10000.0*time)); // alpha_init + time * (alpha_fin - alpha_init) / maxiters;
  //double c = 1.0; // 0.01 - (0.01 - 1e-50)*iter / maxiters; // 0.0003*l2_norm(space_range) - 0.0066; // 0.0003*l2_norm(space_range) - 0.0069; // 0.1 - (0.1 - 0.001)*(iter / maxiters);														// scaling parameter
  //cout << 0.0003*l2_norm(space_range) - 0.0066  << endl;
  //double tau = 0.0;														// location parameter
  //double x = c*pow((cos((1 - alpha)*phi)) / w, 1 / alpha - 1)*(sin(alpha*phi)) / pow(cos(phi), 1 / alpha) + tau;
  double x = distnorm(generator); // + tau;
  //cout << 0.0003*l2_norm(space_range) << endl;
  return x;
};

void Par::update_bpos() {
  for (int j = 0; j < dim; j++) {
    bpos.at(j) = pos.at(j);
  };
};

void Par::update_vel(double inertiaf, double CF, vector < double > globpos, double iter) {

  double r1 = dist1(generator);
  double r2 = dist1(generator);

  for (int v = 0; v < dim; v++) {

    // velocity update with perturbations
    vel.at(v) = CF * (inertiaf * vel.at(v) + c1 * r1 * (bpos.at(v) - pos.at(v)) + c2 * r2 * (globpos.at(v) - pos.at(v)));

  };
};

void Par::update_pos() {
  for (int i = 0; i < dim; i++) {
    pos.at(i) = pos.at(i) + vel.at(i);
    std::uniform_real_distribution < double > dist2(mindomain.at(i), maxdomain.at(i));
    if (pos.at(i) > maxdomain.at(i)) {
      pos.at(i) = dist2(generator);
    } else if (pos.at(i) < mindomain.at(i)) {
      pos.at(i) = dist2(generator);
    };
  };
};

void Par::update_pos_levy(vector < double > globpos, double iter, double inertiaf) {
  //std::normal_distribution <double> normdist2(0.0, 1.0);
  double levystep = abs(get_levy_McCul(iter, maxiters));
  vector < double > direction = get_normdir();

  for (int i = 0; i < dim; i++) {
    std::uniform_real_distribution < double > dist2(mindomain.at(i), maxdomain.at(i));
    pos.at(i) = pos.at(i) + levyscale * get_min_dim() * levystep * direction.at(i);
    //pos.at(i) = pos.at(i) + levyscale*(maxdomain.at(i) - mindomain.at(i))*levystep*direction.at(i);

    if (pos.at(i) > maxdomain.at(i)) {
      pos.at(i) = dist2(generator);
    } else if (pos.at(i) < mindomain.at(i)) {
      pos.at(i) = dist2(generator);
    };
  };
};

double Par::eval_fitness(int cycle, int iter, int parid) {
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_core = std::to_string(core);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);

  // count total # func evaluations
  funceval = funceval + 1;

  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo." + str_cycle + "." + str_parID,
     pwd.string() + "/CPU." + str_core + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);
  // check if ffield is LG or not. execute correct reac accordingly.
  if (lg_yn == true) {
    write_ffield_lg(cycle, iter, parid);
    std::ifstream fin5(("CPU." + str_core + "/reac").c_str());
    if (fin5.fail()) {
      cout << "reac executable not found for CPU " + str_core + ". Aborting! \n";
      fin5.close();
      exit(EXIT_FAILURE);
    }
    fin5.close();
    // prepare mandatory files before executing reac in each CPU directory

    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield",
      pwd.string() + "/CPU." + str_core + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
    
    if (fixcharges == true){
      boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/charges",
        pwd.string() + "/CPU." + str_core + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
    };

    // cd to each CPU.x directory
    string old_path = pwd.string();
    boost::filesystem::path p(pwd.string() + "/CPU." + str_core);
    boost::filesystem::current_path(p);

    // execute reac within each CPU.x directory
    //int status = system("./reac > /dev/null 2>&1");
    int status = boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > boost::process::null);
    if (WIFSIGNALED (status)) {
      cout << "Program exited abnormaly on CPU:" << core << "\n";
    };
    // cd back to main directory
    boost::filesystem::path p2(old_path);
    boost::filesystem::current_path(p2);
  } else {
      write_ffield(cycle, iter, parid);
      std::ifstream fin6(("CPU." + str_core + "/reac").c_str());
      if (fin6.fail()) {
        cout << "reac executable not found for CPU " + str_core + ". Aborting! \n";
        fin6.close();
        exit(EXIT_FAILURE);
      }
      fin6.close();

    // prepare mandatory files before executing reac in each CPU directory
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield",
      pwd.string() + "/CPU." + str_core + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
    if (fixcharges == true){
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/charges",
      pwd.string() + "/CPU." + str_core + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
    };

    // cd to each CPU.x directory
    string old_path = pwd.string();
    boost::filesystem::path p(pwd.string() + "/CPU." + str_core);
    boost::filesystem::current_path(p);

    // execute reac within each CPU.x directory
    //system("./reac > run.log");
    //int status = system("./reac > /dev/null 2>&1");
    int status = boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > boost::process::null);
    if (WIFSIGNALED (status)) {
      cout << "reac exited abnormaly on CPU:" << core << "\n";
    };
    //boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > boost::process::null);
    // cd back to main directory
    boost::filesystem::path p2(old_path);
    boost::filesystem::current_path(p2);
  };

  // read fitness value contained in fort.13 file
  string str;
  if ( !boost::filesystem::exists( "CPU." + str_core + "/fort.13" ) )
  {
    fitness = numeric_limits < double > ::infinity();
  } else {
    boost::filesystem::ifstream file13("CPU." + str_core + "/fort.13");
    stringstream tempstr;
    getline(file13, str);
    // insert str into stringstream tempstr
    tempstr << str;
    // get rid from extra whitespace in stringstream
    tempstr >> std::ws;
    // insert back to str
    tempstr >> str;
    // check if fitness is numeric or ******
    if (str.at(0) == '*') {
      fitness = numeric_limits < double > ::infinity();
    } else {
      // convert to double
      fitness = stod(str);
      file13.close();
    };
    boost::filesystem::remove( "CPU." + str_core + "/fort.13" );
  };


  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/dipole.out" ,
    //pwd.string()+"/CPU."+str_core+"/current_dipole.out",boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
    //pwd.string()+"/CPU."+str_core+"/molgeo.out."+str_cycle+"."+str_iter, 
    //boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.74" ,
    //pwd.string()+"/CPU."+str_core+"/thermo.out."+str_cycle+"."+str_iter+"."+str_parID, 
    //boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.99" ,
    //pwd.string()+"/CPU."+str_core+"/results.out."+str_cycle+"."+str_iter+"."+str_parID, 
    //boost::filesystem::copy_option::overwrite_if_exists);

  // Note: do not update the geometry file during iterations. Each member should use one geo file throughout
  // the training. Assuming we start with DFT_optimized (or sensible structures), and that we use some small
  // number of structural minimizations (3-10), the sensible structures won't change much, so no need to
  // provide previous structure as initial structure for the next round of iteration - since momentarily bad
  // parameters could completely destroy the structure and make the successive minimization work on a crazy
  // structure.
  // boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
  //   pwd.string()+"/CPU."+str_core+"/geo." + str_cycle + "." + str_parID, 
  //     boost::filesystem::copy_option::overwrite_if_exists);

  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.73" ,
    //pwd.string()+"/CPU."+str_core+"/partialE.out."+str_cycle+"."+str_iter+"."+str_parID, 
    //boost::filesystem::copy_option::overwrite_if_exists);
  return fitness;
#endif

#ifndef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);

  // count total # func evaluations
  funceval = funceval + 1;
  boost::filesystem::copy_file(pwd.string() + "/geo." + str_cycle + "." + str_parID,
    pwd.string() + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);
  // check if ffield is LG or not. execute correct reac accordingly.
  if (lg_yn == true) {
    write_ffield_lg(cycle, iter, parid);
    std::ifstream fin5("reac");
    if (fin5.fail()) {
      cout << "reac executable not found. Aborting! \n";
      fin5.close();
      exit(EXIT_FAILURE);
    }
    fin5.close();

    // prepare mandatory files before executing reac

    boost::filesystem::copy_file(pwd.string() + "/ffield",
      pwd.string() + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);

    if (fixcharges == true){
      boost::filesystem::copy_file(pwd.string() + "/charges",
        pwd.string() + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
    };

    // execute reac
    //system("./reac > /dev/null 2>&1"); 
    int status = boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > boost::process::null);
    if (WIFSIGNALED (status)) {
      cout << "reac exited abnormaly." << endl;
    };

  } else {
    write_ffield(cycle, iter, parid);
    std::ifstream fin6("reac");
    if (fin6.fail()) {
      cout << "reac executable not found. Aborting! \n";
      fin6.close();
      exit(EXIT_FAILURE);
    }
    fin6.close();

    // prepare mandatory files before executing reac

    boost::filesystem::copy_file(pwd.string() + "/geo",
      pwd.string() + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

    boost::filesystem::copy_file(pwd.string() + "/ffield",
      pwd.string() + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);

    if (fixcharges == true){
    boost::filesystem::copy_file(pwd.string() + "/charges",
      pwd.string() + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
    };

    // execute reac within each CPU.x directory
    //system("./reac > /dev/null 2>&1");
    int status = boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > boost::process::null); 
    if (WIFSIGNALED (status)) {
      cout << "reac exited abnormaly." << endl;
    };

  };

  // read fitness value contained in fort.13 file
  string str;
  if ( !boost::filesystem::exists("fort.13") )
  {
    fitness = numeric_limits < double > ::infinity();
  } else {
    boost::filesystem::ifstream file13("fort.13");
    stringstream tempstr;
    getline(file13, str);
    // insert str into stringstream tempstr
    tempstr << str;
    // get rid from extra whitespace in stringstream
    tempstr >> std::ws;
    // insert back to str
    tempstr >> str;
    // check if fitness is numeric or ******
    if (str.at(0) == '*') {
      fitness = numeric_limits < double > ::infinity();
    } else {
      // convert to double
      fitness = stod(str);
    };
    file13.close();
    boost::filesystem::remove("fort.13");
  };

  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/dipole.out" ,
    //pwd.string()+"/CPU."+str_core+"/current_dipole.out",boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/fort.90" ,
     //pwd.string() + "/molgeo.out." + str_cycle + "." + str_iter,
     //boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/fort.74" ,
     //pwd.string() + "/thermo.out." + str_cycle+"." + str_iter + "." + str_parID,
     //boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/fort.99" ,
    //pwd.string() + "/results.out." + str_cycle+"." + str_iter + "." + str_parID,
    //boost::filesystem::copy_option::overwrite_if_exists);

    // Note: do not update the geometry file during iterations. Each member should use one geo file throughout
    // the training. Assuming we start with DFT_optimized (or sensible structures), and that we use some small
    // number of structural minimizations (3-10), the sensible structures won't change much, so no need to
    // provide previous structure as initial structure for the next round of iteration - since momentarily bad
    // parameters could completely destroy the structure and make the successive minimization work on a crazy
    // structure.
    // boost::filesystem::copy_file(pwd.string() + "/fort.90" ,
    //   pwd.string() + "/geo."+str_cycle+"." + str_parID ,
    //      boost::filesystem::copy_option::overwrite_if_exists);


   //boost::filesystem::copy_file(pwd.string() + "/fort.73" ,
     //pwd.string() + "/partialE.out." + str_cycle + "." + str_iter + "." + str_parID,
     //boost::filesystem::copy_option::overwrite_if_exists);
  return fitness;
#endif
};

double Par::get_bfit() {
  return bfitness;
};

double Par::get_fitness() {
  return fitness;
};

void Par::set_fitness(double fit) {
  fitness = fit;
};

void Par::set_bfit(double bfit) {
  bfitness = bfit;
};

void Par::set_pos(vector < double > pos_of_best_particle) {
  pos = pos_of_best_particle;
};

void Par::set_vel(vector < double > vel_of_best_particle) {
  vel = vel_of_best_particle;
};

// ---------- Definitions of Swarm class member functions ---------- //
//
//------------------------------------------------------------------ //

Swarm::Swarm() {

};

void Swarm::read_icharg_control() {
  // Read icharg value from reax control file
  string str_core = std::to_string(core);
  std::ifstream control_file("control");
  if (control_file.fail()) {
    cout << "Unable to open control file on CPU " << core << ". \n";
    control_file.close();
    exit(EXIT_FAILURE);
  } else{
    std::string line;
    int numlines = 0;
    // for each line
    while (std::getline(control_file, line)) {
      boost::trim(line);
      if ( line[0] == '#' || line == ""){
        continue;
      }else{
        numlines++;
        // locate icharg line
        if (numlines == 12) {
          if (line[0] == '5') {
            fixcharges = true;
            break;
          };
        };
      };
    };
  };
};


void Swarm::get_userinp(){
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  // prepare dirs for each CPU process
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core);

  if (core == 0){
    vector <string> tempinput;
    std::ifstream fin("inp_flocky.in");
    if (fin.fail()) {
      cout << "Unable to open 'inp_flocky.in' file \n";
      fin.close();
      exit(EXIT_FAILURE);
    } else {
      std::string line;
      // for each line
      tempinput.clear();
      while (std::getline(fin, line)) {
        boost::trim(line);
        if ( line[0] == '#' || line == ""){
          continue;
        }else{
          // create a new line to split
          vector <string> values;
          //std::istringstream lineStream(line);
          using boost::is_any_of;
          boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
          tempinput.push_back(values.at(0));
        };
      };
      fin.close();
    };

    istringstream(tempinput.at(0)) >> lg_yn;
    istringstream(tempinput.at(1)) >> contff;
    istringstream(tempinput.at(2)) >> perc_yn;
    istringstream(tempinput.at(3)) >> perc;
    istringstream(tempinput.at(4)) >> ofit;
    istringstream(tempinput.at(5)) >> uq;
    istringstream(tempinput.at(6)) >> NumP;
    if (NumP < numcores) {
      cout << "Error: Number of swarm members should be bigger than number of allocated processors." << endl;
      exit(EXIT_FAILURE);
    } else {
      NumP = int(floor(NumP / numcores));
    };
    istringstream(tempinput.at(7)) >> c1;
    istringstream(tempinput.at(8)) >> c2;
    istringstream(tempinput.at(9)) >> inertiamax;
    istringstream(tempinput.at(10)) >> inertiamin; 
    istringstream(tempinput.at(11)) >> faili;
    istringstream(tempinput.at(12)) >> levyscale;
    istringstream(tempinput.at(13)) >> freq;
    istringstream(tempinput.at(14)) >> maxiters;
    istringstream(tempinput.at(15)) >> maxcycles;
  };  
  // check if reaxff was set to run with fixed charges and require charges file
  read_icharg_control();
  boost::filesystem::ifstream charge_file("charges");
  //fixcharges = true;
  if (fixcharges == true && charge_file.fail()) {
    cout << "'control' uses icharg=5 (fixed charges) but no 'charges' file was found!" << endl;
    charge_file.close();
    exit(EXIT_FAILURE);
  };
  if (fixcharges == true) {
    boost::filesystem::copy_file(pwd.string() + "/charges", pwd.string() + "/CPU." + str_core + "/charges",
      boost::filesystem::copy_option::overwrite_if_exists);
  // fixcharges = true;
  };
  charge_file.close();

  if (lg_yn == true){
  boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/CPU." + str_core + "/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
  }else{
  boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/CPU." + str_core + "/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file(pwd.string() + "/ffield", pwd.string() + "/CPU." + str_core + "/ffield",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/control", pwd.string() + "/CPU." + str_core + "/control",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/geo", pwd.string() + "/CPU." + str_core + "/geo",
    boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/trainset.in", pwd.string() + "/CPU." + str_core + "/trainset.in",
    boost::filesystem::copy_option::overwrite_if_exists);

  // create fort.20 file needed by reac and reac.lg
  boost::filesystem::ofstream iopt_file(pwd.string() + "/CPU." + str_core + "/fort.20");
  iopt_file << "0";
  iopt_file.close();

  // create fort.35 file needed by reac and reac.lg
  ofstream fort35_file(pwd.string() + "/CPU." + str_core + "/fort.35");
  fort35_file << "23434.1" << endl;
  fort35_file.close();

  // broadcast between all processes from process 0
  // MPI_Bcast must be visible to all processes!
  MPI_Bcast( & lg_yn, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
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
  //MPI_Bcast( & fixcharges, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & ofit, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & uq, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif

#ifndef WITH_MPI
    vector <string> tempinput;
    std::ifstream fin("inp_flocky.in");
    if (fin.fail()) {
      cout << "Unable to open 'inp_flocky.in' file \n";
      fin.close();
      exit(EXIT_FAILURE);
    } else {
      std::string line;
      // for each line
      tempinput.clear();
      while (std::getline(fin, line)) {
        boost::trim(line);
        if ( line[0] == '#' || line == ""){
          continue;
        }else{
          // create a new line to split
          vector <string> values;
          //std::istringstream lineStream(line);
          using boost::is_any_of;
          boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
          tempinput.push_back(values.at(0));
        };
      };
      fin.close();
    };

    istringstream(tempinput.at(0)) >> lg_yn;
    istringstream(tempinput.at(1)) >> contff;
    istringstream(tempinput.at(2)) >> perc_yn;
    istringstream(tempinput.at(3)) >> perc;
    istringstream(tempinput.at(4)) >> ofit;
    istringstream(tempinput.at(5)) >> uq;
    istringstream(tempinput.at(6)) >> NumP;
    istringstream(tempinput.at(7)) >> c1;
    istringstream(tempinput.at(8)) >> c2;
    istringstream(tempinput.at(9)) >> inertiamax;
    istringstream(tempinput.at(10)) >> inertiamin;
    istringstream(tempinput.at(11)) >> faili;
    istringstream(tempinput.at(12)) >> levyscale;
    istringstream(tempinput.at(13)) >> freq;
    istringstream(tempinput.at(14)) >> maxiters;
    istringstream(tempinput.at(15)) >> maxcycles;

  // check if reaxff was set to run with fixed charges and require charges file
  read_icharg_control();
  boost::filesystem::ifstream charge_file("charges");
  if (charge_file.fail()) {
    cout << "'control' uses icharg=5 (fixed charges) but no 'charges' file was found!" << endl;
    charge_file.close();
    exit(EXIT_FAILURE);
  };

  // create fort.20 file needed by reac and reac.lg
  boost::filesystem::ofstream iopt_file("fort.20");
  iopt_file << "0";
  iopt_file.close();

  // create fort.35 file needed by reac and reac.lg
  ofstream fort35_file("fort.35");
  fort35_file << "23434.1" << endl;
  fort35_file.close();

#endif
};

vector <double> Swarm::get_com(Swarm newSwarm) {
 vector <double> swarmcom(dim, 0.0);
 double sumpos = 0.0;
 double totsumpos = 0.0;
#ifdef WITH_MPI
 //NumP = numcores*NumP;
 for (int m = 0; m < dim; m++){
   for (int p = 0; p < NumP; p++){
     sumpos = sumpos + 1.0/(float(NumP)*numcores)*newSwarm.GetPar(p).get_pos(m);
   };
   MPI_Allreduce( & sumpos,  & totsumpos, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   swarmcom.at(m) = totsumpos;
   sumpos = 0.0;
   totsumpos = 0.0;
#endif
#ifndef WITH_MPI
 for (int m = 0; m < dim; m++){
   for (int p = 0; p < NumP; p++){
     sumpos = sumpos + 1.0/(float(NumP)*newSwarm.GetPar(p).get_pos(m));
   };
   swarmcom.at(m) = sumpos;
   sumpos = 0.0;
#endif
 };
/* if (core == 0) {
   cout << "com(0) on CPU: " << core << " is: " << swarmcom.at(0) << endl;
   cout << "com(1) on CPU: " << core << " is: " << swarmcom.at(1) << endl;
 };
*/
 return swarmcom;
};

double Swarm::get_disp(Swarm newSwarm) {
  double swarm_dispersion = 0.0;
  // define deviation vector between a swarm member and the center-of-mass
  vector <double> deviation_vec(dim, 0.0);
  double deviation = 0.0;
  double totdeviation = 0.0;
  // calculate deviation vector
  for (int m = 0; m < dim; m++){
   for (int p = 0; p < NumP; p++){
     deviation = deviation + pow(newSwarm.GetPar(p).get_pos(m) - newSwarm.get_com(newSwarm).at(m),2);
   };
#ifdef WITH_MPI
  MPI_Allreduce(& deviation, & totdeviation, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  deviation_vec.at(m) = sqrt(totdeviation);
  deviation = 0.0;
  totdeviation = 0.0;
#endif
#ifndef WITH_MPI
  deviation_vec.at(m) = sqrt(deviation);
  deviation = 0.0;
#endif
  };

  // To calculate the median separation between swarm members and COM
  // we first sort the array
  int n = deviation_vec.size();
  sort(deviation_vec.begin(), deviation_vec.end()); 
  
  // calculate median
  // check for even case
  if (n % 2 != 0) {
    swarm_dispersion = deviation_vec.at(n/2.0);
  }else{
    swarm_dispersion = (deviation_vec.at((n-1)/2.0) + deviation_vec.at(n/2))/2.0;
  };

  return swarm_dispersion;
};

Par & Swarm::GetPar(int ParID) {
  return AllParticles[ParID];
};

void Swarm::AddPar(Par & newPar) {
  AllParticles.push_back(newPar);
};

int Swarm::get_worse(Swarm newSwarm) {
  double fit = -1.0 * numeric_limits < double > ::infinity();
  int particle_id;
  for (int p = 0; p < NumP; p++) {
    if (newSwarm.GetPar(p).get_fitness() > fit) {
      fit = newSwarm.GetPar(p).get_fitness();
      particle_id = p;
    }
  };
  return particle_id;
};

int Swarm::get_best(Swarm newSwarm) {
  double fit = numeric_limits < double > ::infinity();
  int particle_id = 0;
  for (int p = 0; p < NumP; p++) {
    if (newSwarm.GetPar(p).get_fitness() < fit) {
      fit = newSwarm.GetPar(p).get_fitness();
      particle_id = p;
    }
  };
  return particle_id;
};


void Swarm::Populate(Swarm & newSwarm, int cycle) {
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  string itercount = std::to_string(cycle);
  if (core == 0) {
    boost::filesystem::ofstream log("log.flocky", ofstream::app);
    log << "\n";
    log << "Swarm generation started. Please wait." << endl;
    log.close();
  };
  funceval = 0; // clear counter for cycles
  // ---------------------------------------------- //
  //     POPULATE: MAIN LOOP OVER SWARM MEMBERS
  // ---------------------------------------------- //

  for (int p = 0; p < NumP; p++) {
    string parID = std::to_string(p);
    string str_cycle = std::to_string(cycle);
    // cp geo to geo.parID so each particle works with its own geo file
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo",
      pwd.string() + "/CPU." + str_core + "/geo." + str_cycle + "." + parID,
        boost::filesystem::copy_option::overwrite_if_exists);

    Par NewPar;
    newSwarm.AddPar(NewPar);

    // initial gbfit is INF for all processes and all swarm members
    gbfit = numeric_limits < double > ::infinity();
    // initial gbpos is 0.0 for all processes and all swarm members
    gbpos.clear();
    for (int m=0; m < dim; m++){
      gbpos.push_back(0.0);
    };

    if (core == 0) {
      // If contff == y, then take force field's current values for the position of particle 0 (others are random)
      if (contff == true) {
        vector < double > ffpos;
        ffpos.clear();
        for (int m = 0; m < dim; m++) {
          ffpos.push_back(stod(newSwarm.GetPar(0).ffieldmat.at(newSwarm.GetPar(0).ffline.at(m)).at(newSwarm.GetPar(0).ffcol.at(m))));
        };
        newSwarm.GetPar(0).set_pos(ffpos);
      };
      contff = false;
    };

    // evaluate fitness and set bfit = curfit
    curfit = newSwarm.GetPar(p).eval_fitness(cycle, 0, p);
    newSwarm.GetPar(p).set_fitness(curfit);
    newSwarm.GetPar(p).set_bfit(curfit);
    if (curfit < gbfit) {
      gbfit = curfit;
      gbpos.clear();
      gbpos = newSwarm.GetPar(p).get_pos_vec();
      write_ffield_gbest(core, cycle, 0, p);
    };

    // cleaning ffield.tmp.* files after write_ffield_gbest already
    // copied the correct ffield.tmp.* file as the ffield.gbest.*.0.*
    boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp." + str_cycle+".0." + "." + parID);
  }; // done loop on members
  // pair struct to hold the global best fitness across processes and its core rank
  struct {
    double tmp_fit;
    int tmp_cpu;
  } min_vals_in[1], min_vals_out[1];
  // store current fit on each process
  min_vals_in[0].tmp_fit = gbfit;
  // store core id of that current process
  min_vals_in[0].tmp_cpu = core;
  // get global best fitness *across processes* and corresponding core rank and store them in min_vals_out
  MPI_Allreduce( & min_vals_in, & min_vals_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  // global best fitness across all processes
  gbfit = min_vals_out[0].tmp_fit;
  // core rank the above fitness came from
  cpuid_gbfit = min_vals_out[0].tmp_cpu;
  // broadcast contents of gbpos vector from rank cpuid_gbfit
  MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, MPI_COMM_WORLD);

  if (core == 0) {
    // clean up any old files belonging to previous job
    if ( boost::filesystem::exists( "opti_log.out." + std::to_string(cycle)) ){
      boost::filesystem::remove( "opti_log.out." + std::to_string(cycle) );
    };
    if ( boost::filesystem::exists( "disp_log.out." + std::to_string(cycle)) ){
      boost::filesystem::remove( "disp_log.out." + std::to_string(cycle) );
    };

    newSwarm.printopt(newSwarm, 0, cycle, 1);
    boost::filesystem::ofstream log("log.flocky", ofstream::app);
    log << "\nSwarm generation completed." << endl;
    log << "Initial global best fit: " << gbfit << endl;
    log << "flocky optimization started!" << endl;
    log.close();
  };

  //newSwarm.printdisp(newSwarm, 0, cycle, 1);
  newSwarm.printpos(newSwarm, 0, cycle, 1);
  if (uq == true) {
    newSwarm.printUQFF(newSwarm, 0, cycle, 1);
    newSwarm.printUQQoI(newSwarm, 0, cycle, 1);
  };
#endif

#ifndef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string itercount = std::to_string(cycle);
  boost::filesystem::ofstream log("log.flocky", ofstream::app);
  log << "\n";
  log << "Swarm generation started. Please wait." << endl;
  log.close();

  // ---------------------------------------------- //
  //     POPULATE: MAIN LOOP OVER SWARM MEMBERS
  // ---------------------------------------------- //

  for (int p = 0; p < NumP; p++) {
    string parID = std::to_string(p);
    string str_cycle = std::to_string(cycle);
    
    // cp geo to geo.parID so each particle works with its own geo file
    boost::filesystem::copy_file(pwd.string() + "/geo",
      pwd.string() + "/geo." + str_cycle + "." + parID,
        boost::filesystem::copy_option::overwrite_if_exists);
    
    Par NewPar;
    newSwarm.AddPar(NewPar);

    // initial gbfit is INF for all processes and all swarm members
    gbfit = numeric_limits < double > ::infinity();
    // initial gbpos is 0.0 for all processes and all swarm members
    gbpos.clear();
    for (int m=0; m < dim; m++){
      gbpos.push_back(0.0);
    };

    // If contff == y, then take force field's current values for the position of particle 0 (others are random)
    if (contff == true) {
      vector < double > ffpos;
      ffpos.clear();
      for (int m = 0; m < dim; m++) {
        ffpos.push_back(stod(newSwarm.GetPar(0).ffieldmat.at(newSwarm.GetPar(0).ffline.at(m)).at(newSwarm.GetPar(0).ffcol.at(m))));
      };
      newSwarm.GetPar(0).set_pos(ffpos);
    };
    contff = false;

    // evaluate fitness and set bfit = curfit
    curfit = newSwarm.GetPar(p).eval_fitness(cycle, 0, p);
    newSwarm.GetPar(p).set_fitness(curfit);
    newSwarm.GetPar(p).set_bfit(curfit);

    if (curfit < gbfit) {
      gbfit = curfit;
      gbpos.clear();
      gbpos = newSwarm.GetPar(p).get_pos_vec();
      write_ffield_gbest(0, cycle, 0, p);
    };
    // cleaning ffield.tmp.* files after write_ffield_gbest already
    // copied the correct ffield.tmp.* file as the ffield.gbest.*.0.*
    boost::filesystem::remove("ffield.tmp." + str_cycle + ".0." + "." + parID);
  }; // done loop over members

  //initial_disp = newSwarm.get_disp(newSwarm);

  // clean up any old files belonging to previous job
  if ( boost::filesystem::exists( "opti_log.out." + std::to_string(cycle)) ){
    boost::filesystem::remove( "opti_log.out." + std::to_string(cycle) );
  };
  if ( boost::filesystem::exists( "disp_log.out." + std::to_string(cycle)) ){
    boost::filesystem::remove( "disp_log.out." + std::to_string(cycle) );
  };

  newSwarm.printopt(newSwarm, 0, cycle, 1);
  //boost::filesystem::ofstream log("log.flocky", ofstream::app);
  log << "\nSwarm generation completed." << endl;
  log << "Initial global best fit: " << gbfit << endl;
  log << "flocky optimization started!" << endl;
  log.close();

  //newSwarm.printdisp(newSwarm, 0, cycle, 1);
  newSwarm.printpos(newSwarm, 0, cycle, 1);
  if (uq == true){
    newSwarm.printUQFF(newSwarm, 0, cycle, 1);
    newSwarm.printUQQoI(newSwarm, 0, cycle, 1);
  };
#endif

};



void Swarm::Propagate(Swarm & newSwarm, int cycle) {
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);

  // ---------------------------------------------------------------- //
  //         PROPAGATE: main loop over iterations
  // ---------------------------------------------------------------- //
  for (iter = 1; iter < maxiters; iter++) {
    inertiafac = inertiamax - iter * (inertiamax - inertiamin) / maxiters;
    //
    // main loop over swarm members
    //
    for (int p = 0; p < NumP; p++) {
      // Update velocities and positions
      newSwarm.GetPar(p).update_vel(inertiafac, confac, gbpos, iter);
      if (newSwarm.GetPar(p).fails > faili) {
        newSwarm.GetPar(p).update_pos_levy(gbpos, iter, inertiafac);
        newSwarm.GetPar(p).fails = 0;
      } else {
        newSwarm.GetPar(p).update_pos();
      };
      curfit = newSwarm.GetPar(p).eval_fitness(cycle, iter, p);
      newSwarm.GetPar(p).set_fitness(curfit);
      // Update personal best positions and fitness
      if (newSwarm.GetPar(p).get_fitness() < newSwarm.GetPar(p).get_bfit()) {
        newSwarm.GetPar(p).update_bpos();
        newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());
        newSwarm.GetPar(p).fails = 0;
        if (newSwarm.GetPar(p).get_bfit() < gbfit) {
          gbfitfound = true;
          parid_gbfit = p;
          gbfit = newSwarm.GetPar(p).get_bfit();
          gbpos.clear();
          gbpos = newSwarm.GetPar(p).get_pos_vec();
          write_ffield_gbest(core, cycle, iter, p);
        };
      } else {
        gbfitfound = false;
        newSwarm.GetPar(p).fails = newSwarm.GetPar(p).fails + 1;
      };
      // cleaning ffield.tmp.* files after write_ffield_gbest already
      // copied the correct ffield.tmp.* file as the ffield.gbest.*
      string str_cycle = std::to_string(cycle);
      string str_iter = std::to_string(iter);
      string str_parID = std::to_string(p);
      //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID);
    }; // done loop over swarm members
    //newSwarm.get_com(newSwarm);

    // pair struct to hold the global best fitness across processes and its core rank
    struct {
      double tmp_fit;
      int tmp_cpu;
    } min_vals_in[1], min_vals_out[1];

    // store current fit on each process
    min_vals_in[0].tmp_fit = gbfit;
    // store core id of that current process
    min_vals_in[0].tmp_cpu = core;
    // get global best fitness *across processes* and corresponding core rank and store them in min_vals_out
    MPI_Allreduce( & min_vals_in, & min_vals_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    // pair struct to hold the global best fitness across processes and its parID
    // The parID is required in detection of overfitting to cp the correct ffield.gbest file
    struct {
      double tmp_fit2;
      int tmp_parid;
    } min_vals_in2[1], min_vals_out2[1];

    // store current fit on each process
    min_vals_in2[0].tmp_fit2 = gbfit;
    // store par id of that current process
    min_vals_in2[0].tmp_parid = parid_gbfit;
    // get global best fitness *across processes* and corresponding parID and store them in min_vals_out
    MPI_Allreduce( & min_vals_in2, & min_vals_out2, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    // parid of the gbfit
    parid_gbfit = min_vals_out2[0].tmp_parid;


    // global best fitness across all processes
    gbfit = min_vals_out[0].tmp_fit;
    // core rank the above fitness came from
    cpuid_gbfit = min_vals_out[0].tmp_cpu;
    // broadcast contents of gbpos vector from rank cpuid_gbfit
    MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, MPI_COMM_WORLD);
    if (ofit == true){
      if (gbfitfound == true) {
        // detect overfitting by evaluating fitness on validation set
        if (core == cpuid_gbfit){
          boost::filesystem::ofstream log("log.flocky", ofstream::app);
          newSwarm.detovfit(newSwarm, cpuid_gbfit, cycle, iter, parid_gbfit);
          log.close();
        };
        firstovfit = false;
      };
      MPI_Bcast( &firstovfit, 1, MPI_C_BOOL, cpuid_gbfit, MPI_COMM_WORLD);
    };

    if (core == 0) {
      newSwarm.printopt(newSwarm, iter, cycle, freq);
    };
    newSwarm.printpos(newSwarm, iter, cycle, freq);
    if (uq == true){
      newSwarm.printUQFF(newSwarm, iter, cycle, freq);
      newSwarm.printUQQoI(newSwarm, iter, cycle, freq);
    };
    //newSwarm.printdisp(newSwarm, iter, cycle, freq);
  }; // done loop over iterations
  MPI_Allreduce(& funceval, & funceval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (core == 0) {
    boost::filesystem::ofstream log("log.flocky", ofstream::app);
    log << "\nTraining completed successfuly!\n";
    log << "Total ReaxFF calls: " << funceval << endl;
    log << "Global best ReaxFF fit: " << newSwarm.get_gbfit() << endl;
    log << "Global best ReaxFF parameters:" << endl;
    log << "[ ";
    for (int m = 0; m < dim; m++) {
      log << newSwarm.get_gbpos().at(m) << " ";
    };
    log << "]" << endl;
    log.close();
  };
#endif

#ifndef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);

  // ---------------------------------------------------------------- //
  //         PROPAGATE: main loop over iterations
  // ---------------------------------------------------------------- //
  for (iter = 1; iter < maxiters; iter++) {
    inertiafac = inertiamax - iter * (inertiamax - inertiamin) / maxiters;
    //
    // main loop over swarm members
    //
    for (int p = 0; p < NumP; p++) {
      // Update velocities and positions
      newSwarm.GetPar(p).update_vel(inertiafac, confac, gbpos, iter);

      if (newSwarm.GetPar(p).fails > faili) {
        newSwarm.GetPar(p).update_pos_levy(gbpos, iter, inertiafac);
        newSwarm.GetPar(p).fails = 0;
      } else {
        newSwarm.GetPar(p).update_pos();
      };
      curfit = newSwarm.GetPar(p).eval_fitness(cycle, iter, p);
      newSwarm.GetPar(p).set_fitness(curfit);

      // Update personal best positions and fitness
      if (newSwarm.GetPar(p).get_fitness() < newSwarm.GetPar(p).get_bfit()) {
        newSwarm.GetPar(p).update_bpos();
        newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());
        newSwarm.GetPar(p).fails = 0;
        if (newSwarm.GetPar(p).get_bfit() < gbfit) {
          gbfitfound = true;
          parid_gbfit = p;
          gbfit = newSwarm.GetPar(p).get_bfit();
          gbpos.clear();
          gbpos = newSwarm.GetPar(p).get_pos_vec();
	  write_ffield_gbest(0, cycle, iter, p);
        };
      } else {
        gbfitfound = false;
        newSwarm.GetPar(p).fails = newSwarm.GetPar(p).fails + 1;
      };
      // cleaning ffield.tmp.* files after write_ffield_gbest already
      // copied the correct ffield.tmp.* file as the ffield.gbest.*
      string str_cycle = std::to_string(cycle);
      string str_iter = std::to_string(iter);
      string str_parID = std::to_string(p);
      //boost::filesystem::remove(pwd.string() + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID);
    }; // done loop over swarm members

    //newSwarm.get_com(newSwarm);
   
    if (ofit == true) {
      if (gbfitfound == true){
        // detect overfitting by evaluating fitness on validation set
        newSwarm.detovfit(newSwarm, 0, cycle, iter, parid_gbfit);
        firstovfit = false;
      };
    };

    newSwarm.printopt(newSwarm, iter, cycle, freq);
    newSwarm.printpos(newSwarm, iter, cycle, freq);
    if (uq == true){
      newSwarm.printUQFF(newSwarm, iter, cycle, freq);
      newSwarm.printUQQoI(newSwarm, iter, cycle, freq);
    };
    //newSwarm.printdisp(newSwarm, iter, cycle, freq);
  }; // done loop over iterations

  if (core == 0) {
    boost::filesystem::ofstream log("log.flocky", ofstream::app);
    log << "\nTraining completed successfuly!\n";
    log << "Total ReaxFF calls: " << funceval << endl;
    log << "Global best ReaxFF fit: " << newSwarm.get_gbfit() << endl;
    log << "Global best ReaxFF parameters:" << endl;
    log << "[ ";
    for (int m = 0; m < dim; m++) {
      log << newSwarm.get_gbpos().at(m) << " ";
    };
    log << "]" << endl;
    log.close();
  //cout << "Total skipped particles: " << maxiters*NumP - funceval << endl;
  };
#endif
};

void Swarm::write_ffield_gbest(int core, int cycle, int iter, int par) {
  // cp current ffield to be the global best and current analysis files to global best analysis files
#ifndef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(par);

  boost::filesystem::ifstream test_ffieldtmp(pwd.string() + "ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID);
  if (!test_ffieldtmp.fail()) {
    boost::filesystem::copy_file(pwd.string() + "/ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID,
    "ffield.gbest." + str_cycle + "." + str_iter + "." + str_parID, boost::filesystem::copy_option::overwrite_if_exists);
  };
  test_ffieldtmp.close();

  boost::filesystem::ifstream test_fort99(pwd.string() + "/fort.99");
  if (!test_fort99.fail()) {
    boost::filesystem::copy_file(pwd.string() + "/fort.99",
      "results.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  };
  test_fort99.close();

  //boost::filesystem::copy_file(pwd.string() + "/fort.74",
    //"thermo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/fort.90",
    //"molgeo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/fort.73",
    //"partialE.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/current_dipole.out" ,
  //  "dipole.out"+itercount,boost::filesystem::copy_option::overwrite_if_exists);

  // clean temp files
  //boost::filesystem::remove(pwd.string() + "/ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID);
#endif
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_cycle = std::to_string(cycle);
  string str_core = std::to_string(core);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(par);

  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    "ffield.gbest." + str_cycle + "." + str_iter+"."+str_parID, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.99",
    "results.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.74",
    //"thermo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.90",
    //"molgeo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.73",
    //"partialE.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/current_dipole.out" ,
  //  "dipole.out"+itercount,boost::filesystem::copy_option::overwrite_if_exists);

  // clean temp files
  //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID);
#endif
};

void Swarm::detovfit(Swarm &newSwarm, int cpuid_gbfit, int cycle, int iter, int parid_gbfit) {
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(cpuid_gbfit);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(parid_gbfit);
  // prepare mandatory files for reaxff in 'testovfit' directory
  if (firstovfit == true) {
      boost::filesystem::create_directory("testovfit");
      boost::filesystem::copy_file(pwd.string() + "/ffield.initial." + str_cycle,
        pwd.string() + "/testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    boost::filesystem::copy_file(pwd.string() + "/ffield.gbest." + str_cycle + "." + str_iter + "." + str_parID,
      pwd.string() + "/testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file(pwd.string() + "/geo.val",
    pwd.string() + "/testovfit/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

  if (fixcharges == true){
    boost::filesystem::copy_file(pwd.string() + "/charges_val",
      pwd.string() + "/testovfit/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file(pwd.string() + "/control", pwd.string() + "/testovfit/control",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::copy_file(pwd.string() + "/valset.in", pwd.string() + "/testovfit/trainset.in",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::ofstream iopt_file(pwd.string() + "/testovfit/fort.20");
  iopt_file << "0";
  iopt_file.close();
  ofstream outfile35(pwd.string() + "/testovfit/fort.35");
  outfile35 << "23434.1" << endl;
  outfile35.close();

  // cd to testovfit directory
  string old_path = pwd.string();
  boost::filesystem::path p(pwd.string() + "/testovfit");
  boost::filesystem::current_path(p);

  if (lg_yn == true) {
    // execute reac or reac
    boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/testovfit/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
    system("./reac > /dev/null 2>&1");
  } else {
    boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/testovfit/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
    system("./reac > /dev/null 2>&1");
  };

  // cd back to main directory
  boost::filesystem::path p2(old_path);
  boost::filesystem::current_path(p2);

  // read fitness value contained in newly generated fort.13 file
  string str;
  boost::filesystem::ifstream myfile("testovfit/fort.13");
  stringstream tempstr;
  tempstr.str("");
  getline(myfile, str);
  // insert str into stringstream tempstr
  tempstr << str;
  // get rid from extra whitespace in stringstream
  tempstr >> std::ws;
  // insert back to str
  tempstr >> str;
  // check if fitness is numeric or ******
  if (str.at(0) == '*') {
    currovfitness = numeric_limits < double > ::infinity();
  } else {
    // convert to double
    currovfitness = stod(str);
  };
  myfile.close();
  boost::filesystem::remove( "testovfit/fort.13" );
  if (firstovfit == true) {
  //  cout << "firstovfit is true" << endl;
  ovfitness = currovfitness;
  };
  //cout << "firstovfit on CPU: " << core << " is: " << firstovfit << endl;

    //cout << "current over-fitting fitness: " << currovfitness << endl;
    //cout << "previous over-fitting fitness: " << ovfitness << endl;

  //if (currovfitness > ovfitness) {
    //cout << "\nOverfitting detected! current fitness: " << currovfitness << "> previous fitness: " << ovfitness << endl;
    //cout << "Early-stopping ReaxFF optimization.\n";
    //cout << "Training completed successfuly!\n";
    //cout << "Total ReaxFF calls: " << funceval << endl;
    //cout << "Global best ReaxFF fit: " << newSwarm.get_gbfit() << endl;
    //cout << "Global best ReaxFF parameters:\n";
    //cout << "[ ";
    //for (int m = 0; m < dim; m++) {
    //  cout << newSwarm.get_gbpos().at(m) << " ";
    //};
    //cout << "]" << endl;
    //ierr = MPI_Finalize();
    //exit(EXIT_SUCCESS);
    ofstream outfileovfit("overfit.out." + std::to_string(cycle), ofstream::app);
    stringstream ss3;
    ss3 << boost::format("%5i %25.4f") %iter %currovfitness;
    outfileovfit << ss3.str();
    ss3.str("");
    outfileovfit << endl;
    outfileovfit.close();
  //}else{
   //ovfitness = currovfitness;
  //};
#endif
#ifndef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(parid_gbfit);
  // prepare mandatory files for reaxff in 'testovfit' directory
  if (firstovfit == true) {
      boost::filesystem::create_directory("testovfit");
      boost::filesystem::copy_file(pwd.string() + "/ffield.initial." + str_cycle,
        pwd.string() + "/testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    boost::filesystem::copy_file(pwd.string() + "/ffield.gbest." + str_cycle + "." + str_iter + "." + str_parID,
      pwd.string() + "/testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file(pwd.string() + "/geo.val",
    pwd.string() + "/testovfit/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

  if (fixcharges == true){
    boost::filesystem::copy_file(pwd.string() + "/charges_val",
      pwd.string() + "/testovfit/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file(pwd.string() + "/control", pwd.string() + "/testovfit/control",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::copy_file(pwd.string() + "/valset.in", pwd.string() + "/testovfit/trainset.in",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::ofstream iopt_file(pwd.string() + "/testovfit/fort.20");
  iopt_file << "0";
  iopt_file.close();
  ofstream outfile35(pwd.string() + "/testovfit/fort.35");
  outfile35 << "23434.1" << endl;
  outfile35.close();

  // cd to testovfit directory
  string old_path = pwd.string();
  boost::filesystem::path p(pwd.string() + "/testovfit");
  boost::filesystem::current_path(p);
  if (lg_yn == true) {
    // execute reac or reac
    boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/testovfit/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
    system("./reac > /dev/null 2>&1");
  } else {
    boost::filesystem::copy_file(pwd.string() + "/reac", pwd.string() + "/testovfit/reac",
    boost::filesystem::copy_option::overwrite_if_exists);
    system("./reac > /dev/null 2>&1");
  };

  // cd back to main directory
  boost::filesystem::path p2(old_path);
  boost::filesystem::current_path(p2);

  // read fitness value contained in newly generated fort.13 file
  string str;
  boost::filesystem::ifstream myfile("testovfit/fort.13");
  stringstream tempstr;
  tempstr.str("");
  getline(myfile, str);
  // insert str into stringstream tempstr
  tempstr << str;
  // get rid from extra whitespace in stringstream
  tempstr >> std::ws;
  // insert back to str
  tempstr >> str;
  // check if fitness is numeric or ******
  if (str.at(0) == '*') {
    currovfitness = numeric_limits < double > ::infinity();
  } else {
    // convert to double
    currovfitness = stod(str);
  };
  myfile.close();
  boost::filesystem::remove( "testovfit/fort.13" );
  if (firstovfit == true) {
  //  cout << "firstovfit is true" << endl;
  ovfitness = currovfitness;
  };
  //cout << "firstovfit on CPU: " << core << " is: " << firstovfit << endl;

    //cout << "current over-fitting fitness: " << currovfitness << endl;
    //cout << "previous over-fitting fitness: " << ovfitness << endl;

  //if (currovfitness > ovfitness) {
    //cout << "\nOverfitting detected! current fitness: " << currovfitness << "> previous fitness: " << ovfitness << endl;
    //cout << "Early-stopping ReaxFF optimization.\n";
    //cout << "Training completed successfuly!\n";
    //cout << "Total ReaxFF calls: " << funceval << endl;
    //cout << "Global best ReaxFF fit: " << newSwarm.get_gbfit() << endl;
    //cout << "Global best ReaxFF parameters:\n";
    //cout << "[ ";
    //for (int m = 0; m < dim; m++) {
    //  cout << newSwarm.get_gbpos().at(m) << " ";
    //};
    //cout << "]" << endl;
    //ierr = MPI_Finalize();
    //exit(EXIT_SUCCESS);
    ofstream outfileovfit("overfit.out." + std::to_string(cycle), ofstream::app);
    stringstream ss3;
    ss3 << boost::format("%5i %25.4f") %iter %currovfitness;
    outfileovfit << ss3.str();
    ss3.str("");
    outfileovfit << endl;
    outfileovfit.close();
  //}else{
   //ovfitness = currovfitness;
  //};
#endif
};

vector < double > Swarm::get_gbpos() {
  return gbpos;
};

double Swarm::get_gbfit() {
  return gbfit;
};

void Swarm::set_gbfit(double fit) {
  gbfit = fit;
};

void Swarm::update_gbpos(Par & newPar) {
  gbpos.clear();
  for (int j = 0; j < dim; j++) {
    gbpos.push_back(newPar.get_bpos(j));
  };
};

void Swarm::printdisp(Swarm & newSwarm, int iter, int cycle, int fr) {
#ifdef WITH_MPI
 ofstream outfiledisp("disp_log.out." + std::to_string(cycle), ofstream::app);
 if (mod(iter, fr) == 0.0) {
   stringstream ss;
   ss << boost::format("%5i %15.4f") %iter %(newSwarm.get_disp(newSwarm)/initial_disp);
   if (core == 0){
     outfiledisp << ss.str();
     outfiledisp << endl;
     outfiledisp.close();
   };
 };
#endif
#ifndef WITH_MPI
 ofstream outfiledisp("disp_log.out." + std::to_string(cycle), ofstream::app);
 if (mod(iter, fr) == 0.0) {
   stringstream ss;
   ss << boost::format("%5i %15.4f") %iter %(newSwarm.get_disp(newSwarm)/initial_disp);
   outfiledisp << ss.str();
   outfiledisp << endl;
 };
 outfiledisp.close();
#endif
};

void Swarm::printopt(Swarm & newSwarm, int iter, int cycle, int fr) {
  ofstream outfileopt("opti_log.out." + std::to_string(cycle), ofstream::app);
  if (mod(iter, fr) == 0.0) {
    stringstream ss;
    ss << boost::format("%5i %25.4f") %iter %newSwarm.get_gbfit(); 
    outfileopt << ss.str();
    outfileopt << endl;
    
  };
  outfileopt.close();
};


void Swarm::printUQQoI(Swarm & newSwarm, int iter, int cycle, int fr) {
  string str_core = std::to_string(core);
  string line;
  // define temp sums of partial errors from fort.99 file for each core
  double serr_bonds = 0.0;
  double serr_angles = 0.0;
  double serr_torsions = 0.0;
  double serr_energies = 0.0;
  double serr_charges = 0.0;
  double serr_heats = 0.0;
  std::smatch m;
  std::regex e1 ("Heat of formation:");
  std::regex e2 ("Bond distance:");
  std::regex e3 ("Valence angle:");
  std::regex e4 ("Charge atom:");
  std::regex e5 ("Torsion angle:");
  std::regex e6 ("Energy");

  vector <string> values; 

  if (mod(iter, fr) == 0.0) {
    for (int p = 0; p < NumP; p++){
      boost::filesystem::ifstream myfile("CPU." + str_core +"/fort.99");
      // only include swarm members whose fitness is not higher than (gbfit + 15%)
      if (newSwarm.GetPar(p).get_fitness() <= 1.15*gbfit){
        // for each line
        while (getline(myfile, line)) {
          if (regex_search (line,m,e1)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQHEATS.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQHEATS("UQHEATS." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQHEATS << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQHEATS.close();
            continue;
          };
          if (regex_search (line,m,e2)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQBONDS.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQBONDS("UQBONDS." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQBONDS << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQBONDS.close();
            continue;
          };
          if (regex_search (line,m,e3)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQANGLES.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQANGLES("UQANGLES." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQANGLES << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQANGLES.close();
            continue;
          };
          if (regex_search (line,m,e4)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQCHARGES.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQCHARGES("UQCHARGES." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQCHARGES << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQCHARGES.close();
            continue;
          };
          if (regex_search (line,m,e5)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQTORSIONS.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQTORSIONS("UQTORSIONS." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQTORSIONS << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQTORSIONS.close();
            continue;
          };
          if (regex_search (line,m,e6)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQENER.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQENER("UQENER." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQENER << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQENER.close();
            continue;
          };
        }; // done reading all lines
      };
      myfile.close();
    }; // done loop on members
  };
};

void Swarm::printUQFF(Swarm & newSwarm, int iter, int cycle, int fr) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  if (mod(iter, fr) == 0.0) {
    stringstream ss1;
    for (int p = 0; p < NumP; p++){
      // only include swarm members whose fitness is not higher than (gbfit + 15%)
      if (newSwarm.GetPar(p).get_fitness() <= 1.15*gbfit){
          for (int m = 0; m < dim; m++){
            // UQFF.* file will be generated for each m dimension (ffield parameter) for every cycle
            ofstream outfileUQFF("UQFF." + std::to_string(m) + "." + std::to_string(cycle), ofstream::app);
            ss1 << boost::format("%10.4f") %newSwarm.GetPar(p).get_pos_vec().at(m);
            outfileUQFF << ss1.str() << endl;
            ss1.str("");
            outfileUQFF.close();
          };
      };
    };
  };
};


void Swarm::printpos(Swarm & newSwarm, int iter, int cycle, int fr) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
#ifdef WITH_MPI
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core + "/pos");
  const char * path = "/pos/pos_log.out.";
  ofstream outfilepos("CPU." + str_core + path + std::to_string(cycle), ofstream::app);
#endif
#ifndef WITH_MPI
  boost::filesystem::create_directory("pos");
  const char * path = "pos/pos_log.out.";
  ofstream outfilepos(path + std::to_string(cycle), ofstream::app);
#endif
  if (mod(iter, fr) == 0.0) {
    outfilepos << "#Timestep: " << iter << endl;
    stringstream ss1;
    stringstream ss2;
    stringstream ss3;
    for (int p = 0; p < NumP; p++) {
      ss1 << boost::format("%3i") %p;
      outfilepos << ss1.str();
      ss1.str("");
      for (int m = 0; m < dim; m++) {
        ss2 << boost::format("%10.4f") %newSwarm.GetPar(p).get_pos(m); 
        outfilepos << ss2.str();
        ss2.str("");
      };
      ss3 << boost::format("%25.4f") %newSwarm.GetPar(p).get_fitness();
      outfilepos << ss3.str();
      ss3.str("");
      outfilepos << endl;
    };
  };
  outfilepos.close();
};

void Swarm::printvel(Swarm & newSwarm, int iter, int cycle, int fr) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
#ifdef WITH_MPI
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core + "/vel");
  const char * pathvel = "./vel/vel_opti.out.";
  ofstream outfilevel(pathvel + std::to_string(cycle), ofstream::app);
#endif
#ifndef WITH_MPI
  boost::filesystem::create_directory("vel");
  const char * pathvel = "./vel/vel_opti.out.";
  ofstream outfilevel(pathvel + std::to_string(cycle), ofstream::app);
#endif
  if (mod(iter, fr) == 0.0) {
    outfilevel << "#iter: " << iter << endl;

    for (int n = 0; n < NumP; n++) {
      outfilevel << n << "  ";
      for (int q = 0; q < dim; q++) {
        outfilevel << newSwarm.GetPar(n).get_vel(q) << "  ";
      };
      outfilevel << endl;
    };
  };
  outfilevel.close();
};

void Swarm::printdeg(Swarm & newSwarm, int iter, int cycle, int fr) {
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core + "/deg");
  const char * pathdeg = "./deg/deg_opti.out.";
  ofstream outfiledeg(pathdeg + std::to_string(cycle), ofstream::app);
#endif
#ifndef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::create_directory("deg");
  const char * pathdeg = "./deg/deg_opti.out.";
  ofstream outfiledeg(pathdeg + std::to_string(cycle), ofstream::app);
#endif
  if (dim >= 2) {
    if (mod(iter, fr) == 0.0) {
      outfiledeg << "#iter: " << iter << endl;

      for (int n = 0; n < NumP; n++) {
        outfiledeg << n << "  ";
        // print out the angle (deg) of each particle with respect to x-axis
        double deg = 0.0;
        deg = (180.0 / pi) * atan2(newSwarm.GetPar(n).get_vel(1), newSwarm.GetPar(n).get_vel(0));
        outfiledeg << deg << "  ";
        outfiledeg << endl;

      };
    };
  };
  outfiledeg.close();
};
