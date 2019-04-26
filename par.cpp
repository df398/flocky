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
bool fixcharges = false;
bool lg_yn = false;
bool contff = false;
bool chang = false;
bool perc_yn = true;
double perc = 0.2;
double c1 = 2.0;
double c2 = 2.0;
double inertiafac = 0.9;
double inertiamax = 0.9;
double inertiamin = 0.4;
double levyscale = 1.0;
double confac = 1.0;
int faili = 1;

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
    cout << "Unable to open 'ffield' file on CPU" << core << " \n";
#endif
#ifndef WITH_MPI
    cout << "Unable to open 'ffield' file \n"; 
#endif
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
  ffield_file.close();
  // replace ffield file with the new ffield (ffield.tmp.cycle.iter.parid)
  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    pwd.string() + "/CPU." + str_core + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp");
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
  ffield_file.close();
  // replace ffield file with the new ffield (ffield.tmp.cycle.iter.parid)
  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::copy_file(pwd.string() + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID, 
    pwd.string() + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  //boost::filesystem::remove(pwd.string() + "/ffield.tmp");
};

void Par::read_bounds() {
  // Reading params columns into allData array split by entries
  // following: https://stackoverflow.com/questions/10521658/reading-files-columns-into-array

  std::vector < std::vector < double > > allData;
  string str_core = std::to_string(core);
  std::ifstream fin(("CPU." + str_core + "/params.mod").c_str());
  if (fin.fail()) {
    cout << "Unable to open parameters file 'params.mod' on CPU " << core << ". \n";
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

    }else{
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
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_core = std::to_string(core);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);

  // count total # func evaluations
  funceval = funceval + 1;
#ifdef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo." + str_parID,
    pwd.string() + "/CPU." + str_core + "/geo", boost::filesystem::copy_option::overwrite_if_exists);
#endif
#ifndef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/geo." + str_parID,
    pwd.string() + "/geo", boost::filesystem::copy_option::overwrite_if_exists);
#endif
  // check if ffield is LG or not. execute correct reac accordingly.
  if (lg_yn == true) {
    write_ffield_lg(cycle, iter, parid);
#ifdef WITH_MPI
    std::ifstream fin5(("CPU." + str_core + "/reac_lg").c_str());
    if (fin5.fail()) {
      cout << "reac_lg executable not found for CPU " + str_core + ". Aborting! \n";
      exit(EXIT_FAILURE);
    }
#endif
#ifdef WITH_MPI
    // prepare mandatory files before executing reac_lg in each CPU directory
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo",
      pwd.string() + "/CPU." + str_core + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

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
#endif
    // execute reac within each CPU.x directory
    boost::process::system("./reac_lg", boost::process::std_out > boost::process::null, boost::process::std_err > stderr);
#ifdef WITH_MPI
    // cd back to main directory
    boost::filesystem::path p2(old_path);
    boost::filesystem::current_path(p2);
#endif
  } else {
    write_ffield(cycle, iter, parid);
#ifdef WITH_MPI
    std::ifstream fin6(("CPU." + str_core + "/reac").c_str());
    if (fin6.fail()) {
      cout << "reac executable not found for CPU " + str_core + ". Aborting! \n";
      exit(EXIT_FAILURE);
    }

    // prepare mandatory files before executing reac in each CPU directory
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo",
      pwd.string() + "/CPU." + str_core + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield",
      pwd.string() + "/CPU." + str_core + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);

    if (fixcharges == true){
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/charges",
      pwd.string() + "/CPU." + str_core + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
    };

    double st2, et2, df2, smf2, avf2;
    st2 = MPI_Wtime();

    // cd to each CPU.x directory
    string old_path = pwd.string();
    boost::filesystem::path p(pwd.string() + "/CPU." + str_core);
    boost::filesystem::current_path(p);
#endif
    // execute reac within each CPU.x directory
    boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > stderr);
#ifdef WITH_MPI
    // cd back to main directory
    boost::filesystem::path p2(old_path);
    boost::filesystem::current_path(p2);

    et2 = MPI_Wtime();
    df2 = et2 - st2;
    //cout << "CPU time of reac: " << df2 << endl;
#endif
  };

  // read fitness value contained in fort.13 file
  string str;
#ifdef WITH_MPI
  boost::filesystem::ifstream myfile("CPU." + str_core + "/fort.13");
#endif
#ifndef WITH_MPI
  boost::filesystem::ifstream myfile("fort.13");
#endif
  stringstream tempstr;
  getline(myfile, str);
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
#ifdef WITH_MPI
  boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/fort.13");
#endif
#ifndef WITH_MPI
  boost::filesystem::remove(pwd.string() + "fort.13");
#endif
  //	// saving relevant files
  /*	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/summary.txt");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/moldyn.vel");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/run.log");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/xmolout");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/bond*");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/molfra.out");		
  */
  	//boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/dipole.out" ,
  	  //pwd.string()+"/CPU."+str_core+"/current_dipole.out",boost::filesystem::copy_option::overwrite_if_exists);
#ifdef WITH_MPI
  	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
  	  pwd.string()+"/CPU."+str_core+"/molgeo.out."+str_cycle+"."+str_iter, 
            boost::filesystem::copy_option::overwrite_if_exists);
  	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.74" ,
  	  pwd.string()+"/CPU."+str_core+"/thermo.out."+str_cycle+"."+str_iter, 
            boost::filesystem::copy_option::overwrite_if_exists);
  	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.99" ,
  	  pwd.string()+"/CPU."+str_core+"/results.out."+str_cycle+"."+str_iter, 
           boost::filesystem::copy_option::overwrite_if_exists);
  	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
  	  pwd.string()+"/CPU."+str_core+"/geo."+str_cycle+"."+str_iter+"."+str_parID, 
           boost::filesystem::copy_option::overwrite_if_exists);
  	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.73" ,
  	  pwd.string()+"/CPU."+str_core+"/partialE.out."+str_cycle+"."+str_iter, 
           boost::filesystem::copy_option::overwrite_if_exists);
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/fort.*");
#endif
  return fitness;
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
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  string itercount = std::to_string(cycle);
#ifdef WITH_MPI
  if (core == 0) {
#endif
    cout << "\n";
    cout << "Swarm generation started. Please wait." << endl;
#ifdef WITH_MPI
  };
#endif
  // initial gbfit is INF for all processes and all swarm members
  gbfit = numeric_limits < double > ::infinity();
  // initial gbpos is 0.0 for all processes and all swarm members
  gbpos.clear();
  for (int m=0; m<dim; m++){
    gbpos.push_back(0.0);
  };
#ifdef WITH_MPI
  double starttime, endtime, diff, sumdiff, avgdiff;
  starttime = MPI_Wtime();
#endif
  // ---------------------------------------------- //
  //     POPULATE: MAIN LOOP OVER SWARM MEMBERS
  // ---------------------------------------------- //
  for (int p = 0; p < NumP; p++) {
    string parID = std::to_string(p);
#ifdef WITH_MPI
    // cp geo to geo.parID so each particle works with its own geo file
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo", 
      pwd.string() + "/CPU." + str_core + "/geo." + parID, 
        boost::filesystem::copy_option::overwrite_if_exists);
#endif
    Par NewPar;
    newSwarm.AddPar(NewPar);
#ifdef WITH_MPI
    if (core == 0) {
#endif
      // If contff == y, then take force field's current values for the position of particle 0 (the rest are just random)
      if (contff == true) {
        vector < double > ffpos;
        ffpos.clear();
        for (int m = 0; m < dim; m++) {
          ffpos.push_back(stod(newSwarm.GetPar(0).ffieldmat.at(newSwarm.GetPar(0).ffline.at(m)).at(newSwarm.GetPar(0).ffcol.at(m))));
        };
        newSwarm.GetPar(0).set_pos(ffpos);
      };
      contff = false;
#ifdef WITH_MPI
    };
#endif
    // evaluate fitness and set personal best equal to fitness
    curfit = newSwarm.GetPar(p).eval_fitness(cycle, 0, p);
    newSwarm.GetPar(p).set_fitness(curfit);
    newSwarm.GetPar(p).set_bfit(curfit);

    //cout << "curfit for particle " << p << " for CPU " << core << " is: " << curfit << endl;
    //cout << "pos for particle " << p << " for CPU " << core << " is: " << endl;
    //for (int m=0; m<dim; m++){
    //  cout << newSwarm.GetPar(p).get_pos(m) << " ";
    //};
    //cout << endl;

    if (curfit < gbfit) {
      gbfit = curfit;
      gbpos.clear();
      gbpos = newSwarm.GetPar(p).get_pos_vec();
#ifndef WITH_MPI
      core = 0;
#endif
      write_ffield_gbest(core, cycle, 0, p);
    };
  }; // done loop over particles
  //cout << "done with all particles" << endl;
#ifdef WITH_MPI
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
  //cout << "gbfit for CPU " << core << " is: " << gbfit << endl;
  // broadcast contents of gbpos vector from rank cpuid_gbfit
  MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, MPI_COMM_WORLD);
  //cout << "gbpos for CPU " << core << " is: " << endl;
  //for (int m=0; m<dim; m++){
  //  cout << newSwarm.get_gbpos().at(m) << " ";
  //};
  //cout << endl;

  //newSwarm.printopt(newSwarm, 0, cycle, 1);

  endtime = MPI_Wtime();
  diff = endtime - starttime;
  MPI_Reduce( & diff, & sumdiff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  avgdiff = sumdiff / numcores;

  if (core == 0) {
    cout << "Time in Populate NumP: " << avgdiff << " seconds" << endl;
  };
#endif
#ifdef WITH_MPI
  if (core == 0) {
#endif
    newSwarm.printopt(newSwarm, 0, cycle, 1);
    cout << "Swarm generation completed." << endl;
    cout << "Initial global best fit: " << gbfit << endl;
    cout << "flocky optimization started!" << endl;
#ifdef WITH_MPI
  };
#endif
};

void Swarm::Propagate(Swarm & newSwarm, int cycle) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);

  // ---------------------------------------------------------------- //
  //         PROPAGATE: main loop over iterations
  // ---------------------------------------------------------------- //
  for (iter = 1; iter < maxiters+1; iter++) {
    inertiafac = inertiamax - iter * (inertiamax - inertiamin) / maxiters;
    //
    // main loop over swarm members
    //
    for (int p = 0; p < NumP; p++) {
      // Update velocities and positions
      newSwarm.GetPar(p).update_vel(inertiafac, confac, newSwarm.get_gbpos(), iter);
      if (newSwarm.GetPar(p).fails > faili) {
        newSwarm.GetPar(p).update_pos_levy(gbpos, iter, inertiafac);
        newSwarm.GetPar(p).fails = 0;
      } else {
        newSwarm.GetPar(p).update_pos();
      };

      curfit = newSwarm.GetPar(p).eval_fitness(cycle, iter, p);
      newSwarm.GetPar(p).set_fitness(curfit);
      //cout << "curfit for particle " << p << " for CPU " << core << " is: " << curfit << endl;
      //cout << "pos for particle " << p << " for CPU " << core << " is: " << endl;
      //for (int m=0; m<dim; m++){
      //  cout << newSwarm.GetPar(p).get_pos(m) << " ";
      //};
      //cout << endl;

      // Update personal best positions and fitness
      if (newSwarm.GetPar(p).get_fitness() < newSwarm.GetPar(p).get_bfit()) {
        newSwarm.GetPar(p).update_bpos();
        newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());
        newSwarm.GetPar(p).fails = 0;
        if (newSwarm.GetPar(p).get_bfit() < gbfit) {
          gbfit = newSwarm.GetPar(p).get_bfit();
          gbpos.clear();
          gbpos = newSwarm.GetPar(p).get_pos_vec();
#ifndef WITH_MPI
          core = 0;
#endif
	  write_ffield_gbest(core, cycle, iter, p);
        };
      } else {
        newSwarm.GetPar(p).fails = newSwarm.GetPar(p).fails + 1;
      };
    }; // done loop over swarm members
#ifdef WITH_MPI
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
    //cout << "gbfit for CPU " << core << " is: " << gbfit << endl;
    // broadcast contents of gbpos vector from rank cpuid_gbfit
    MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, MPI_COMM_WORLD);
    //cout << "gbpos for CPU " << core << " is: " << endl;
    //for (int m=0; m<dim; m++){
    //  cout << newSwarm.get_gbpos().at(m) << " ";
    //};
    //cout << endl;
#endif
#ifdef WITH_MPI
    if (core == 0) {
#endif
      newSwarm.printopt(newSwarm, iter, cycle, freq);
#ifdef WITH_MPI
    };
#endif
    newSwarm.printpos(newSwarm, iter, cycle, freq);
  }; // done loop over iterations
#ifdef WITH_MPI
  MPI_Allreduce(& funceval, & funceval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef WITH_MPI
  if (core == 0) {
#endif
    cout << "Training completed successfuly!\n";
    cout << "Total ReaxFF calls: " << funceval << endl;
    cout << "Global best ReaxFF fit: " << newSwarm.get_gbfit() << endl;
    cout << "Global best ReaxFF parameters:" << endl;
    cout << "[ ";
    for (int m = 0; m < dim; m++) {
      cout << newSwarm.get_gbpos().at(m) << " ";
    };
    cout << " ]" << endl;
#ifdef WITH_MPI
  };
#endif
  //cout << "Total skipped particles: " << maxiters*NumP - funceval << endl;

  // general clean-up
  	/*boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_results.out");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_thermo.out");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_molgeo.out");
  	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_partialE.out");
  	//boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_dipole.out");
       */
};

void Swarm::write_ffield_gbest(int core, int cycle, int iter, int par) {
  //cout << "I'm in write_ffield_gbest!" << endl;
  // cp current ffield to be the global best and current analysis files to global best analysis files
#ifndef WITH_MPI
  boost::filesystem::copy_file(pwd.string() + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    "ffield.gbest." + str_cycle + "." + str_iter+"."+str_parID, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/fort.99",
    "results.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/fort.74",
    "thermo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/fort.90",
    "molgeo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/fort.73",
    "partialE.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string() + "/current_dipole.out" ,
  //  "dipole.out"+itercount,boost::filesystem::copy_option::overwrite_if_exists);
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
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.74",
    "thermo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.90",
    "molgeo.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.73",
    "partialE.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  //boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/current_dipole.out" ,
  //  "dipole.out"+itercount,boost::filesystem::copy_option::overwrite_if_exists);
#endif
};

vector < double > Swarm::get_gbpos() {
  return gbpos;
};

double Swarm::get_gbfit() {
  return gbfit;
};

void Swarm::set_gbfit(double bfit) {
  gbfit = bfit;
};

void Swarm::update_gbpos(Par & newPar) {
  gbpos.clear();
  for (int j = 0; j < dim; j++) {
    gbpos.push_back(newPar.get_bpos(j));
  };
};
void Swarm::printopt(Swarm & newSwarm, int iter, int cycle, int fr) {
  //boost::filesystem::path pwd(boost::filesystem::current_path());
  //string str_core = std::to_string(core);
  //boost::filesystem::create_directory("CPU." + str_core + "/pos");
  //const char * path = "opti_log.out.";
  ofstream outfileopt("opti_log.out." + std::to_string(cycle), ofstream::app);

  if (mod(iter, fr) == 0.0) {
    outfileopt << "#Timestep: " << iter << endl;
    outfileopt << newSwarm.get_gbfit() << " ";
    outfileopt << endl;
  };
  outfileopt.close();
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
  const char * path = "/pos/pos_log.out.";
  ofstream outfilepos(path + std::to_string(cycle), ofstream::app);
#endif
  if (mod(iter, fr) == 0.0) {
    outfilepos << "#Timestep: " << iter << endl;

    for (int p = 0; p < NumP; p++) {
      outfilepos << p << "  ";
      for (int m = 0; m < dim; m++) {
        outfilepos << newSwarm.GetPar(p).get_pos(m) << "  ";
      };
      outfilepos << newSwarm.GetPar(p).get_fitness() << "  ";
      outfilepos << endl;
    };

    outfilepos.close();
  };
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
    outfilevel.close();
  };
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
    outfiledeg.close();
  };
};
