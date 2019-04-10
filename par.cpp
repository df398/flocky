/*  RiPSOGM v 1.1 Copyright (C) 2019 David Furman, PhD. df398@cam.ac.uk
    Department of Chemistry, University of Cambridge, UK.
    
    RiPSOGM: 
    Rotation Invariant Particle Swarm Optimization with Gaussian
    Mutations.
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
*/

#include "par.h"

std::random_device seed; // initialize non-determinstic seed
std::mt19937 generator(seed() + core); // Mersene Twister: Good quality random number generator
std::uniform_real_distribution < double > dist1(0.0, 1.0); // Type of distribution to sample from and given range

int ierr, core, numcores = 0;
int funceval = 0; // number of function evaluations
double curfit;
//bool inside = true;							// tag if particle position in dimension lies inside [mindomain, maxdomain]
int dim = 0; // initialization for extern variable 'int dim' declared in par.h
int maxtime = 0; // initialization for extern variable 'int maxtime' declared in par.h
int NumP = 0; // initialization for extern variable 'int NumP' declared in par.h
int freq = 1;
int iter = 0;
int repeatn = 1;
int timestep = 0;
bool lg_yn = false;
bool contff = false;
bool chang = false;
double c1 = 2.0; // constant 1
double c2 = 2.0; // constant 2
double inertiafac = 0.9;
double inertiamax = 0.9; // initial inertia weight
double inertiamin = 0.4; // final inertia weight
double levyscale = 1.0;
double confac = 1.0; //0.72984;							// constriction factor multiplier. Shi and Eberhart (2000)
int faili = 1;

int mod(int x, int m) { // definition of modulo function
  return (x % m + m) % m;
}

double l2_norm(vector < double >
  const & u) { // return L2-norm (magnitude) of a vector. uses a range-based loop.
  double accum = 0.;
  for (double x: u) {
    accum += x * x;
  }
  return sqrt(accum);
};

Par::Par() {
  // set min/max domains from params.mod file and set dim = numlines in params.mod file
  read_bounds();

  // read ffield file into matrix ffieldmat. split by each entry.
  read_ffield();

  for (int m = 0; m < dim; m++) {
    std::uniform_real_distribution < double > dist2(mindomain.at(m), maxdomain.at(m));
    //	std::normal_distribution <double> distnorm(0.0,20.0);
    double x = dist2(generator); // uniform positions initialization
    double v = 0.5 * (dist2(generator) - x);
    pos.push_back(x); // initialize particle's position vector
    vel.push_back(v); // initialize particle's velocity vector
    bpos = pos; // initialize particle's best own position vector
  };
};

void Par::read_ffield() {
  // read ffield file into ffieldmat matrix split by entries
  // following: https://stackoverflow.com/questions/10521658/reading-files-columns-into-array
  ffieldmat.clear();
  string str_core = std::to_string(core);
  std::ifstream fin("CPU." + str_core + "/ffield");
  if (fin.fail()) {
    cout << "Unable to open 'ffield' file on CPU" << core << " \n";
    exit(EXIT_FAILURE);
  } else {
    std::string line;
    int numlines = 0;

    while (std::getline(fin, line)) { // for each line
      numlines++;
      std::vector < string > lineData; // create a new row
      string val;
      std::istringstream lineStream(line);
      while (lineStream >> val) { // for each value in line
        lineData.push_back(val); // add to the current row
      };
      ffieldmat.push_back(lineData); // add row to ffieldmat matrix

    };
  };
};

void Par::write_ffield() {
  // update particle ffield
  int index = 0;
  for (int line: ffline) { // for each line value in ffline and corresponding column value in ffcol cp corresponding pos value into ffieldmat. Range-based for statement --> Following: https://msdn.microsoft.com/en-us/library/jj203382.aspx
    char buffer[50]; // buffer for sprintf
    double n;
    n = sprintf(buffer, "%9.4f", pos.at(index));
    ffieldmat[line][ffcol.at(index)] = buffer;
    index = index + 1;

  };
  string str_core = std::to_string(core);
  // write updated ffield to ffield.tmp file
  ofstream output_file;
  // current ffield file stream
  ifstream ffield_file;
  string comment;
  output_file.open("CPU." + str_core + "/ffield.tmp", ios::out);
  ffield_file.open("CPU." + str_core + "/ffield", ios:: in );
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
  /* 
          06.03.19 Check how many elements are there in the ffield then limit the write of atom params
          to only those lines that belong to number of elements. This is to avoid writing lines
          that belong to the next section (bonds) in the atom section.
        */
  using boost::is_any_of;
  string numel_line = ffield_lines.at(41); // the line that contains the number of elements
  vector < string > results; // vector to store the words after split
  boost::trim(numel_line);
  boost::split(results, numel_line, is_any_of("\t "));
  int numel = stoi(results.at(0));
  int max_line_atompar = 4 * numel + 45;
  /* -------------------------------------
   * WRITE ATOM PARAMS SECTION
   * -------------------------------------
   */

  for (int m = 45; m < 46; m++) {
    if (m < max_line_atompar) {
      //cout << "m: " << m << endl;
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 46; m < 49; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 49; m < 50; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 50; m < 53; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 53; m < 54; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 54; m < 57; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 57; m < 58; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 58; m < 61; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 61; m < 62; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 62; m < 65; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 65; m < 66; m++) {
    if (m < max_line_atompar) {
      boost::format f("% 2s%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 66; m < 69; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 69; m < 70; m++) {
    if (m < max_line_atompar) {
      boost::format f("% s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };
  for (int m = 70; m < 73; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
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
  string numbty_line = ffield_lines.at(max_line_atompar); // the line the contains the number of bond types
  vector < string > results_bonds; // vector to store the words after split
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

  /* 
    06.03.19 Check how many elements are there in the ffield then limit the write of bond params
    to only those lines that belong to number of bond types. This is to avoid writing lines
    that belong to the next section (off-dia) in the bond section.
  */
  using boost::is_any_of;
  string numodty_line = ffield_lines.at(max_line_bondpar + 2); // the line that contains the number of off-diag types
  //cout << "numodty_line: " << numodty_line << endl;
  vector < string > results_offdiag; // vector to store the words after split
  boost::trim(numodty_line);
  boost::split(results_offdiag, numodty_line, is_any_of("\t "));
  int numodty = stoi(results_offdiag.at(0));
  //cout << "results_offdiag.at(0):" << results_offdiag.at(0) << endl;
  int max_line_offdpar = max_line_bondpar + 3 + numodty;

  for (int m = max_line_bondpar + 3; m < max_line_offdpar; m++) {
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f"); // the last entry is 10.4f because dispersion coeff. can get to 4 digits long. So, to prevent it sticking to the left column
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
  string numaty_line = ffield_lines.at(max_line_offdpar); // the line that contains the number of angle types
  //cout << "numaty_line: " << numaty_line << endl;
  vector < string > results_angle; // vector to store the words after split
  boost::trim(numaty_line);
  boost::split(results_angle, numaty_line, is_any_of("\t "));
  int numaty = stoi(results_angle.at(0));
  //cout << "numaty:" << numaty << endl;
  //cout << "max_line_offdpar+1" << max_line_offdpar+2;
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
  string numtoty_line = ffield_lines.at(max_line_angles); // the line that contains the number of torsion types
  //cout << "numtoty_line: " << numtoty_line << endl;
  vector < string > results_tors; // vector to store the words after split
  boost::trim(numtoty_line);
  boost::split(results_tors, numtoty_line, is_any_of("\t "));
  int numtoty = stoi(results_tors.at(0));
  //cout << "numtoty:" << numtoty << endl;
  //cout << "max_line_offdpar+1" << max_line_offdpar+2;
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
  string numhbty_line = ffield_lines.at(max_line_tors); // the line that contains the number of Hbond types
  //cout << "numhbty_line: " << numhbty_line << endl;
  vector < string > results_hb; // vector to store the words after split
  boost::trim(numhbty_line);
  boost::split(results_hb, numhbty_line, is_any_of("\t "));
  int numhbty = stoi(results_hb.at(0));
  //cout << "numhbty:" << numhbty << endl;
  //cout << "max_line_offdpar+1" << max_line_offdpar+2;
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
  // replace ffield file with the new ffield (ffield.tmp)
  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield.tmp",
    pwd.string() + "/CPU." + str_core + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp");
};

void Par::write_ffield_lg() {
  // update particle ffield
  int index = 0;
  for (int line: ffline) { // for each line value in ffline and corresponding column value in ffcol cp corresponding pos value into ffieldmat. Range-based for statement --> Following: https://msdn.microsoft.com/en-us/library/jj203382.aspx
    char buffer[50]; // buffer for sprintf
    double n;
    n = sprintf(buffer, "%9.4f", pos.at(index));
    ffieldmat[line][ffcol.at(index)] = buffer;
    index = index + 1;

  };
  string str_core = std::to_string(core);
  // write updated ffield to ffield file
  ofstream output_file;
  // current ffield file stream
  ifstream ffield_file;
  string comment;

  output_file.open("CPU." + str_core + "/ffield.tmp", ios::out);
  ffield_file.open("CPU." + str_core + "/ffield", ios:: in );
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

  /* 
          06.03.19 Check how many elements are there in the ffield then limit the write of atom params
          to only those lines that belong to number of elements. This is to avoid writing lines
          that belong to the next section (bonds) in the atom section.
        */
  using boost::is_any_of;
  string numel_line = ffield_lines.at(41); // the line that contains the number of elements
  vector < string > results; // vector to store the words after split
  boost::trim(numel_line);
  boost::split(results, numel_line, is_any_of("\t "));
  int numel = stoi(results.at(0));
  int max_line_atompar = 5 * numel + 45;

  /* -------------------------------------
   * WRITE ATOM PARAMS SECTION
   * -------------------------------------
   */

  for (int m = 45; m < 46; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 46; m < 49; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 49; m < 50; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 50; m < 51; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 51; m < 54; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 54; m < 55; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 55; m < 56; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 56; m < 59; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 59; m < 60; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 60; m < 61; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 61; m < 64; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 64; m < 65; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 65; m < 66; m++) {
    if (m < max_line_atompar) {
      boost::format f(" %s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 66; m < 69; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 69; m < 70; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 70; m < 71; m++) {
    if (m < max_line_atompar) {
      boost::format f("% 2s%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 71; m < 74; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 74; m < 75; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 75; m < 76; m++) {
    if (m < max_line_atompar) {
      boost::format f("% s %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 76; m < 79; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
  };

  for (int m = 79; m < 80; m++) {
    if (m < max_line_atompar) {
      boost::format f("   %9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
      for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    };
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
  string numbty_line = ffield_lines.at(max_line_atompar); // the line the contains the number of bond types
  vector < string > results_bonds; // vector to store the words after split
  boost::trim(numbty_line);
  boost::split(results_bonds, numbty_line, is_any_of("\t "));
  int numbty = stoi(results_bonds.at(0));
  int max_line_bondpar = 2 * numbty + max_line_atompar;

  /* ------------------------------------
   * WRITE BONDS PARAMS SECTION
   * ------------------------------------
   */

  for (int m = max_line_atompar + 2; m < max_line_bondpar + 2; m++) {
    //if (mod(m,2) == 0) {
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
    m = m + 1;
    //} else {
    boost::format f2("      %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f2.exceptions(f2.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));
    for (std::vector < std::string > ::iterator it = ffieldmat.at(m).begin(); it != ffieldmat.at(m).end(); ++it) {
      f2 = f2 % ( * it);
    };
    output_file << f2 << endl;
    //};
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
  string numodty_line = ffield_lines.at(max_line_bondpar + 2); // the line that contains the number of off-diag types
  //cout << "numodty_line: " << numodty_line << endl;
  vector < string > results_offdiag; // vector to store the words after split
  boost::trim(numodty_line);
  boost::split(results_offdiag, numodty_line, is_any_of("\t "));
  int numodty = stoi(results_offdiag.at(0));
  //cout << "results_offdiag.at(0):" << results_offdiag.at(0) << endl;
  int max_line_offdpar = max_line_bondpar + 3 + numodty;

  for (int m = max_line_bondpar + 3; m < max_line_offdpar; m++) {
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%10.4f"); // the last entry is 10.4f because dispersion coeff. can get to 4 digits long. So, to prevent it sticking to the left column
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
  string numaty_line = ffield_lines.at(max_line_offdpar); // the line that contains the number of angle types
  //cout << "numaty_line: " << numaty_line << endl;
  vector < string > results_angle; // vector to store the words after split
  boost::trim(numaty_line);
  boost::split(results_angle, numaty_line, is_any_of("\t "));
  int numaty = stoi(results_angle.at(0));
  //cout << "numaty:" << numaty << endl;
  //cout << "max_line_offdpar+1" << max_line_offdpar+2;
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
  string numtoty_line = ffield_lines.at(max_line_angles); // the line that contains the number of torsion types
  //cout << "numtoty_line: " << numtoty_line << endl;
  vector < string > results_tors; // vector to store the words after split
  boost::trim(numtoty_line);
  boost::split(results_tors, numtoty_line, is_any_of("\t "));
  int numtoty = stoi(results_tors.at(0));
  //cout << "numtoty:" << numtoty << endl;
  //cout << "max_line_offdpar+1" << max_line_offdpar+2;
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
  string numhbty_line = ffield_lines.at(max_line_tors); // the line that contains the number of Hbond types
  //cout << "numhbty_line: " << numhbty_line << endl;
  vector < string > results_hb; // vector to store the words after split
  boost::trim(numhbty_line);
  boost::split(results_hb, numhbty_line, is_any_of("\t "));
  int numhbty = stoi(results_hb.at(0));
  //cout << "numhbty:" << numhbty << endl;
  //cout << "max_line_offdpar+1" << max_line_offdpar+2;
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
  // replace ffield file with the new ffield (ffield.tmp)
  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::copy_file(pwd.string() + "/ffield.tmp", pwd.string() + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  boost::filesystem::remove(pwd.string() + "/ffield.tmp");
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

    while (std::getline(fin, line)) { // for each line
      numlines++;
      std::vector < double > lineData; // create a new row
      double val;
      std::istringstream lineStream(line);
      while (lineStream >> val) { // for each value in line
        lineData.push_back(val); // add to the current row
      };
      allData.push_back(lineData); // add row to allData
    };

    dim = numlines;

    for (int i = 0; i < dim; i++) {
      mindomain.push_back(allData[i][5]); // read bounds from modified params file
      maxdomain.push_back(allData[i][6]);
      ffline.push_back(allData[i][0] - 1); // read line number of parameter from modified params file
      ffcol.push_back(allData[i][1] - 1); // read column number of parameter from modified params file
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
  // see: http://mathworld.wolfram.com/HyperspherePointPicking.html
  //
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
    dir.at(n) = dir.at(n) / length; // uniformly positioned vector on a hypersphere
  };

  return dir;
};

double Par::get_levy_McCul(double time, double maxtime) {
  // Generation of Levy symmetric distribution (beta=0) using
  // McCulloch's Algorithm
  // When alpha approaches 2, the distribution becomes Gaussian with mean tau
  // and variance 2 * c ^ 2, and beta has no effect.

  //std::uniform_real_distribution <double> dist1(0.0, 1.0);				// uniform distribution
  std::normal_distribution < double > distnorm(0.0, 1.0);
  //double w = -log(dist1(generator));							// generate variable from exponential distribution
  //double phi = (dist1(generator) - 0.5)*pi;						// generate variable from uniform distribution
  //double alpha_init = 1.0;
  //double alpha_fin = 2.0;
  //double alpha = 2.0; // alpha_init + time * (alpha_fin - alpha_init) / maxtime; //min(2.0,exp(log(2.0)/10000.0*time)); // alpha_init + time * (alpha_fin - alpha_init) / maxtime;
  //double c = 1.0; // 0.01 - (0.01 - 1e-50)*time / maxtime; // 0.0003*l2_norm(space_range) - 0.0066; // 0.0003*l2_norm(space_range) - 0.0069; // 0.1 - (0.1 - 0.001)*(time / maxtime);														// scaling parameter
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

void Par::update_vel(double inertiaf, double CF, vector < double > globpos, double time) {

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

void Par::update_pos_levy(vector < double > globpos, double time, double inertiaf) {
  //std::normal_distribution <double> normdist2(0.0, 1.0);
  double levystep = abs(get_levy_McCul(time, maxtime));
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

double Par::eval_fitness(int parid) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string parID = std::to_string(parid);
  string str_core = std::to_string(core);

  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo." + parID,
    pwd.string() + "/CPU." + str_core + "/geo", boost::filesystem::copy_option::overwrite_if_exists);

  // check if ffield is LG or not. execute correct reac accordingly.
  if (lg_yn == true) {
    write_ffield_lg();
    std::ifstream fin5(("CPU." + str_core + "/reac_lg").c_str());
    if (fin5.fail()) {
      cout << "reac_lg executable not found for CPU " + str_core + ". Aborting! \n";
      exit(EXIT_FAILURE);
    }

    // prepare mandatory files before executing reac_lg in each CPU directory
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo",
      pwd.string() + "/CPU." + str_core + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield",
      pwd.string() + "/CPU." + str_core + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);

    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/charges",
      pwd.string() + "/CPU." + str_core + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);

    // cd to each CPU.x directory
    string old_path = pwd.string();
    boost::filesystem::path p(pwd.string() + "/CPU." + str_core);
    boost::filesystem::current_path(p);
    // execute reac within each CPU.x directory
    boost::process::system("./reac_lg", boost::process::std_out > boost::process::null, boost::process::std_err > stderr);
    // cd back to main directory
    boost::filesystem::path p2(old_path);
    boost::filesystem::current_path(p2);

  } else {
    write_ffield();
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

    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/charges",
      pwd.string() + "/CPU." + str_core + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);

    double st2, et2, df2, smf2, avf2;
    st2 = MPI_Wtime();

    // cd to each CPU.x directory
    string old_path = pwd.string();
    boost::filesystem::path p(pwd.string() + "/CPU." + str_core);
    boost::filesystem::current_path(p);
    // execute reac within each CPU.x directory
    boost::process::system("./reac", boost::process::std_out > boost::process::null, boost::process::std_err > stderr);
    // cd back to main directory
    boost::filesystem::path p2(old_path);
    boost::filesystem::current_path(p2);

    et2 = MPI_Wtime();
    df2 = et2 - st2;
    cout << "CPU time of reac: " << df2 << endl;

  };

  // read fitness value contained in fort.13 file
  string str;
  boost::filesystem::ifstream myfile("CPU." + str_core + "/fort.13");
  stringstream tempstr;
  getline(myfile, str);
  tempstr << str; // insert str into stringstream tempstr
  tempstr >> std::ws; // get rid from extra whitespace in stringstream
  tempstr >> str; // insert back to str
  // check if fitness is numeric or ****** 
  if (str.at(0) == '*') {
    fitness = numeric_limits < double > ::infinity();
  } else {
    fitness = stod(str); // convert to double  
  };
  boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/fort.13");

  //	// saving relevant files
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/summary.txt");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/moldyn.vel");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/run.log");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/xmolout");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/bond*");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/molfra.out");		
  //
  //	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/dipole.out" ,
  //	pwd.string()+"/CPU."+str_core+"/current_dipole.out",boost::filesystem::copy_option::overwrite_if_exists);
  //	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
  //	pwd.string()+"/CPU."+str_core+"/current_molgeo.out",boost::filesystem::copy_option::overwrite_if_exists);
  //	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.74" ,
  //	pwd.string()+"/CPU."+str_core+"/current_thermo.out",boost::filesystem::copy_option::overwrite_if_exists);
  //	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.99" ,
  //	pwd.string()+"/CPU."+str_core+"/current_results.out",boost::filesystem::copy_option::overwrite_if_exists);
  //	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
  //	pwd.string()+"/CPU."+str_core+"/geo."+parID,boost::filesystem::copy_option::overwrite_if_exists);
  //	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.73" ,
  //	pwd.string()+"/CPU."+str_core+"/current_partialE.out",boost::filesystem::copy_option::overwrite_if_exists);
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/fort.*");

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
      //cout << "particle id: " << p << endl;
      //cout << "fitness: " << fit << endl;
    }
  };
  return particle_id;
};

void Swarm::Populate(Swarm & newSwarm, int iter) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  string itercount = std::to_string(iter);

  if (core == 0) {
    cout << "\n";
    cout << "Swarm generation started. Please wait." << endl;
  };

  // initial gbfit is INF for all processes
  gbfit = numeric_limits < double > ::infinity();

  double starttime, endtime, diff, sumdiff, avgdiff;
  starttime = MPI_Wtime();

  for (int p = 0; p < NumP; p++) {
    string parID = std::to_string(p);
    // cp geo to geo.parID so each particle works with its own geo file
    boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/geo", pwd.string() + "/CPU." + str_core + "/geo." + parID, boost::filesystem::copy_option::overwrite_if_exists);

    Par NewPar;
    newSwarm.AddPar(NewPar);

    if (core == 0) {
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
    };

    // Initialize global best fitness and corresponding global best position and personal best
    curfit = newSwarm.GetPar(p).eval_fitness(p);

    //cout << "timing of eval_fit: " << curfit << " in Populate for CPU: " << core << " is: " << diff << endl;

    newSwarm.GetPar(p).set_fitness(curfit);
    newSwarm.GetPar(p).set_bfit(curfit);

    if (curfit < gbfit) {
      gbfit = curfit;
    };

    // define a struct (pair) to hold the minimum fitness across processes and its rank for particle p
    struct {
      double tmp_fit;
      int tmp_cpu;
    }
    min_vals_in[1], min_vals_out[1];

    min_vals_in[0].tmp_fit = gbfit; // store current fit on each process
    min_vals_in[0].tmp_cpu = core; // store core id of that current process

    // get minimum fitness across all processes and the corresponding core id and store them in min_vals_out
    MPI_Allreduce( & min_vals_in, & min_vals_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    // global best fitness across all processes for particle p
    gbfit = min_vals_out[0].tmp_fit;
    // core rank the above fitness came from
    cpuid_gbfit = min_vals_out[0].tmp_cpu;
    // store particle id the above fitness came from in parid_gbfit
    parid_gbfit = p;

    // broadcast the global best fitness data: gbfit, cpuid_gbfit and parid_gbfit
    MPI_Bcast( & cpuid_gbfit, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( & parid_gbfit, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( & gbfit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    gbpos.clear();
    gbpos = newSwarm.GetPar(parid_gbfit).get_pos_vec();
    MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, MPI_COMM_WORLD);

  }; // done with all particles

  endtime = MPI_Wtime();
  diff = endtime - starttime;
  MPI_Reduce( & diff, & sumdiff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  avgdiff = sumdiff / numcores;

  if (core == 0) {
    cout << "Time in Populate NumP: " << avgdiff << " seconds" << endl;
  };

  //cout << "gbfit over all particles: " << gbfit << " on CPU " << core << endl;
  if (core == 0) {
    cout << "Swarm generation completed." << endl;
    cout << "Initial global best fit: " << gbfit << endl;
    cout << "RiPSOGM optimization started!" << endl;
  };
};

void Swarm::Propagate(Swarm & newSwarm, int iter) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);

  boost::filesystem::ofstream outfile("opti_log.out." + std::to_string(iter));
  outfile << "#Timestep, Global best fitness" << endl;
  timestep = 0;
  funceval = 0;
  //if (core == 0) {cout << "Progress:" << endl;};

  for (int time = 0; time < maxtime; time++) {
    //	if (core == 0){
    //		cout << time << endl;
    //	};
    inertiafac = inertiamax - time * (inertiamax - inertiamin) / maxtime;
    // Print out particle information
    //if (core == 0){
    //newSwarm.printpos(newSwarm, timestep, iter, freq);
    //newSwarm.printvel(newSwarm, timestep, iter, freq);
    //newSwarm.printdeg(newSwarm, timestep, iter, freq);
    //};
    for (int p = 0; p < NumP; p++) {
      // Update velocities and positions
      newSwarm.GetPar(p).update_vel(inertiafac, confac, newSwarm.get_gbpos(), time);
      if (newSwarm.GetPar(p).fails > faili) {
        newSwarm.GetPar(p).update_pos_levy(gbpos, time, inertiafac);
        newSwarm.GetPar(p).fails = 0;
      } else {
        newSwarm.GetPar(p).update_pos();
      };

      curfit = newSwarm.GetPar(p).eval_fitness(p);
      newSwarm.GetPar(p).set_fitness(curfit);
      funceval = funceval + 1;

      // Update personal best positions and fitness
      if (newSwarm.GetPar(p).get_fitness() < newSwarm.GetPar(p).get_bfit()) {
        newSwarm.GetPar(p).update_bpos();
        newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());
        newSwarm.GetPar(p).fails = 0;

        if (newSwarm.GetPar(p).get_bfit() < gbfit) {
          gbfit = newSwarm.GetPar(p).get_bfit();
          gbpos.clear();
          gbpos = newSwarm.GetPar(p).get_pos_vec();
        };

      } else {
        newSwarm.GetPar(p).fails = newSwarm.GetPar(p).fails + 1;
      };

      // define a struct to hold the global best fit and cpuid in each process for particle p
      struct {
        double tmp_fit;
        int tmp_cpu;
      }
      min_vals_in[1], min_vals_out[1];
      // store best fit on each process
      min_vals_in[0].tmp_fit = gbfit;
      // store core id of that best fitness process
      min_vals_in[0].tmp_cpu = core;
      // get global best fitness across all best fitnesses over processes and the corresponding core id and store them in min_vals_out
      MPI_Reduce( & min_vals_in, & min_vals_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
      MPI_Bcast( & min_vals_out, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);
      gbfit = min_vals_out[0].tmp_fit;
      cpuid_gbfit = min_vals_out[0].tmp_cpu;
      // broadcast global best position from cpuid_gbfit to everybody
      MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, MPI_COMM_WORLD);
      // write gbest data
      if (core == cpuid_gbfit) {
        write_ffield_gbest();
      };
    }; // done loop on particles
    // write opti file from core == 0
    if (core == 0) {
      outfile << timestep << " " << newSwarm.get_gbfit() << endl;
    };
    newSwarm.printpos(newSwarm, timestep, iter, freq);
    timestep = timestep + 1;
  }; // done loop on iterations
  outfile.close();

  MPI_Barrier(MPI_COMM_WORLD);
  if (core == 0) {
    cout << "Optimization completed successfuly!\n";
    cout << "Global best ReaxFF fit: " << newSwarm.get_gbfit() << endl;
    cout << "Total #function eval: " << funceval << endl;
    cout << " " << endl;
    cout << "Global best ReaxFF parameters:" << endl;
    cout << "[ ";
    for (int m = 0; m < dim; m++) {
      cout << newSwarm.get_gbpos().at(m) << " ";
    };
    cout << " ]" << endl;
  };
  //cout << "Total skipped particles: " << maxtime*NumP - funceval << endl;
  //	outfile.close();

  // general clean-up
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_results.out");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_thermo.out");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_molgeo.out");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_partialE.out");
  //	boost::filesystem::remove(pwd.string()+"/CPU."+str_core+"/current_dipole.out");
};

void Swarm::write_ffield_gbest() {
  // cp current ffield to be the global best and current analysis files to global best analysis files
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string itercount = std::to_string(iter);
  string str_core = std::to_string(core);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/ffield",
    "ffield.gbest." + itercount, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.99",
    "results.out." + itercount, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.74",
    "thermo.out." + itercount, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.90",
    "molgeo.out." + itercount, boost::filesystem::copy_option::overwrite_if_exists);
  boost::filesystem::copy_file(pwd.string() + "/CPU." + str_core + "/fort.73",
    "partialE.out." + itercount, boost::filesystem::copy_option::overwrite_if_exists);
  //     	boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/current_dipole.out" ,
  //     	"dipole.out"+itercount,boost::filesystem::copy_option::overwrite_if_exists);
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

void Swarm::printpos(Swarm & newSwarm, int timestep, int iter, int fr) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core + "/pos");
  const char * path = "/pos/pos_log.out.";
  ofstream outfilepos("CPU." + str_core + path + std::to_string(iter), ofstream::app);

  if (mod(timestep, fr) == 0.0) {
    outfilepos << "#Timestep: " << timestep << endl;

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

void Swarm::printvel(Swarm & newSwarm, int timestep, int iter, int fr) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core + "/vel");
  const char * pathvel = "./vel/vel_opti.out.";
  ofstream outfilevel(pathvel + std::to_string(iter), ofstream::app);

  //ofstream outfilevel;
  //outfilevel.open("vel.txt", ofstream::app);
  if (mod(timestep, fr) == 0.0) {
    outfilevel << "#timestep: " << timestep << endl;

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

void Swarm::printdeg(Swarm & newSwarm, int timestep, int iter, int fr) {
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  boost::filesystem::create_directory("CPU." + str_core + "/deg");
  const char * pathdeg = "./deg/deg_opti.out.";
  ofstream outfiledeg(pathdeg + std::to_string(iter), ofstream::app);
  if (dim >= 2) {
    if (mod(timestep, fr) == 0.0) {
      outfiledeg << "#timestep: " << timestep << endl;

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
