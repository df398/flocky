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

// PCG family of PRNGS: seed with a real random value, if available
pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number generator
pcg32 generator(seed_source);

std::uniform_real_distribution <double> dist1(0.0, 1.0);

int    mycore_swarmcore = 0;
int    size_swarmcores = 0;
int    mycore_reaxffcore = 0;
int    size_reaxffcores = 0;
int    funceval = 0;
double curfit;
double localminfit;
// tag if particle position in dimension lies inside [mindomain, maxdomain]
//bool inside = true;
int    ptrainset=1;
int    totmembers=2;
vector <int> swarmcores(numcores/ptrainset,0);
vector <int> reaxffcores((numcores - numcores/ptrainset),0);
int    dim = 2;
int    NumP = 0;
int    freq = 1;
int    cycle = 0;
int    maxcycles = 1;
int    iter = 0;
int    maxiters = 0;
int    parid_gbfit = 0;
int    localmin = 0;
int    lm_iter_max = 10;
int    numffbags = 0;
double lm_err_tol = 1.0E-6;
bool   fixcharges = false;
bool   lg_yn = false;
bool   contff = false;
bool   verbose = false;
int    regular = 0;
double dropvalue = 0.5;
double hlambda = 0.01;
bool   chang = false;
bool   perc_yn = true;
bool   ofit = false;
bool   uq = false;
bool   gbfitfound = false;
bool   firstovfit = true;
bool   ensembleave = false;
int    preppath = 0;
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
int    faili = 1;
int    ninf = 0;
double econv = 627.50947406;
double frqconv = 2721.2220202;
int    count_lmin_success = 0;
int    lmin_success = 0;
/* define global variables for start
   and end lines of ffield sections
   -------------------------------- */
int    numel = 1;
int    max_line_atompar = 1;
int    numbty = 1;
int    max_line_bondpar = 1;
int    numodty = 1;
int    max_line_offdpar = 1;
int    numaty = 1;
int    max_line_angles = 1;
int    numtoty = 1;
int    max_line_tors = 1;
int    numhbty = 1;
int    max_line_hbs = 1;
vector <vector <int>> inversep;
/* -------------------------------- */

// --------------------------------------------------------------- //
//             General functions definitions                       //
// --------------------------------------------------------------- //


// define int of modulo
int mod(int x, int m) {
  return (x % m + m) % m;
}

// check if str is numeric
bool is_numeric(std::string const& str)
{
    std::string::const_iterator first(str.begin()), last(str.end());
    return boost::spirit::qi::parse(first, last, boost::spirit::double_)
        && first == last;
}

// return L2-norm (magnitude) of a vector. uses a range-based loop.
double l2_norm(vector < double > const & u) {
  double accum = 0.;
  for (double x: u) {
    accum += x * x;
  }
  return sqrt(accum);
};


void get_inversep() {
  vector <int> temp;

  // add line number of inverse parameter from general section
  // p_boc2
  temp = {4, 1};
  inversep.push_back(temp);

  // add line numbers of inverse parameters from atom section
  // r0_sigma
  for (int i = 1; i < max_line_atompar+1 ; i++) {
    temp = {numel +4*i, 2};
    inversep.push_back(temp);
  };
  // r0_pi
  for (int i = 1; i < max_line_atompar+1 ; i++) {
    temp = {numel + 4*i, 8};
    inversep.push_back(temp);
  };
  // r0_pipi
  for (int i = 1; i < max_line_atompar+1 ; i++) {
    temp = {numel + 6*i, 1};
    inversep.push_back(temp);
  };
  // rvdw
  for (int i = 1; i < max_line_atompar+1 ; i++) {
    temp = {numel + 4*i, 5};
    inversep.push_back(temp);
  };
  // gammaw
  for (int i = 1; i < max_line_atompar+1 ; i++) {
    temp = {numel + 5*i, 2};
    inversep.push_back(temp);
  };
  // gammaij
  for (int i = 1; i < max_line_atompar+1 ; i++) {
    temp = {numel + 4*i, 7};
    inversep.push_back(temp);
  };

  // add line numbers of inverse parameters from off-diagonal section
  // Ro
  for (int i = 1; i < max_line_offdpar+1 ; i++) {
    temp = {numodty + i, 4};
    inversep.push_back(temp);
  };
  // rsigma
  for (int i = 1; i < max_line_offdpar+1 ; i++) {
    temp = {numodty + i, 8};
    inversep.push_back(temp);
  };
  // rpi
  for (int i = 1; i < max_line_offdpar+1 ; i++) {
    temp = {numodty + i, 7};
    inversep.push_back(temp);
  };
  // rpipi
  for (int i = 1; i < max_line_offdpar+1 ; i++) {
    temp = {numodty + i, 6};
    inversep.push_back(temp);
  };
  // add line numbers of inverse parameters from hbond section
  for (int i = 1; i < max_line_hbs+1 ; i++) {
    temp = {numhbty + i, 4};
    inversep.push_back(temp);
  };
};


double fitness_wrapper(const vector <double> &x, vector <double> &numgrad, void *data) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered wrapper" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered wrapper" << endl;
};
#endif
  // force-stop the lmin for passive swarm members if all active swarm members completed their lmin
  if (count_lmin_success == swarmcores.size()) {
     throw nlopt::forced_stop();
  };

  // *** test function (2D Rosenbrock. x_opt[1,1] = 0.0. Search domain = [-5,5] should be defined in params.mod
  //double xp, yp;
  //double minfunc;
  //Par * helper = static_cast <Par *> (data);
  //vector <double> mind = helper->mindomain;
  //vector <double> maxd = helper->maxdomain;
  //xp = x.at(0);//*(maxd.at(0) - mind.at(0)) + mind.at(0);
  //yp = x.at(1);//*(maxd.at(1) - mind.at(1)) + mind.at(1);

  //f = 100.0*pow((yp - xp*xp),2) + pow((1.0 - xp),2);
  //if (localmin == 1) {
  //   numgrad.at(0) = 400.0*pow(xp,3) - 400.0*xp*yp + 2.0*xp - 2.0;
  //   numgrad.at(1) = 200.0*(yp-xp*xp);
  //};

  //cout << "WRAPPER x in CPU: " << core << " is: " << x[0] << ", " << x[1] << endl;
  //cout << "WRAPPER f in CPU: " << core << " is: " << boost::format("%10.4f") %minfunc << endl;
  //if (localmin == 1) {
  //   cout << "WRAPPER gradients in CPU: " << core << " is: " << boost::format("%10.4f") %numgrad.at(0) << ", " << boost::format("%10.4f") %numgrad.at(1) << endl; 
  //};
  //cout << endl;

  // *** Tapered ReaxFF cost function
  // helper is a pointer to the calling member (which is an object of the Par class)
  Par * helper = static_cast <Par *> (data);

  if (localmin == 1 && !numgrad.empty()) {
    numgrad = helper->eval_numgrad(x,data);
  };

  double minfunc = helper->eval_fitness(x,data);

  // print standardized data
  //cout << "------------------ flocky iter = " << iter << " ------------------------------------" << endl;
  //cout << "WRAPPER x in CPU: " << core << " is: " << boost::format("%10.4f") %x[0] << ", " << boost::format("%10.4f") %x[1] << endl;
  //cout << "WRAPPER minfunc in CPU: " << core << " is: " << boost::format("%20.10f") %minfunc << endl;
  //if (localmin == 1 && !numgrad.empty()) {
  //   cout << "WRAPPER gradients in CPU: " << core << " is: " << boost::format("%30.10f") %numgrad.at(0) << ", " << boost::format("%30.10f") %numgrad.at(1) << endl;
  //};
  //cout << endl;

  // print physical data
  int cyc = helper->state.cycle;
  #ifdef WITH_MPI
  ofstream lmin_rep("lmin_rep.out." + std::to_string(cyc)+"."+std::to_string(core), ofstream::app);
  #endif
  #ifndef WITH_MPI
  ofstream lmin_rep("lmin_rep.out." + std::to_string(cyc), ofstream::app);
  #endif
  lmin_rep << "------------------ flocky iter = " << iter << " ------------------------------------" << endl;
  vector <double> mind = helper->mindomain;
  vector <double> maxd = helper->maxdomain;
  vector <double> xphys;
  xphys.clear();
  for (int i=0; i < dim; i++) {
     xphys.push_back(x.at(i)*(maxd.at(i) - mind.at(i)) + mind.at(i));
  };
  lmin_rep << "WRAPPER x in CPU: " << core << " is: " << endl;
  for (int i=0; i < dim; i++) {
      lmin_rep << boost::format("%10.4f") %xphys[i] << " ";
  };
  lmin_rep << endl;
  lmin_rep << "WRAPPER minfunc in CPU: " << core << " is: " << boost::format("%1.10e") %minfunc << endl;
  if (localmin == 1 && !numgrad.empty()) {
     lmin_rep << "WRAPPER gradients in CPU: " << core << " is: " << endl;
     for (int i=0; i < dim; i++) {
        lmin_rep << boost::format("%30.10f") %numgrad.at(i) << " ";
     };
  };
  lmin_rep << endl;
  lmin_rep.close();
  return minfunc;
};

// --------------- End general functions definitions -------------- //


Par::Par() {

  if  (ensembleave == false && preppath ==0) {
     // read ffield file into matrix ffieldmat. split by each entry.
     read_ffield();

     // set min/max domains from params.mod file
     read_bounds();

     for (int m = 0; m < dim; m++) {
       // physical domains
       std::uniform_real_distribution < double > dist2(mindomain.at(m), maxdomain.at(m));
       double x = dist2(generator);
       double v = 0.5 * (dist2(generator) - x);
       // initialize particle's position vector
       pos.push_back(x);
       // standardize the position
       pos.at(m) = (pos.at(m) - mindomain.at(m))/(maxdomain.at(m) - mindomain.at(m));
       // initialize particle's velocity vector
       vel.push_back(v);
       // initialize particle's best own position vector
       bpos = pos;
     };
  };

};

void Par::read_ffield() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered read_ffield()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered read_ffield()" << endl;
};
#endif

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
    fin.close();
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
#ifndef WITH_MPI
    cout << "Unable to open 'ffield' file \n"; 
    fin.close();
    exit(EXIT_FAILURE);
#endif
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

void Par::write_geo() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered write_geo()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered write_geo()" << endl;
};
#endif

#ifdef WITH_MPI
  // note: pwd is not portable for mounted partitions that are not accessible directly
  // so it's better to omit it as it's not really needed(?)
  //boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  std::ifstream fin("CPU."+str_core+"/geo");
  if (fin.fail()) {
    cout << "Error: unable to open 'geo' file in CPU" << core << " \n";
    fin.close();
    MPI_Abort(MPI_COMM_WORLD,1);
  } else {
    std::string line;
    int numlines = 0;
    vector <string> keywords;
    vector <string> mykeywords;
    int numstructs = 0;
    int blockID = 0;
    vector <vector <string>> geodata;
    vector <string> geodata_nonsplit;
    struct bgfblock {
       int myblockID;
       string mydescrp;
       int mystart;
       int myend;
       vector <vector <string>> mygeodata;
       vector <string> mygeodata_nonsplit;
    };

    // find total number of structures in geo file
    while (std::getline(fin, line)) {
       boost::split(keywords, line, boost::is_any_of(" "),boost::token_compress_on);
       geodata_nonsplit.push_back(line);
       geodata.push_back(keywords);
       if (keywords.at(0) == "DESCRP")  {
          numstructs++;
       };
    };
    // create a total of #numstructs myblock structs of type 'bgfblock' to hold blocks separately
    bgfblock myblock[numstructs];

    // store geo lines in geodata
    // and for each line in geo find start and end line numbers of blocks
    for (auto line : geodata_nonsplit) {
       numlines++;
       boost::split(keywords, line, boost::is_any_of(" "),boost::token_compress_on);
 
       if (keywords.at(0) == "BIOGRF" || keywords.at(0) == "XTLGRF") {
         myblock[blockID].mystart = numlines;
       };

       if (keywords.at(0) == "DESCRP") {
         myblock[blockID].mydescrp = keywords.at(1);
       };

       if (keywords.at(0) == "END") {
         myblock[blockID].myend = numlines;
         blockID++;
       };
       
    };
    for (int i=0; i < numstructs; i++) {
        for (int j=myblock[i].mystart-1; j < myblock[i].myend; j++) {
            if (!geodata_nonsplit.at(j).empty()) {
               myblock[i].mygeodata_nonsplit.push_back(geodata_nonsplit.at(j));
            };
        };
    };


    // sanity check: blockID should be equal to numstructs
    if (blockID != numstructs) {
#ifdef WITH_MPI
    cout << "Error: blockID " << blockID << " not equal to numstructs " << numstructs << " on CPU" << core << "\n";
    fin.close();
    MPI_Abort(MPI_COMM_WORLD,9);
#endif
#ifndef WITH_MPI
    cout << "Error: blockID not equal to numstructs \n";
    fin.close();
    exit(EXIT_FAILURE);
#endif
    };

    // compose unique DESCRP entries from trainset.in
    string str_core = std::to_string(core);
    std::ifstream fin("CPU."+str_core+"/trainset.in");
    std::string geoline;
    vector <string> traindata;
    vector <string> traindata_nonsplit;
    vector <string> geokeywords;
    if (fin.fail()) {
      cout << "Error: unable to open 'trainset.in' file on CPU" << core << " \n";
      fin.close();
      MPI_Abort(MPI_COMM_WORLD,1);
    } else {
        // store only structure names from trainset.in lines in traindata vector
        geokeywords.clear();
        traindata.clear();
        traindata_nonsplit.clear();
        while (std::getline(fin, geoline)) {
           // ignore lines with zero length or start with a comment
           if (geoline.length() != 0 && geoline[0] != '#') {
             traindata_nonsplit.push_back(geoline);
             // split to separate words with delimeters
             boost::split(geokeywords, geoline, boost::is_any_of("=!@$*&^%:;)(/+- "),boost::token_compress_on);
             // filter out words
             for (int i=0; i < geokeywords.size(); i++) {
               boost::trim(geokeywords.at(i));
               if (geokeywords.at(i) != "CHARGE"      && geokeywords.at(i) != "GEOMETRY"  && 
                   geokeywords.at(i) != "FREQUENCIES" && geokeywords.at(i) != "ENDFREQUENCIES"  &&
                   geokeywords.at(i) != "ENDGEOMETRY" && geokeywords.at(i) != "FORCES"    &&
                   geokeywords.at(i) != "ENDENERGY"   && geokeywords.at(i) != "ENERGY"    &&
                   geokeywords.at(i) != "HEATFO"      && geokeywords.at(i) != "CELL"      &&
                   geokeywords.at(i) != "ENDCHARGE"   && geokeywords.at(i) != "END"       &&
                   geokeywords.at(i) != "ENDFORCES"   && geokeywords.at(i) != "ENDHEATFO" &&
                   geokeywords.at(i) != "PARAMETERS"  && geokeywords.at(i) != "ENDCELL"   &&
                   geokeywords.at(i) != "RMS"         && geokeywords.at(i) != "CRMS"      &&
                   geokeywords.at(i) != "Weight"      && geokeywords.at(i) != "Atom"      &&
                   geokeywords.at(i) != "\n"          && geokeywords.at(i) != ""          &&
                   geokeywords.at(i) != " "           && !is_numeric(geokeywords.at(i))   &&
                   !boost::contains(geokeywords.at(i),"#") && !boost::contains(geokeywords.at(i),",")) {
                   traindata.push_back(geokeywords.at(i));
               };
             };
           };
        };
        // remove duplicate entries in keywords
        std::sort(traindata.begin(), traindata.end());
        traindata.erase(std::unique(traindata.begin(), traindata.end()), traindata.end());
    };
    // print all blocks that pertain to unique trainset entries
    // into respective newgeo files
    ofstream newgeofile;
    newgeofile.open("CPU." + str_core + "/newgeo", ios::out);
    for (int i=0; i < numstructs; i++) {
       if ( find(traindata.begin(), traindata.end(), myblock[i].mydescrp) != traindata.end() ) {
          for (int m=0; m < myblock[i].mygeodata_nonsplit.size(); m++) {
             newgeofile << myblock[i].mygeodata_nonsplit.at(m) << endl;
          };
          newgeofile << endl;
       };
    };
    newgeofile << std::flush;
    newgeofile.close();
    // rename the new geo file
    //boost::filesystem::path pwd(boost::filesystem::current_path());
    boost::filesystem::copy_file("CPU." + str_core + "/newgeo",
    "CPU." + str_core + "/geo", boost::filesystem::copy_option::overwrite_if_exists);
    boost::filesystem::remove("CPU." + str_core + "/newgeo");
  };
#endif
};

void Par::write_trainset() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered write_trainset()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered write_trainset()" << endl;
};
#endif

#ifdef WITH_MPI
  //boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);
  std::ifstream fin("trainset.in");
  if (fin.fail()) {
    cout << "Error: unable to open 'trainset.in' file on CPU" << core << " \n ";
    fin.close();
    MPI_Abort(MPI_COMM_WORLD,1);
  } else {
    std::string line;
    int numlines = 0;
    int charge_begin = 0;
    int charge_end = 0;
    int cell_begin = 0;
    int cell_end = 0;
    int heat_begin = 0;
    int heat_end = 0;
    int geo_begin = 0;
    int geo_end = 0;
    int energy_begin = 0;
    int energy_end = 0;
    int forces_begin = 0;
    int forces_end = 0;
    int freqs_begin = 0; 
    int freqs_end = 0;
    vector < vector <string> > traindata;
    vector <string> traindata_nonsplit;
    vector <string> keywords;

    // store trainset.in lines in traindata
    // and for each line in trainset.in find begin and end line numbers of sections
    keywords.clear();
    traindata.clear();
    traindata_nonsplit.clear();
    while (std::getline(fin, line)) {
      numlines++;
      boost::split(keywords, line, boost::is_any_of(" "));
      traindata_nonsplit.push_back(line);
      traindata.push_back(keywords);

      if (keywords.at(0) == "FREQUENCIES") {
         freqs_begin = numlines;
      };

      if (keywords.at(0) == "ENDFREQUENCIES" || (keywords.at(0) == "END" && keywords.at(1) == "FREQUENCIES")) {
         freqs_end = numlines;
      };


      if (keywords.at(0) == "FORCES") {
         forces_begin = numlines;
      };

      if (keywords.at(0) == "ENDFORCES" || (keywords.at(0) == "END" && keywords.at(1) == "FORCES")) {
         forces_end = numlines;
      };

      if (keywords.at(0) == "CHARGE") {
         charge_begin = numlines;
      };

      if (keywords.at(0) == "ENDCHARGE" || (keywords.at(0) == "END" && keywords.at(1) == "CHARGE")) {
         charge_end = numlines;
      };

      if (keywords.at(0) == "CELL" && keywords.at(1) == "PARAMETERS") {
         cell_begin = numlines;
      };

      if ((keywords.at(0) == "ENDCELL" && keywords.at(1) == "PARAMETERS") || (keywords.at(0) == "END" && keywords.at(1) == "CELL" && keywords.at(2) == "PARAMETERS")) {
         cell_end = numlines;
      };

      if (keywords.at(0) == "HEATFO") {
         heat_begin = numlines;
      };

      if (keywords.at(0) == "ENDHEATFO" || (keywords.at(0) == "END" && keywords.at(1) == "HEATFO")) {
         heat_end = numlines;
      };

      if (keywords.at(0) == "GEOMETRY") {
         geo_begin = numlines;
      };

      if (keywords.at(0) == "ENDGEOMETRY" || (keywords.at(0) == "END" && keywords.at(1) == "GEOMETRY")) {
         geo_end = numlines;
      };

      if (keywords.at(0) == "ENERGY") {
         energy_begin = numlines;
      };

      if (keywords.at(0) == "ENDENERGY" || (keywords.at(0) == "END" && keywords.at(1) == "ENERGY")) {
         energy_end = numlines;
      };
    };
    fin.close();

    // write new trainset.in file with just subsets of data for each swarmcore
    string str_core = std::to_string(core);
    ofstream trainset_file;
    trainset_file.open("CPU." + str_core + "/trainset.in", ios::out);

    // FREQUENCIES SECTION
    if (freqs_end - freqs_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = freqs_begin-1;
       int keyword_end_traindata = freqs_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((freqs_end - freqs_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ( (freqs_end - freqs_begin) - mod((freqs_end - freqs_begin), ptrainset) )/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0;
       int reset_end = 0;
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;

       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((freqs_end - freqs_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ( (freqs_end - freqs_begin) - mod((freqs_end - freqs_begin), ptrainset) )/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((freqs_end - freqs_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ( (freqs_end - freqs_begin) - mod((freqs_end - freqs_begin), ptrainset) )/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    // FORCES SECTION
    if (forces_end - forces_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = forces_begin-1;
       int keyword_end_traindata = forces_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((forces_end - forces_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ((forces_end - forces_begin) - mod((forces_end - forces_begin), ptrainset))/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0;
       int reset_end = 0;
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;
       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((forces_end - forces_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ((forces_end - forces_begin) - mod((forces_end - forces_begin), ptrainset))/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((forces_end - forces_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ((forces_end - forces_begin) - mod((forces_end - forces_begin), ptrainset))/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    // CHARGES SECTION
    if (charge_end - charge_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = charge_begin-1;
       int keyword_end_traindata = charge_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((charge_end - charge_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ((charge_end - charge_begin) - mod((charge_end - charge_begin), ptrainset))/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0; 
       int reset_end = 0;
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;
       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((charge_end - charge_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ((charge_end - charge_begin) - mod((charge_end - charge_begin), ptrainset))/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((charge_end - charge_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ((charge_end - charge_begin) - mod((charge_end - charge_begin), ptrainset))/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    // CELL PARAMETERS SECTION
    if (cell_end - cell_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = cell_begin-1;
       int keyword_end_traindata = cell_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((cell_end - cell_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ((cell_end - cell_begin) - mod((cell_end - cell_begin), ptrainset))/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0;
       int reset_end = 0;
 
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;

       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((cell_end - cell_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ((cell_end - cell_begin) - mod((cell_end - cell_begin), ptrainset))/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((cell_end - cell_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ((cell_end - cell_begin) - mod((cell_end - cell_begin), ptrainset))/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    // HEATFO SECTION
    if (heat_end - heat_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = heat_begin-1;
       int keyword_end_traindata = heat_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((heat_end - heat_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ((heat_end - heat_begin) - mod((heat_end - heat_begin), ptrainset))/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0;
       int reset_end = 0;
 
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;

       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((heat_end - heat_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ((heat_end - heat_begin) - mod((heat_end - heat_begin), ptrainset))/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((heat_end - heat_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ((heat_end - heat_begin) - mod((heat_end - heat_begin), ptrainset))/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    // GEOMETRY SECTION
    if (geo_end - geo_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = geo_begin-1;
       int keyword_end_traindata = geo_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((geo_end - geo_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ( (geo_end - geo_begin) - mod((geo_end - geo_begin), ptrainset) )/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0;
       int reset_end = 0;
 
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;

       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((geo_end - geo_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ( (geo_end - geo_begin) - mod((geo_end - geo_begin), ptrainset) )/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((geo_end - geo_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ( (geo_end - geo_begin) - mod((geo_end - geo_begin), ptrainset) )/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    // ENERGY SECTION
    if (energy_end - energy_begin > 1) {
       // calculate where to start and end
       int keyword_begin_traindata = energy_begin-1;
       int keyword_end_traindata = energy_end-1;
       int substart = keyword_begin_traindata + 1;
       //int subend = min(substart + int(floor((energy_end - energy_begin)/ptrainset)), keyword_end_traindata);
       int subend = min(substart + ( (energy_end - energy_begin) - mod((energy_end - energy_begin), ptrainset) )/ptrainset, keyword_end_traindata);
       int tries = 0;
       int reset_start = 0;
       int reset_end = 0;
       // print keyword
       trainset_file << traindata_nonsplit.at(keyword_begin_traindata) << endl;

       // print data for swarmcores
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
         for (int i = substart; i < subend; i++) {
            trainset_file << traindata_nonsplit.at(i) << endl;
         };
       };
       substart = subend;
       //subend = min(substart + int(floor((energy_end - energy_begin)/ptrainset)), keyword_end_traindata);
       subend = min(substart + ( (energy_end - energy_begin) - mod((energy_end - energy_begin), ptrainset) )/ptrainset, keyword_end_traindata);
       reset_start = substart;
       reset_end = subend;
       // print data for reaxffcores
       for (const int& reaxffcore : reaxffcores) {
           tries = tries + 1;
           if (core == reaxffcore) {
              if (tries == ptrainset-1) {subend = keyword_end_traindata;};
              for (int i = substart; i < subend; i++) {
                 trainset_file << traindata_nonsplit.at(i) << endl;
              };
           };
           //tries = tries + 1;
           substart = subend;
           //subend = min(substart + int(floor((energy_end - energy_begin)/ptrainset)), keyword_end_traindata);
           subend = min(substart + ( (energy_end - energy_begin) - mod((energy_end - energy_begin), ptrainset) )/ptrainset, keyword_end_traindata);
           if (tries == ptrainset-1) {
              substart = reset_start;
              subend = reset_end;
              tries = 0;
           };
       };
       // print end keyword
       trainset_file << traindata_nonsplit.at(keyword_end_traindata) << endl;
    };

    trainset_file.close();
  };
 
#endif
};

void Par::write_aveffield_lg() {
  // write distribution of parameters
  ofstream histo_file;
  // write averaged ffield
  ofstream output_file;
  // current ffield file stream
  ifstream ffield_file;
  string comment;
#ifdef WITH_MPI
if (core == 0) {
  output_file.open("ffield.ensemble", ios::out);
  ffield_file.open("ffield", ios:: in );
  histo_file.open("ffield.histo",ios::out);
};
#endif
#ifndef WITH_MPI
  output_file.open("ffield.ensemble", ios::out);
  ffield_file.open("ffield", ios:: in );
  histo_file.open("ffield.histo",ios::out);
#endif

  // store all ffield lines
  vector < string > ffield_lines;
  // store temporary lines of averaged ffield
  vector < string > averageffieldline;
  // stringstream object to convert double to formatted string
  std::stringstream ss;

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
    averageffield.push_back(ffensemble.at(0).at(m));
  };

  /* -------------------------------------
   * WRITE GENERAL PARAMS SECTION
   * -------------------------------------
   */

  histo_file << "#general parameters:\n";
  for (int m = 2; m < 41; m++) {
    boost::format f("%10.4f %s");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           ss << boost::format("%10.4f") %average;
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
      f = f % ( * it); // feeding strings into formatter
    };
    output_file << f << endl; // print the result of formatter
  };

  /* ------------------------------------
   * WRITE ATOM PARAMS HEADERS
   * ------------------------------------
   */
  for (int m = 41; m < 45; m++) {
    output_file << ffield_lines.at(m) << endl;
    averageffield.push_back(ffensemble.at(0).at(m));
  };

  using boost::is_any_of;
  // the line that contains the number of elements
  string numel_line = ffield_lines.at(41);
  // vector to store the words after split
  vector < string > results;
  boost::trim(numel_line);
  boost::split(results, numel_line, is_any_of("\t "));
  numel = stoi(results.at(0));
  max_line_atompar = 5 * numel + 45;

  /* -------------------------------------
   * WRITE ATOM PARAMS SECTION
   * -------------------------------------
   */

  histo_file << "#atom parameters:\n";
  int m = 45;
  while (m < max_line_atompar) {
    boost::format f("% s%   9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions( f.exceptions() &
    ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           ss << boost::format("%9.4f") %average;
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
           averageffieldline.push_back( ss.str() );
           ss.str("");
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector<std::string>::iterator it=averageffield.at(m).begin();it!=averageffield.at(m).end();++it){
            f = f % (*it);
    };
    output_file << f << endl;

    for (int i = 1; i < 4; i++){
      m = m + 1;
      boost::format f("%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions( f.exceptions() &
      ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );

      double average = 0;
      ss << boost::format("   ");
      for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
          if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
             for (int k = 0; k < ffensemble.size(); k++) {
                 histo_file << ffensemble.at(k).at(m).at(col) << " ";
                 average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
             };
             ss << boost::format("%9.4f") %average;
             averageffieldline.push_back(ss.str());
             average = 0;
             ss.str("");
          }else{
             ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
             averageffieldline.push_back( ss.str() );
             ss.str("");
          };
          histo_file << "\n";
      };
      averageffield.push_back(averageffieldline);
      averageffieldline.clear();

      for (std::vector<std::string>::iterator it=averageffield.at(m).begin();it!=averageffield.at(m).end();++it){
              f = f % (*it);
      };
      output_file << f << endl;
    };

    m = m + 1;
    boost::format f2("   %9.4f%9.4f");
    f2.exceptions( f2.exceptions() &
    ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );

    average = 0;
    //ss << boost::format("   ");
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           ss << boost::format("%9.4f") %average;
           averageffieldline.push_back(ss.str());
           average = 0;
           ss.str("");
        }else{
           ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
           averageffieldline.push_back( ss.str() );
           ss.str("");
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector<std::string>::iterator it=averageffield.at(m).begin();it!=averageffield.at(m).end();++it){
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
    averageffield.push_back(ffensemble.at(0).at(m));
 };

  using boost::is_any_of;
  // the line the contains the number of bond types
  string numbty_line = ffield_lines.at(max_line_atompar);
  // vector to store the words after split
  vector < string > results_bonds;
  boost::trim(numbty_line);
  boost::split(results_bonds, numbty_line, is_any_of("\t "));
  numbty = stoi(results_bonds.at(0));
  max_line_bondpar = 2 * numbty + max_line_atompar;

  /* ------------------------------------
   * WRITE BONDS PARAMS SECTION
   * ------------------------------------
   */

  histo_file << "#bonds parameters:\n";
  for (int m = max_line_atompar + 2; m < max_line_bondpar + 2; m++) {
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
           averageffieldline.push_back( ss.str() );
           ss.str("");
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();


    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
      f = f % ( * it);
    };

    output_file << f << endl;
    m = m + 1;
    boost::format f2("      %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f2.exceptions(f2.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));


    average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           ss << boost::format("%9.4f") %average;
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           ss << boost::format("%i") %ffensemble.at(0).at(m).at(col);
           averageffieldline.push_back( ss.str() );
           ss.str("");

        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();


    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numodty = stoi(results_offdiag.at(0));
  max_line_offdpar = max_line_bondpar + 3 + numodty;
  
  histo_file << "#off-diag parameters:\n";
  for (int m = max_line_bondpar + 3; m < max_line_offdpar; m++) {
    // the last entry is 10.4f because dispersion coeff. can get to 4 digits long.
    // So, to prevent it sticking to the left column
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           if (col == ffensemble.at(0).at(m).size()-1) {
              ss << boost::format("%10.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numaty = stoi(results_angle.at(0));
  max_line_angles = max_line_offdpar + 1 + numaty;

  histo_file << "#angle parameters:\n";
  for (int m = max_line_offdpar + 1; m < max_line_angles; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));


    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1) || (col == 2)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();


    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numtoty = stoi(results_tors.at(0));
  max_line_tors = max_line_angles + 1 + numtoty;

  histo_file << "#torsion parameters:\n";
  for (int m = max_line_angles + 1; m < max_line_tors; m++) {
    boost::format f("  %i  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1) || (col == 2) || (col == 3)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numhbty = stoi(results_hb.at(0));
  max_line_hbs = max_line_tors + 1 + numhbty;

  histo_file << "#hbond parameters:\n";
  for (int m = max_line_tors + 1; m < max_line_hbs; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1) || (col == 2)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
    histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  output_file.close();
  histo_file.close();

};

void Par::write_aveffield() {
  // distribution of parameters
  ofstream histo_file;
  // write averaged ffield
  ofstream output_file;
  // current ffield file stream
  ifstream ffield_file;
  string comment;
#ifdef WITH_MPI
if (core == 0) {
  output_file.open("ffield.ensemble", ios::out);
  ffield_file.open("ffield", ios:: in );
  histo_file.open("ffield.histo",ios::out);
};
#endif
#ifndef WITH_MPI
  output_file.open("ffield.ensemble", ios::out);
  ffield_file.open("ffield", ios:: in );
  histo_file.open("ffield.histo",ios::out);
#endif

  // store all ffield lines
  vector < string > ffield_lines;
  // store temporary lines of averaged ffield
  vector < string > averageffieldline;
  // stringstream object to convert double to formatted string
  std::stringstream ss;

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
    averageffield.push_back(ffensemble.at(0).at(m));
  };

  /* -------------------------------------
   * WRITE GENERAL PARAMS SECTION
   * -------------------------------------
   */

  histo_file << "#general parameters:\n";
  for (int m = 2; m < 41; m++) {
    boost::format f("%10.4f %s");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           ss << boost::format("%10.4f") %average;
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
      f = f % ( * it); // feeding strings into formatter
    };
    output_file << f << endl; // print the result of formatter
  };
  /* ------------------------------------
   * WRITE ATOM PARAMS HEADERS
   * ------------------------------------
   */
  for (int m = 41; m < 45; m++) {
    output_file << ffield_lines.at(m) << endl;
    averageffield.push_back(ffensemble.at(0).at(m));
  };

  using boost::is_any_of;
  // the line that contains the number of elements
  string numel_line = ffield_lines.at(41);
  // vector to store the words after split
  vector < string > results;
  boost::trim(numel_line);
  boost::split(results, numel_line, is_any_of("\t "));
  numel = stoi(results.at(0));
  max_line_atompar = 4 * numel + 45;

  /* -------------------------------------
   * WRITE ATOM PARAMS SECTION
   * -------------------------------------
   */

  histo_file << "#atom parameters:\n";
  int m = 45;
  while (m < max_line_atompar) {
    boost::format f("% s%   9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions( f.exceptions() &
    ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           ss << boost::format("%9.4f") %average;
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
           averageffieldline.push_back( ss.str() );
           ss.str("");
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector<std::string>::iterator it=averageffield.at(m).begin();it!=averageffield.at(m).end();++it){
            f = f % (*it);
    };
    output_file << f << endl;

    for (int i = 1; i < 4; i++){
      m = m + 1;
      boost::format f("%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions( f.exceptions() &
      ~ ( boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );

      double average = 0;
      ss << boost::format("   ");
      for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
          if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
             for (int k = 0; k < ffensemble.size(); k++) {
                 histo_file << ffensemble.at(k).at(m).at(col) << " ";
                 average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
             };
             ss << boost::format("%9.4f") %average;
             averageffieldline.push_back(ss.str());
             average = 0;
             ss.str("");
          }else{
             ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
             averageffieldline.push_back( ss.str() );
             ss.str("");
          };
          histo_file << "\n";
      };
      averageffield.push_back(averageffieldline);
      averageffieldline.clear();

      for (std::vector<std::string>::iterator it=averageffield.at(m).begin();it!=averageffield.at(m).end();++it){
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
    averageffield.push_back(ffensemble.at(0).at(m));
 };

  using boost::is_any_of;
  // the line the contains the number of bond types
  string numbty_line = ffield_lines.at(max_line_atompar);
  // vector to store the words after split
  vector < string > results_bonds;
  boost::trim(numbty_line);
  boost::split(results_bonds, numbty_line, is_any_of("\t "));
  numbty = stoi(results_bonds.at(0));
  max_line_bondpar = 2 * numbty + max_line_atompar;

  /* ------------------------------------
   * WRITE BONDS PARAMS SECTION
   * ------------------------------------
   */

  histo_file << "#bonds parameters:\n";
  for (int m = max_line_atompar + 2; m < max_line_bondpar + 2; m++) {
    if (mod(m, 2) != 0) {
      boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

      double average = 0;
      for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
          if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
             for (int k = 0; k < ffensemble.size(); k++) {
                 histo_file << ffensemble.at(k).at(m).at(col) << " ";
                 average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
             };
             if ((col == 0) || (col == 1)) {
                ss << boost::format("%i") %average;
             }else{
                ss << boost::format("%9.4f") %average;
             };
             averageffieldline.push_back(ss.str());
             ss.str("");
             average = 0;
          }else{
             ss << boost::format("%s") %ffensemble.at(0).at(m).at(col);
             averageffieldline.push_back( ss.str() );
             ss.str("");
          };
          histo_file << "\n";
      };
      averageffield.push_back(averageffieldline);
      averageffieldline.clear();

      for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
        f = f % ( * it);
      };
      output_file << f << endl;
    } else {
      boost::format f("      %9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
      f.exceptions(f.exceptions() &
        ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

      double average = 0;
      for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
          if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
             for (int k = 0; k < ffensemble.size(); k++) {
                 histo_file << ffensemble.at(k).at(m).at(col) << " ";
                 average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
             };
             ss << boost::format("%9.4f") %average;
             averageffieldline.push_back( ss.str() );
             ss.str("");
             average = 0;
          }else{
             ss << boost::format("%i") %ffensemble.at(0).at(m).at(col);
             averageffieldline.push_back( ss.str() );
             ss.str("");
             
          };
          histo_file << "\n";
      };
      averageffield.push_back(averageffieldline);
      averageffieldline.clear();

      for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numodty = stoi(results_offdiag.at(0));
  max_line_offdpar = max_line_bondpar + 3 + numodty;

  histo_file << "#off-diag parameters:\n";
  for (int m = max_line_bondpar + 3; m < max_line_offdpar; m++) {
    // the last entry is 10.4f because dispersion coeff. can get to 4 digits long.
    // So, to prevent it sticking to the left column
    boost::format f("  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numaty = stoi(results_angle.at(0));
  max_line_angles = max_line_offdpar + 1 + numaty;

  histo_file << "#angle parameters:\n";
  for (int m = max_line_offdpar + 1; m < max_line_angles; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));


    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1) || (col == 2)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();


    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numtoty = stoi(results_tors.at(0));
  max_line_tors = max_line_angles + 1 + numtoty;

  histo_file << "#torsion parameters:\n";
  for (int m = max_line_angles + 1; m < max_line_tors; m++) {
    boost::format f("  %i  %i  %i  %i%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1) || (col == 2) || (col == 3)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
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
    averageffield.push_back(ffensemble.at(0).at(m));
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
  numhbty = stoi(results_hb.at(0));
  max_line_hbs = max_line_tors + 1 + numhbty;

  histo_file << "#hbond parameters:\n";
  for (int m = max_line_tors + 1; m < max_line_hbs; m++) {
    boost::format f("  %i  %i  %i%9.4f%9.4f%9.4f%9.4f");
    f.exceptions(f.exceptions() &
      ~(boost::io::too_many_args_bit | boost::io::too_few_args_bit));

    double average = 0;
    for (int col = 0; col < ffensemble.at(0).at(m).size(); col++) {
        if ( is_numeric(ffensemble.at(0).at(m).at(col)) ) {
           for (int k = 0; k < ffensemble.size(); k++) {
               histo_file << ffensemble.at(k).at(m).at(col) << " ";
               average += stod(ffensemble.at(k).at(m).at(col))/numffbags;
           };
           if ((col == 0) || (col == 1) || (col == 2)) {
              ss << boost::format("%i") %average;
           }else{
              ss << boost::format("%9.4f") %average;
           };
           averageffieldline.push_back(ss.str());
           ss.str("");
           average = 0;
        }else{
           averageffieldline.push_back( ffensemble.at(0).at(m).at(col) );
        };
        histo_file << "\n";
    };
    averageffield.push_back(averageffieldline);
    averageffieldline.clear();

    for (std::vector < std::string > ::iterator it = averageffield.at(m).begin(); it != averageffield.at(m).end(); ++it) {
      f = f % ( * it);
    };
    output_file << f << endl;
  };

  output_file.close();
  histo_file.close();
};

void Par::write_ffield(const vector <double> &active_params, int cycle, int iter, int parid) {
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(parid);
  // update particle ffield
  int index = 0;
  double tempphys;
  // for each line value in ffline and corresponding column value
  // in ffcol cp corresponding pos value into ffieldmat.
  for (int line: ffline) {
    // buffer for sprintf
    char buffer[50];
    double n;
    // transform back to physical positions before writing ffield
    //tempphys = pos.at(index)*(maxdomain.at(index) - mindomain.at(index)) + mindomain.at(index);
    tempphys = active_params.at(index)*(maxdomain.at(index) - mindomain.at(index)) + mindomain.at(index);
    // use physical positions
    n = sprintf(buffer, "%9.4f", tempphys);
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
  numel = stoi(results.at(0));
  max_line_atompar = 4 * numel + 45;
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
  numbty = stoi(results_bonds.at(0));
  max_line_bondpar = 2 * numbty + max_line_atompar;

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
  numodty = stoi(results_offdiag.at(0));
  max_line_offdpar = max_line_bondpar + 3 + numodty;

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
  numaty = stoi(results_angle.at(0));
  max_line_angles = max_line_offdpar + 1 + numaty;

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
  numtoty = stoi(results_tors.at(0));
  max_line_tors = max_line_angles + 1 + numtoty;

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
  numhbty = stoi(results_hb.at(0));
  max_line_hbs = max_line_tors + 1 + numhbty;

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
  //boost::filesystem::path pwd(boost::filesystem::current_path());
#ifdef WITH_MPI
  boost::filesystem::copy_file("CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    "CPU." + str_core + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
#endif
#ifndef WITH_MPI
  boost::filesystem::copy_file("ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    "ffield", boost::filesystem::copy_option::overwrite_if_exists);
#endif
};


void Par::write_ffield_lg(const vector <double> &active_params, int cycle, int iter, int parid) {
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(parid);

  // update particle ffield
  int index = 0;
  double tempphys;
  // for each line value in ffline and corresponding column value in ffcol 
  // cp corresponding pos value into ffieldmat. Range-based for statement 
  // following: https://msdn.microsoft.com/en-us/library/jj203382.aspx
  for (int line: ffline) {
    // buffer for sprintf
    char buffer[50];
    double n;
    // transform back to physical positions before writing ffield
    //tempphys = pos.at(index)*(maxdomain.at(index) - mindomain.at(index)) + mindomain.at(index);
    tempphys = active_params.at(index)*(maxdomain.at(index) - mindomain.at(index)) + mindomain.at(index);
    // use physical positions
    n = sprintf(buffer, "%9.4f", tempphys);
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
  numel = stoi(results.at(0));
  max_line_atompar = 5 * numel + 45;

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
  numbty = stoi(results_bonds.at(0));
  max_line_bondpar = 2 * numbty + max_line_atompar;

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
  numodty = stoi(results_offdiag.at(0));
  max_line_offdpar = max_line_bondpar + 3 + numodty;

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
  numaty = stoi(results_angle.at(0));
  max_line_angles = max_line_offdpar + 1 + numaty;

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
  numtoty = stoi(results_tors.at(0));
  max_line_tors = max_line_angles + 1 + numtoty;

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
  numhbty = stoi(results_hb.at(0));
  max_line_hbs = max_line_tors + 1 + numhbty;

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
  //boost::filesystem::path pwd(boost::filesystem::current_path());
#ifdef WITH_MPI
  boost::filesystem::copy_file("CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    "CPU." + str_core + "/ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp.*");
#endif
#ifndef WITH_MPI
  boost::filesystem::copy_file("ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
    "ffield", boost::filesystem::copy_option::overwrite_if_exists);
  // remove temporary ffield so other particles do not append to the file
  //boost::filesystem::remove(pwd.string() + "/ffield.tmp.*");
#endif

};

void Par::check_bounds_contff() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered check_bounds_contff()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered check_bounds_contff()" << endl;
};
#endif

  // check physical bounds (since params in ffieldmat are always physical)
  for (int i = 0; i < dim; i++) { 
      if (pos.at(i) > maxdomain.at(i)) {
          pos.at(i) = maxdomain.at(i) - 0.0001;
      };
      if (pos.at(i) < mindomain.at(i)) {
          pos.at(i) = mindomain.at(i) + 0.0001;
      };
      if (pos.at(i) == 0.0) {
          pos.at(i) = 0.0001;
      };
  };                                         

};

void Par::read_bounds() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered read_bounds()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered read_bounds()" << endl;
};
#endif

  // Reading params columns into allData array split by entries
  // following: https://stackoverflow.com/questions/10521658/reading-files-columns-into-array

  std::vector < std::vector < double > > allData;
  string str_core = std::to_string(core);
  std::ifstream fin("params.mod");
  if (fin.fail()) {
#ifdef WITH_MPI
    cout << "Error: Unable to open parameters file 'params.mod' on CPU " << core << ". \n";
    fin.close();
    MPI_Abort(MPI_COMM_WORLD,2);
#endif
#ifndef WITH_MPI
    cout << "Error: Unable to open parameters file 'params.mod' on CPU " << core << ". \n";
    fin.close();
    exit(EXIT_FAILURE);
#endif
  } else {
    std::string line;
    int numlines = 0;

    // for each line
    while (std::getline(fin, line)) {
      if (line[0] == '#') {continue;};
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
        maxdomain.at(i) = (1.0 + perc)*maxdomain.at(i);
        // for negative values exchange bounds
        if (maxdomain.at(i) < mindomain.at(i)) {
           double temp;
           temp = maxdomain.at(i);
           maxdomain.at(i) = mindomain.at(i);
           mindomain.at(i) = temp;
        };
      };

    };
if (core == 0 && verbose == true) {
    cout << "parameters to train with their bounds:" << endl;
    for (int i = 0; i < dim; i++) {
        cout << ffline.at(i)+1 << " " << ffcol.at(i)+1 << " " << mindomain.at(i) << " " << maxdomain.at(i) << " " << endl;
    };
};
  };
};

double Par::get_min_dim() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_min_dim()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_min_dim()" << endl;
};
#endif

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
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_normdir()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_normdir()" << endl;
};
#endif

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
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_levy_McCul()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_levy_McCul()" << endl;
};
#endif

  // ****************************************************** //
  // Gaussian mutation (no Levy)
   std::normal_distribution < double > distnorm(0.0, 1.0);
   double x = distnorm(generator); // + tau;
  // ****************************************************** //

  // ****************************************************** //
  // Generation of Levy symmetric distribution (beta=0) using
  // McCulloch's Algorithm
  // When alpha approaches 2, the distribution becomes Gaussian with mean tau
  // and variance 2 * c ^ 2, and beta has no effect.

  //std::uniform_real_distribution <double> dist1(0.0,1.0);
  //// generate variable from exponential distribution
  //double w = -log(dist1(generator));
  //// generate variable from uniform distribution
  //double phi = (dist1(generator) - 0.5)*pi;
  ////double alpha_init = 1.0;
  ////double alpha_fin = 2.0;
  //double alpha = 1.5; // alpha_init + iter * (alpha_fin - alpha_init) / maxiters; //min(2.0,exp(log(2.0)/10000.0*time)); // alpha_init + time * (alpha_fin - alpha_init) / maxiters;
  //double c = 1.0; // 0.01 - (0.01 - 1e-50)*iter / maxiters; // 0.0003*l2_norm(space_range) - 0.0066; // 0.0003*l2_norm(space_range) - 0.0069; // 0.1 - (0.1 - 0.001)*(iter / maxiters);														// scaling parameter
  ////cout << 0.0003*l2_norm(space_range) - 0.0066  << endl;
  //double tau = 0.0;														// location parameter
  //double x = c*pow((cos((1.0 - alpha)*phi)) / w, 1.0 / alpha - 1.0)*(sin(alpha*phi)) / pow(cos(phi), 1.0 / alpha) + tau;
  ////cout << 0.0003*l2_norm(space_range) << endl;
  // ***************************************************** //
  return x;
};

void Par::update_bpos() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered update_bpos()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered update_bpos()" << endl;
};
#endif

  for (int j = 0; j < dim; j++) {
    bpos.at(j) = pos.at(j);
  };
};

void Par::update_vel(double inertiaf, double CF, vector < double > globpos, double iter) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered update_vel()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered update_vel()" << endl;
};
#endif

  double r1 = dist1(generator);
  double r2 = dist1(generator);

  for (int v = 0; v < dim; v++) {

    // velocity update with perturbations
    vel.at(v) = CF * (inertiaf * vel.at(v) + c1 * r1 * (bpos.at(v) - pos.at(v)) + c2 * r2 * (globpos.at(v) - pos.at(v)));

  };
};

void Par::update_pos() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered update_pos()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered update_pos()" << endl;
};
#endif

  for (int i = 0; i < dim; i++) {
    pos.at(i) = pos.at(i) + vel.at(i);
    //check physical bounds
    //std::uniform_real_distribution < double > dist2(mindomain.at(i), maxdomain.at(i));
    //if (pos.at(i) > maxdomain.at(i)) {
    //  pos.at(i) = dist2(generator);
    //};
    //if (pos.at(i) < mindomain.at(i)) {
    //  pos.at(i) = dist2(generator);
    //};
    // check standardized bounds
    std::uniform_real_distribution < double > dist2(0.0001,1.0);
    if (pos.at(i) > 1.0) {
      pos.at(i) = dist2(generator);
    };
    if (pos.at(i) < 0.0001) {
      pos.at(i) = dist2(generator);
    };
  };
};

void Par::update_pos_levy(vector < double > globpos, double iter, double inertiaf) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered update_pos_levy()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered update_pos_levy()" << endl;
};
#endif

  //std::normal_distribution <double> normdist2(0.0, 1.0);
  double levystep = fabs(get_levy_McCul(iter, maxiters));
  vector < double > direction = get_normdir();

  for (int i = 0; i < dim; i++) {
    // physical positions
    //std::uniform_real_distribution < double > dist2(mindomain.at(i), maxdomain.at(i));
    // standardized positions
    std::uniform_real_distribution < double > dist2(0.0001,1.0);
    // levy step on physical position
    //pos.at(i) = pos.at(i) + levyscale * get_min_dim() * levystep * direction.at(i);
    // old step
    //pos.at(i) = pos.at(i) + levyscale*(maxdomain.at(i) - mindomain.at(i))*levystep*direction.at(i);
    // levy step on standardized position
    pos.at(i) = pos.at(i) + levyscale * levystep * direction.at(i);

    // check bounds on physical positions
    //if (pos.at(i) > maxdomain.at(i)) {
    //  pos.at(i) = dist2(generator);
    //};
    //if (pos.at(i) < mindomain.at(i)) {
    //  pos.at(i) = dist2(generator);
    //};
    // check bounds on standardized positions
    if (pos.at(i) > 1.0) {
      pos.at(i) = dist2(generator);
    }; 
    if (pos.at(i) < 0.0001) {
      pos.at(i) = dist2(generator);
    };

  };
};


void Par::iterate() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered iterate()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered iterate()" << endl;
};
#endif

 std::vector <double> standard_mindomain(dim, 0.0);
 std::vector <double> standard_maxdomain(dim, 1.0);

 double temphys;

  // if we parallelize the training set, each swarmcore sets the positions of its reaxffcores
  if (ptrainset > 1) {
     for (const int& swarmcore : swarmcores) {
           for (int i = 1; i < reaxffcores.size()/swarmcores.size() + 1; i++) {
             int reaxffcore = swarmcore + i;
             if (core == swarmcore ) {
               MPI_Send( pos.data(), pos.size(), MPI_DOUBLE, reaxffcore, 11, MPI_COMM_WORLD );
               //cout << "swarmcore " << swarmcore << " sending pos[0] = " << active_params.at(0) << " pos[1] = " << active_params.at(1) << endl;
             };
             if (core == reaxffcore) {
               MPI_Recv( pos.data(), pos.size(), MPI_DOUBLE, swarmcore, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
               //cout << "reaxffcore " << reaxffcore << " receiving pos[0] = " << pos.at(0) << " pos[1] = " << pos.at(1) << endl;
             };
         };
     };
  };

 // note: for localmin, dropout means ranges that belong to dropped dimensions, are set equal
 // so to exclude them from local minimization
 if (regular == 3) {
    for (int j : dropped_dimns) {
        mindomain.at(j) = 0.0001;
        maxdomain.at(j) = 0.0001;
        pos.at(j) = 0.0001;
    };
 };


 if (localmin == 2) {
   nlopt::opt opt(nlopt::LN_SBPLX, dim);
   opt.set_min_objective(fitness_wrapper, this);
   // use physical bounds
   //opt.set_lower_bounds(mindomain);
   //opt.set_upper_bounds(maxdomain);
   // use standardized bounds
   opt.set_lower_bounds(standard_mindomain);
   opt.set_upper_bounds(standard_maxdomain);
   // note: relative tolerance should be used rather than absolute
   // so to take into account different magnitudes of different parameters
   opt.set_xtol_rel(lm_err_tol);
   opt.set_maxeval(lm_iter_max);
   x = pos;
   //reset counter for swarm core lmin successes
   count_lmin_success=0;

   try{
       nlopt::result result = opt.optimize(x, minf);
       if (ptrainset > 1) {
          if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
               lmin_rep << "local min completed! --> flocky iter: " << state.iter << endl;
               lmin_rep << "new local min: " << boost::format("        %1.10e") %minf << endl;
               lmin_rep << "new parameters:" << endl;
               lmin_rep << "[ ";
               // mark lmin as success and increase overall counter across active swarm members
               // if count_lmin_success == swarmcores.size(), this will cause a forced-termination as given in the wrapper.
               lmin_success=1;
               //MPI_Allreduce( & lmin_success, & count_lmin_success, 1, MPI_INT, MPI_SUM, ACTIVESWARM);

               for (int m = 0; m < dim; m++) {
                 // print physical local min positions
                 temphys = x[m]*(maxdomain.at(m) - mindomain.at(m)) + mindomain.at(m);
                 lmin_rep << boost::format("    %8.4f") %temphys << " ";
               };
               lmin_rep << "]\n" << endl;
               lmin_rep.close();
          };
          // sum the success counter throughout MPI_COMM_WORLD. Passive swarm members will contribute 0 each,
          // while the active swarm members will contribute 1 each.
          MPI_Allreduce( & lmin_success, & count_lmin_success, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
       }else {
               #ifdef WITH_MPI
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
               #endif
               #ifndef WITH_MPI
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
               #endif
               lmin_rep << "local min completed! --> flocky iter: " << state.iter << endl;
               lmin_rep << "new local min: " << boost::format("        %1.10e") %minf << endl;
               lmin_rep << "new parameters:" << endl;
               lmin_rep << "[ ";
               for (int m = 0; m < dim; m++) {
                 // print physical local min positions
                 temphys = x[m]*(maxdomain.at(m) - mindomain.at(m)) + mindomain.at(m);
                 lmin_rep << boost::format("%8.4f") %temphys << " ";
               };
               lmin_rep << "]\n" << endl;
               lmin_rep.close();
       };
   }
   catch(std::exception &e) {
       if (ptrainset > 1) {
          if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              #ifdef WITH_MPI
              ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
              #endif
              #ifndef WITH_MPI
              ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
              #endif
              lmin_rep << "Warning: local minimization failed for CPU " << core << "--> " << e.what() << endl;
              lmin_rep.close();
          };
       }else { 
           #ifdef WITH_MPI
           ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
           #endif
           #ifndef WITH_MPI
           ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
           #endif
           lmin_rep << "Warning: local minimization failed for CPU " << core << "--> " << e.what() << endl;
           lmin_rep.close();
       };
   };
 };
 

 if (localmin == 1) {
   nlopt::opt opt(nlopt::LD_LBFGS, dim);
   opt.set_min_objective(fitness_wrapper, this);  
   opt.set_vector_storage(10);
   // use physical bounds
   //opt.set_lower_bounds(mindomain);
   //opt.set_upper_bounds(maxdomain);
   // use standardized bounds
   opt.set_lower_bounds(standard_mindomain);
   opt.set_upper_bounds(standard_maxdomain);
   // note: relative tolerance should be used rather than absolute
   // so to take into account different magnitudes of different parameters
   opt.set_xtol_rel(lm_err_tol);
   opt.set_maxeval(lm_iter_max);
   x = pos;

   try{
       nlopt::result result = opt.optimize(x, minf);
       if (ptrainset > 1) {
          if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
               #ifdef WITH_MPI
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
               #endif
               #ifndef WITH_MPI
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
               #endif
               lmin_rep << "local min completed! --> flocky iter: " << state.iter << endl;
               lmin_rep << "new local min: " << boost::format("        %1.10e") %minf << endl;
               lmin_rep << "new parameters:" << endl;
               lmin_rep << "[ ";
               for (int m = 0; m < dim; m++) {
                 // print physical local min positions
                 temphys = x[m]*(maxdomain.at(m) - mindomain.at(m)) + mindomain.at(m);
                 lmin_rep << boost::format("%8.4f") %temphys << " ";
               };
               lmin_rep << "]\n" << endl;
               lmin_rep.close();
          };
       }else {
               #ifdef WITH_MPI
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
               #endif
               #ifndef WITH_MPI
               ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
               #endif
               lmin_rep << "local min completed! --> flocky iter: " << state.iter << endl;
               lmin_rep << "new local min: " << boost::format("        %1.10e") %minf << endl;
               lmin_rep << "new parameters:" << endl;
               lmin_rep << "[ ";
               for (int m = 0; m < dim; m++) {
                 // print physical local min positions
                 temphys = x[m]*(maxdomain.at(m) - mindomain.at(m)) + mindomain.at(m);
                 lmin_rep << boost::format("%8.4f") %temphys << " ";
               };
               lmin_rep << "]\n" << endl;
               lmin_rep.close();
       };
   }
   catch(std::exception &e) {
       if (ptrainset > 1) {
          if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              #ifdef WITH_MPI
              ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
              #endif
              #ifndef WITH_MPI
              ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
              #endif
              lmin_rep << "Warning: local minimization failed for CPU " << core << "--> " << e.what() << endl;
              lmin_rep.close();
          };
       }else { 
           #ifdef WITH_MPI
           ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle)+"."+std::to_string(core), ofstream::app);
           #endif
           #ifndef WITH_MPI
           ofstream lmin_rep("lmin_rep.out." + std::to_string(state.cycle), ofstream::app);
           #endif
           lmin_rep << "Warning: local minimization failed for CPU " << core << "--> " << e.what() << endl;
           lmin_rep.close();
       };      
   };
 };

};

double Par::eval_fitness(const vector <double> &active_params, void *my_func_data) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered eval_fitness()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered eval_fitness()" << endl;
};
#endif

#ifdef WITH_MPI
  int cycle;
  int iter;
  int parid;

  long double evalfit;
  long double pfitness;
  double evalfitdouble;

  cycle = state.cycle;
  iter = state.iter;
  parid = state.parid;

  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_core = std::to_string(core);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  // count total # func evaluations
  funceval = funceval + 1;

  // prepare fort.3 files
  boost::filesystem::copy_file("CPU." + str_core + "/geo." + str_cycle + "." + str_parID,
     "CPU." + str_core + "/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

  // if we parallelize the training set, each swarmcore sets the positions of its reaxffcores
  if (ptrainset > 1) {
     for (const int& swarmcore : swarmcores) {
           for (int i = 1; i < reaxffcores.size()/swarmcores.size() + 1; i++) {
             int reaxffcore = swarmcore + i;
             if (core == swarmcore ) {
               pos = active_params;
               MPI_Send( active_params.data(), active_params.size(), MPI_DOUBLE, reaxffcore, 0, MPI_COMM_WORLD );
               //cout << "swarmcore " << swarmcore << " sending pos[0] = " << boost::format("%18.4f") %active_params.at(0) << " pos[1] = " << boost::format("%18.4f") %active_params.at(1) << endl;
             };
             if (core == reaxffcore) {
               MPI_Recv( pos.data(), pos.size(), MPI_DOUBLE, swarmcore, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
               //cout << "reaxffcore " << reaxffcore << " receiving pos[0] = " << boost::format("%18.4f") %pos.at(0) << " pos[1] = " << boost::format("%18.4f") %pos.at(1) << endl;
             };
         };
     };
  };

  // check if ffield is LG or not. execute correct tapreaxff accordingly
  if (lg_yn == true) {
     if (ptrainset > 1 && find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
        write_ffield_lg(pos, cycle, iter, parid);
     }else{
         write_ffield_lg(active_params, cycle, iter, parid);
     };
  } else {
      if (ptrainset > 1 && find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
         write_ffield(pos, cycle, iter, parid);
      }else{
          write_ffield(active_params, cycle, iter, parid);
      };
  };

  // check if tapreaxff exec is present
  std::ifstream fin5(("CPU." + str_core + "/tapreaxff").c_str());
  if (fin5.fail()) {
     cout << "tapreaxff executable not found for CPU " + str_core + ". Aborting! \n";
     fin5.close();
     MPI_Abort(MPI_COMM_WORLD,3);
  };
  fin5.close();
  boost::filesystem::copy_file("CPU." + str_core + "/ffield",
      "CPU." + str_core + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  if (fixcharges == true){
     boost::filesystem::copy_file("CPU." + str_core + "/charges",
       "CPU." + str_core + "/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
  };
  // cd to each CPU.x directory
  string old_path = pwd.string();
  boost::filesystem::path p("CPU." + str_core);
  boost::filesystem::current_path(p);
  //arguments for tapreaxff, will run: tapreaxff                                                                                                                                  
  char *args[3] = { "./tapreaxff", "", NULL} ;
  pid_t c_pid, pid;
  int status;

  /* create a child process */
  c_pid = fork();

  if (c_pid == 0){
    /* CHILD */
    char *filename  = "run.log";
    int outfile = open(filename, O_CREAT | O_WRONLY, S_IRWXU);
    if (outfile == -1){
          fprintf(stderr, "Error: failed to create file %s\n", filename);
    }else{
          /* redirect the standard output from this process to the file. */
          if(dup2(outfile, STDOUT_FILENO) != STDOUT_FILENO){
            fprintf(stderr, "Error: failed to redirect standard output\n");
          }
          /* redirect tapreaxff stdout to outfile*/
          dup2 (outfile, STDOUT_FILENO);
          /* redirect tapreaxff stderr to /dev/null */
          dup2(open("/dev/null", 0), 2);

        /* printf("Child: executing args\n"); */
        // execute args                                                                                                                                                               
        execvp( args[0], args);
        // only get here if exec failed                                                                                                                                             
        perror("execve failed");
        wait(&status);
    };
  }else if (c_pid > 0){
    /* PARENT */
    if( (pid = wait(&status)) < 0){
      perror("wait");
      _exit(1);
    };
    //printf("Parent: finished\n");
  }else{
    perror("fork failed");
    _exit(1);
  };
  if (WIFSIGNALED (status)) {
    cout << "tapreaxff exited abnormaly on CPU:" << core << "\n";
    MPI_Abort(MPI_COMM_WORLD,4);
  };

  // cd back to main directory
  boost::filesystem::path p2(old_path);
  boost::filesystem::current_path(p2);
  // read fitness value contained in fort.13 file
  string str;
  if ( !boost::filesystem::exists( "CPU." + str_core + "/fort.13" ) )
  {
    evalfit = numeric_limits < double > ::infinity();
    if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
         pfitness = evalfit;
    };
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
    if (str.at(0) == '*' || !isdigit(str.at(0))) {
      evalfit = numeric_limits < double > ::infinity();
      if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
           pfitness = evalfit;
      };
    } else {
      // convert to double
      if (regular == 1 || regular == 2 ) {
        evalfit = stod(str) + get_reg();
        if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
           pfitness = evalfit;
        };
      }else{
        evalfit = stod(str);
        if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
           pfitness = evalfit;
        };
      };
      file13.close();
    };
    boost::filesystem::remove( "CPU." + str_core + "/fort.13" );
  };
  // If we parallelize the training set, each reaxffcore adds its fitness to the fitness of its swarmcore
  // Note: The following MPI_Send/Recv logic of assigning swarmcore and reaxffcore ranks is equivalent to
  // the logic implemented in the top of this function
  // that exploits the fact that number of cores in reaxffcores that belong to a swarmcore, is half the 
  // total number of reaxffcores. 
  if (ptrainset > 1) {
     int j = 1;
     int swarmcore;

     for (int reaxffcore : reaxffcores) {
         swarmcore = reaxffcore - j;
         if (core == reaxffcore ) {
           MPI_Send( &pfitness, 1, MPI_LONG_DOUBLE, swarmcore, 1, MPI_COMM_WORLD );
           //cout << "reaxffcore " << reaxffcore << " sent " << boost::format("%1.10e") %pfitness << " to swarmcore " << swarmcore << endl;
         };
         if (core == swarmcore) {
           MPI_Recv( &pfitness, 1, MPI_LONG_DOUBLE, reaxffcore, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
           //cout << "swarmcore " << swarmcore << "received " << boost::format("%1.10e") %pfitness << "from reaxffcore " << reaxffcore << endl;
         };
         if (core == swarmcore) {
            evalfit = evalfit + pfitness;
            //cout << "swarmcore " << swarmcore << " newfitness: " << boost::format("%1.10e") %evalfit << endl;
         };

         if (j == ptrainset-1) {
           j = 1;
         }else{
           j = j + 1;
         };
     };
  };

  if (ptrainset > 1) {
    // swarmcores setting reaxffcores to trick localmin
    for (const int& swarmcore : swarmcores) {
          for (int i = 1; i < reaxffcores.size()/swarmcores.size() + 1; i++) {
            int reaxffcore = swarmcore + i;
            if (core == swarmcore ) {
              MPI_Send( &evalfit, 1, MPI_LONG_DOUBLE, reaxffcore, 3, MPI_COMM_WORLD );
              //if (verbose == true) {
              //        cout << ">>iter = " << iter << " CPU: " << core << " swarmcore sending fitness to CPU: " << core << " reaxffcore " << endl;
              //    cout << endl;
              //};
            };
            if (core == reaxffcore) {
              MPI_Recv( &evalfit, 1, MPI_LONG_DOUBLE, swarmcore, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
              //if (verbose == true) {
              //       cout << ">>iter = " << iter << " CPU: " << core << " reaxffcore receiving fitness = " << evalfit << " from CPU: " << core << " swarmcore " << endl;
              //   cout << endl;
              //};
            };
        };
    };
  };


//MPI_Barrier(MPI_COMM_WORLD);
  // Note: do not update the geometry file during iterations. Each member should use one geo file throughout
  // the training. Assuming we start with DFT_optimized (or sensible structures), and that we use some small
  // number of structural minimizations (3-10), the sensible structures won't change much, so no need to
  // provide previous structure as initial structure for the next round of iteration - since momentarily bad
  // parameters could completely destroy the structure and make the successive minimization work on a crazy
  // structure.
  // boost::filesystem::copy_file(pwd.string()+"/CPU."+str_core+"/fort.90" ,
  //   pwd.string()+"/CPU."+str_core+"/geo." + str_cycle + "." + str_parID, 
  //     boost::filesystem::copy_option::overwrite_if_exists);

  // LOCAL MIN: return gradient
  //if (!grad_out.empty()) {
  //    grad_out = eval_numgrad(active_params, cycle, iter, p);
  //};
  evalfitdouble=evalfit;
  fitness=evalfit;
  return evalfitdouble;
#endif

#ifndef WITH_MPI

  int cycle;
  int iter;
  int parid;

  cycle = state.cycle;
  iter = state.iter;
  parid = state.parid;

  //boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);

  // count total # func evaluations
  funceval = funceval + 1;

  // prepare fort.3 files
  boost::filesystem::copy_file("geo." + str_cycle + "." + str_parID,
    "fort.3", boost::filesystem::copy_option::overwrite_if_exists);

  // check if ffield is LG or not. execute correct tapreaxff accordingly
  if (lg_yn == true) {
      write_ffield_lg(active_params, cycle, iter, parid);
  } else {
      write_ffield(active_params, cycle, iter, parid);
  };

  // check if tapreaxff is present
  std::ifstream fin5("tapreaxff");
  if (fin5.fail()) {
      cout << "tapreaxff executable not found. Aborting! \n";
      fin5.close();
      exit(EXIT_FAILURE);
  }
  fin5.close();

  // prepare mandatory files before executing tapreaxff
  boost::filesystem::copy_file("ffield",
      "fort.4", boost::filesystem::copy_option::overwrite_if_exists);

  if (fixcharges == true){
      boost::filesystem::copy_file("charges",
      "fort.26", boost::filesystem::copy_option::overwrite_if_exists);
  };

  /* execute tapreaxff */
  //arguments for tapreaxff, will run: tapreaxff                                                                                                                                  
  char *args[3] = { "./tapreaxff", "", NULL} ;
  pid_t c_pid, pid;
  int status;

  /* create a child process */
  c_pid = fork();

  if (c_pid == 0){
    /* CHILD */
    char * filename = "run.log";
    int outfile = open(filename, O_CREAT | O_WRONLY, S_IRWXU);
    if (outfile == -1){
          fprintf(stderr, "Error: failed to create file %s\n", filename);
    }else{
          /* redirect the standard output from this process to the file. */
          if(dup2(outfile, STDOUT_FILENO) != STDOUT_FILENO){
            fprintf(stderr, "Error: failed to redirect standard output\n");
          }
          /* redirect tapreaxff stdout to outfile*/
          dup2 (outfile, STDOUT_FILENO);
          /* redirect tapreaxff stderr to /dev/null */
          dup2(open("/dev/null", 0), 2);

        /* printf("Child: executing args\n"); */
        // execute args                                                                                                                                                               
        execvp( args[0], args);
        // only get here if exec failed                                                                                                                                             
        perror("execve failed");
        wait(&status);
    };
  }else if (c_pid > 0){
    /* PARENT */

    if( (pid = wait(&status)) < 0){
      perror("wait");
      _exit(1);
    };
    //printf("Parent: finished\n");

  }else{
    perror("fork failed");
    _exit(1);
  };

  if (WIFSIGNALED (status)) {
    cout << "Error: tapreaxff exited abnormaly on CPU:" << core << "\n";
    exit(EXIT_FAILURE);
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
      if (str.at(0) == '*' || !isdigit(str.at(0))) {
        evalfit = numeric_limits < double > ::infinity();
        if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
             pfitness = evalfit;
        };
      } else {
          if (regular == 1 || regular == 2) {
             // convert to double
             fitness = stod(str) + get_reg();
          } else {
             fitness = stod(str);
          };
      file13.close();
      boost::filesystem::remove("fort.13");
      };
  };
    // Note: do not update the geometry file during iterations. Each member should use one geo file throughout
    // the training. Assuming we start with DFT_optimized (or sensible structures), and that we use some small
    // number of structural minimizations (3-10), the sensible structures won't change much, so no need to
    // provide previous structure as initial structure for the next round of iteration - since momentarily bad
    // parameters could completely destroy the structure and make the successive minimization work on a crazy
    // structure.
    // boost::filesystem::copy_file(pwd.string() + "/fort.90" ,
    //   pwd.string() + "/geo."+str_cycle+"." + str_parID ,
    //      boost::filesystem::copy_option::overwrite_if_exists);

    // LOCAL MIN: return gradient
    //if (!grad_out.empty()) {
    //    grad_out = numgrad;
    //};


  return fitness;
#endif
};

vector <double> Par::eval_numgrad(const vector <double> &active_params, void * my_func_data) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered eval_numgrad()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered eval_numgrad()" << endl;
};
#endif

#ifdef WITH_MPI
  int cycle;
  int iter;
  int parid;

  cycle = state.cycle;
  iter = state.iter;
  parid = state.parid;

  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_core = std::to_string(core);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);

  // note: diff cannot be lower than 1e-4 because it determines the written parameters in the ffield file
  // which are then read with an assumed fixed format of 4 decimal digits in tapreaxff!
  double diff = 0.0001;
  double fitplus;
  double fitminus;
  vector <double> local_params;
  vector <double> grad(dim, 0.0);
  local_params = active_params;

  struct {
    int cyc;
    int it;
    int p;
  } myfdata;
  myfdata.cyc = cycle;
  myfdata.it = iter;
  myfdata.p = parid;

  // save fort.4 before writing new ffield
  if (boost::filesystem::exists("CPU." + str_core + "/fort.4")) {
      boost::filesystem::copy_file("CPU." + str_core + "/fort.4",
         "CPU." + str_core + "/fort.4.save", boost::filesystem::copy_option::overwrite_if_exists);
  };

  // perform central finite difference
  for (int i = 0; i < dim; i++) {
     local_params.at(i) = active_params.at(i) + diff;
  
     fitplus = eval_fitness(local_params, &myfdata);
     //if (verbose == true) {
     //   cout << "CPU: " << core << " fitplus[" << i << "] = " << boost::format("%18.6f") %fitplus << endl;
     //};

     local_params.at(i) = active_params.at(i) - diff;

     fitminus = eval_fitness(local_params, &myfdata);
     //if (verbose == true) {
     //   cout << "CPU: " << core << " fitminus[" << i << "] = " << boost::format("%18.6f") %fitminus << endl;
     //};

     grad.at(i) = (fitplus - fitminus)/(2.0*diff);
     // transform grad.at(i) to physical again
     //grad.at(i) = grad.at(i)*(maxdomain.at(i) - mindomain.at(i)) + mindomain.at(i);

     //if (verbose == true) {
     //   cout << "CPU: " << core << " grad[" << i << "] = " << boost::format("%18.6f") %grad.at(i) << endl;
     //};
  };
  
  // count total # func evaluations
  funceval = funceval + 2*dim;

  // retrive back saved ffield
  if (boost::filesystem::exists("CPU." + str_core + "/fort.4.save")) {
      boost::filesystem::copy_file("CPU." + str_core + "/fort.4.save",
        "CPU." + str_core + "/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  };

  return grad;
#endif

#ifndef WITH_MPI
  int cycle;
  int iter;
  int parid;

  cycle = state.cycle;
  iter = state.iter;
  parid = state.parid;

  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_parID = std::to_string(parid);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);

  // note: diff cannot be lower than 1e-4 because it determines the written parameters in the ffield file
  // which are then read with an assumed fixed format of 4 decimal digits in tapreaxff!
  double diff = 0.0001;
  double fitplus;
  double fitminus;
  vector <double> local_params;
  vector <double> grad(dim, 0.0);
  local_params = active_params;

  struct {
    int cyc;
    int it;
    int p;
  } myfdata;
  myfdata.cyc = cycle;
  myfdata.it = iter;
  myfdata.p = parid;

  // save fort.4 before writing new ffield
  boost::filesystem::copy_file("fort.4",
    "fort.4.save", boost::filesystem::copy_option::overwrite_if_exists);

  // perform central finite difference
  for (int i = 0; i < dim; i++) {
     local_params.at(i) = active_params.at(i) + diff;

     fitplus = eval_fitness(local_params, &myfdata);
     //if (verbose == true) {
     //   cout << "fitplus[" << i << "] = " << boost::format("%18.6f") %fitplus << endl;
     //};

     local_params.at(i) = active_params.at(i) - diff;

     fitminus = eval_fitness(local_params, &myfdata);
     //if (verbose == true) {
     //   cout << "fitminus[" << i << "] = " << boost::format("%18.6f") %fitminus << endl;
     //};

     grad.at(i) = (fitplus - fitminus)/(2.0*diff);
     //if (verbose == true) {
     //   cout << "grad[" << i << "] = " << boost::format("%18.6f") %grad.at(i) << endl;
     //};
  };

  // count total # func evaluations
  funceval = funceval + 2*dim;

  // retrive back saved ffield
  boost::filesystem::copy_file("fort.4.save",
    "fort.4", boost::filesystem::copy_option::overwrite_if_exists);

  return grad;
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

void Par::set_posdim(int i, double posx){
  pos.at(i) = posx;
};

void Par::set_vel(vector < double > vel_of_best_particle) {
  vel = vel_of_best_particle;
};

void Par::set_veldim(int i, double velx){
  vel.at(i) = velx;
};

void Par::dropout (double dropprobability) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered dropout()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered dropout()" << endl;
};
#endif


  dropped_dimns.clear();
  for (int i=0; i < dim; i++) {
     std::uniform_real_distribution <double> unidist(0.0, 1.0);
     if (unidist(generator) < dropprobability) {
        dropped_dimns.push_back(i); // stores dimensions to be dropped
        pos.at(i) = 0.0001;  // dropped dimension
        // standardize again to counter the transformation to physical params when writing ffield and printpos
        // since during localmin (iterate) once-standardized positions are necessary, a separate initialization (0.0001)
        // of drooped positions made in 'iterate'. 
        pos.at(i) = ( pos.at(i) - mindomain.at(i) ) / (maxdomain.at(i) - mindomain.at(i));
        //cout << "CPU " << core << " dropped standardized pos.at(" << i << ") = " << pos.at(i) << endl;
     }else{
        // multiply by 1/dropprobability all the other params so as to preserve total params number in the ffield.
        // now it should be possible to use the gbest ffield for the validation/test sets as-is.
        // step 1: multiply by 1/dropprobability the *physical* params
        pos.at(i) = (1.0/dropprobability)*( pos.at(i)*(maxdomain.at(i) - mindomain.at(i)) + mindomain.at(i) );
        // step 2: update domains
        mindomain.at(i) = (1.0/dropprobability)*mindomain.at(i);
        maxdomain.at(i) = (1.0/dropprobability)*maxdomain.at(i);
        // step 3: convert back to standardized form based on new domains
        pos.at(i) = ( pos.at(i) - mindomain.at(i) ) / (maxdomain.at(i) - mindomain.at(i));
        //cout << "CPU " << core << " non-dropped standardized pos.at(" << i << ") = " << pos.at(i) << endl;
        // step 4: convert new domains back to original before next dropout
        mindomain.at(i) = dropprobability*mindomain.at(i);
        maxdomain.at(i) = dropprobability*maxdomain.at(i);
        // NOTE: should multiplication by 1/dropprobability moved to printpos and writeffield sections??
     };
  };

};

double Par::get_reg() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_reg()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_reg()" << endl;
};
#endif
  double reg = 0.0;

  // vector to store the positions that will be converted below to their inverses to use in regularization
  vector <double> pos_for_reg;
  pos_for_reg.clear();

  // if position belongs to the subset of positions that are used as inverses by ReaxFF (i.e. 1/parameter),
  // use the inverse of that parameter into the sum of regularization
  for (vector <int> inversevalue : inversep) {
    for (int i = 0; i < ffline.size(); i++) {
      if ( (ffline.at(i)+1 == inversevalue.at(0)) && (ffcol.at(i)+1 == inversevalue.at(1)) ) {
         pos_for_reg.push_back( 1.0/(pos.at(i)*(maxdomain.at(i) - mindomain.at(i)) + mindomain.at(i)) );
      }else {
         pos_for_reg.push_back( pos.at(i)*(maxdomain.at(i) - mindomain.at(i)) + mindomain.at(i) );
      };
    };
  };

  // calculate L1 penalty
  if (regular == 1) {
    for (int i = 0; i < dim; i++) {
        reg = reg + fabs( pos_for_reg.at(i) );
    };
    reg = hlambda*reg;
  };

  // calculate L2 penalty
  if (regular == 2) {
    for (int i = 0; i < dim; i++) {
      reg = reg + pow(pos_for_reg.at(i),2);
    };
    reg = hlambda*reg;
  };
  return reg;

};




// ---------- Definitions of Swarm class member functions ---------- //
//
//------------------------------------------------------------------ //

Swarm::Swarm() {

};

void Swarm::read_icharg_control() {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered read_icharg_control()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered read_icharg_control()" << endl;
};
#endif
  // Read icharg value from reax control file
  string str_core = std::to_string(core);
  std::ifstream control_file("control");
  if (control_file.fail()) {
#ifdef WITH_MPI
    cout << "Unable to open control file on CPU " << core << ". \n";
    control_file.close();
    MPI_Abort(MPI_COMM_WORLD,5);
#endif
#ifndef WITH_MPI
    cout << "Unable to open control file. \n";
    control_file.close();
    exit(EXIT_FAILURE);
#endif
  } else{
    std::string line;
    int numlines = 0;
    // for each line
    while (std::getline(control_file, line)) {
      boost::trim(line);
      string val;
      std::istringstream lineStream(line);
      // for each value in line
      while (lineStream >> val) {
        if (val == "icharg") {
           continue;
        };
        if (val == "5") {
           fixcharges = true;
           break;
        };
      };
     };
  };
  control_file.close();
};

void Swarm::get_userinp(){
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_userinp()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_userinp()" << endl;
};
#endif

#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(core);

  if (core == 0){
    vector <string> tempinput;
    std::ifstream fin("inp_flocky.in");
    if (fin.fail()) {
      cout << "Error: Unable to open 'inp_flocky.in' file \n";
      fin.close();
      MPI_Abort(MPI_COMM_WORLD,6);
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
          using boost::is_any_of;
          boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
          tempinput.push_back(values.at(0));
        };
      };
      fin.close();
    };
    istringstream(tempinput.at(0)) >> verbose;
    istringstream(tempinput.at(1)) >> lg_yn;
    istringstream(tempinput.at(2)) >> ptrainset;
    if (ptrainset == 0) { 
    //note: 18.02.20 letting ptrainset be equal to 0 causes hangs in the generation of ffield.gbest
    //      it appears that setting it to 1 instead, solves the issue (and has the same effect
    //      of not parallelising the training set. 
        ptrainset = 1;
    };
    if (ptrainset > numcores && ensembleave == false && preppath == 0) {
      cout << "Error: Number of allocated processors < number of training set sub-units" << endl;
      MPI_Abort(MPI_COMM_WORLD,7);
    }else if (ptrainset > 1) {
      swarmcores.clear();
      reaxffcores.clear();
      for (int i = 0; i < numcores; i=i+ptrainset) {
          swarmcores.push_back(i);
          for (int j = 1; j < ptrainset; j++) {
            reaxffcores.push_back(i+j);
          };
      };
    };
    istringstream(tempinput.at(3)) >> contff;
    istringstream(tempinput.at(4)) >> perc_yn;
    istringstream(tempinput.at(5)) >> perc;
    istringstream(tempinput.at(6)) >> localmin;
    istringstream(tempinput.at(7)) >> lm_iter_max;
    istringstream(tempinput.at(8)) >> lm_err_tol;
    istringstream(tempinput.at(9)) >> regular;
    istringstream(tempinput.at(10)) >> dropvalue;
    istringstream(tempinput.at(11)) >> hlambda;
    istringstream(tempinput.at(12)) >> ofit;
    istringstream(tempinput.at(13)) >> uq;
    istringstream(tempinput.at(14)) >> NumP;
    if (ensembleave == false && preppath == 0) {
       if (NumP < numcores) {
         cout << "Error: Number of swarm members < number of allocated processors." << endl;
         MPI_Abort(MPI_COMM_WORLD,7);
       } else {
         totmembers = NumP;
         NumP = int(floor(NumP / numcores));
       };

       if (ptrainset < 0 || (ptrainset > 1 && mod(totmembers,ptrainset) != 0)) {
          cout << "Error: Number of training set sub-units < 0 or not a divisor of number of allocated processors." << endl;
          MPI_Abort(MPI_COMM_WORLD,7);
       };
    };
    istringstream(tempinput.at(15)) >> c1;
    istringstream(tempinput.at(16)) >> c2;
    istringstream(tempinput.at(17)) >> inertiamax;
    istringstream(tempinput.at(18)) >> inertiamin; 
    istringstream(tempinput.at(19)) >> faili;
    istringstream(tempinput.at(20)) >> levyscale;
    istringstream(tempinput.at(21)) >> freq;
    istringstream(tempinput.at(22)) >> maxiters;
    istringstream(tempinput.at(23)) >> maxcycles;
    istringstream(tempinput.at(24)) >> ensembleave;
    istringstream(tempinput.at(25)) >> preppath;
    istringstream(tempinput.at(26)) >> econv;
    istringstream(tempinput.at(27)) >> frqconv;
  };  // close if core==0
 

  // broadcast between all processes from process 0
  // MPI_Bcast must be visible to all processes!
  MPI_Bcast( & lg_yn, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & ptrainset, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & verbose, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
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
  MPI_Bcast( & localmin, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & lm_iter_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & lm_err_tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & regular, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & dropvalue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & hlambda, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & ofit, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & uq, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & ensembleave, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast( & preppath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( & econv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( & frqconv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // prepare dirs for each CPU process
  if (ensembleave == false && preppath == 0) {
     boost::filesystem::create_directory("CPU." + str_core);
     boost::filesystem::copy_file("tapreaxff", "CPU." + str_core + "/tapreaxff",
         boost::filesystem::copy_option::overwrite_if_exists);
     boost::filesystem::copy_file("ffield", "CPU." + str_core + "/ffield",
         boost::filesystem::copy_option::overwrite_if_exists);
     boost::filesystem::copy_file("control", "CPU." + str_core + "/control",
         boost::filesystem::copy_option::overwrite_if_exists);
     boost::filesystem::copy_file("geo", "CPU." + str_core + "/geo",
         boost::filesystem::copy_option::overwrite_if_exists);
  };

  // after broadcasting ensembleave all CPUs have the value of ensembleave
  // now it is possible to perform the following on every CPU
  if (ptrainset == 1 && (ensembleave == false) && preppath == 0) {
    boost::filesystem::copy_file("trainset.in", "CPU." + str_core + "/trainset.in",
        boost::filesystem::copy_option::overwrite_if_exists);
  };
  if (ptrainset > 1 && (ensembleave == false) && preppath == 0) {
    // define size of vectors for all the rest of the cores
    swarmcores.resize(numcores/ptrainset);
    reaxffcores.resize(numcores - (numcores/ptrainset));
    // now is possible to bcast swarmcores.data() and reaxffcores.data()
    MPI_Bcast( swarmcores.data(), swarmcores.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( reaxffcores.data(), reaxffcores.size(), MPI_INT, 0, MPI_COMM_WORLD);
  };

  if (ensembleave == false && preppath == 0) {
     // check if reaxff was set to run with fixed charges and require charges file
     read_icharg_control();

     boost::filesystem::ifstream charge_file("charges");
     //fixcharges = true;
     if (fixcharges == true && charge_file.fail()) {
       cout << "Error: 'control' uses icharg=5 (fixed charges) but no 'charges' file was found!" << endl;
       charge_file.close();
       MPI_Abort(MPI_COMM_WORLD,8);
     };
     charge_file.close();
  };

  if (fixcharges == true && (ensembleave == false) && preppath == 0) {
    boost::filesystem::copy_file("charges", "CPU." + str_core + "/charges",
      boost::filesystem::copy_option::overwrite_if_exists);
  };

  // create fort.20 file needed by tapreaxff
  if (ensembleave == false && preppath == 0) {
     boost::filesystem::ofstream iopt_file("CPU." + str_core + "/fort.20");
     iopt_file << "0";
     iopt_file.close();
  };

  // create fort.35 file needed by tapreaxff
  if (ensembleave == false && preppath == 0) {
     ofstream fort35_file("CPU." + str_core + "/fort.35");
     fort35_file << "23434.1" << endl;
     fort35_file.close();
  };
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

    istringstream(tempinput.at(0)) >> verbose;
    istringstream(tempinput.at(1)) >> lg_yn;
    istringstream(tempinput.at(2)) >> ptrainset;
    istringstream(tempinput.at(3)) >> contff;
    istringstream(tempinput.at(4)) >> perc_yn;
    istringstream(tempinput.at(5)) >> perc;
    istringstream(tempinput.at(6)) >> localmin;
    istringstream(tempinput.at(7)) >> lm_iter_max;
    istringstream(tempinput.at(8)) >> lm_err_tol;
    istringstream(tempinput.at(9)) >> regular;
    istringstream(tempinput.at(10)) >> dropvalue;
    istringstream(tempinput.at(11)) >> hlambda;
    istringstream(tempinput.at(12)) >> ofit;
    istringstream(tempinput.at(13)) >> uq;
    istringstream(tempinput.at(14)) >> NumP;
    istringstream(tempinput.at(15)) >> c1;
    istringstream(tempinput.at(16)) >> c2;
    istringstream(tempinput.at(17)) >> inertiamax;
    istringstream(tempinput.at(18)) >> inertiamin;
    istringstream(tempinput.at(19)) >> faili;
    istringstream(tempinput.at(20)) >> levyscale;
    istringstream(tempinput.at(21)) >> freq;
    istringstream(tempinput.at(22)) >> maxiters;
    istringstream(tempinput.at(23)) >> maxcycles;
    istringstream(tempinput.at(24)) >> ensembleave;
    istringstream(tempinput.at(25)) >> preppath;
    istringstream(tempinput.at(26)) >> econv;
    istringstream(tempinput.at(27)) >> frqconv;

  // check if reaxff was set to run with fixed charges and require charges file
  read_icharg_control();
  boost::filesystem::ifstream charge_file("charges");
  if (charge_file.fail()) {
    cout << "Error: 'control' uses icharg=5 (fixed charges) but no 'charges' file was found!" << endl;
    charge_file.close();
    exit(EXIT_FAILURE);
  };
  charge_file.close();

  // create fort.20 file needed by tapreaxff
  boost::filesystem::ofstream iopt_file("fort.20");
  iopt_file << "0";
  iopt_file.close();

  // create fort.35 file needed by tapreaxff
  ofstream fort35_file("fort.35");
  fort35_file << "23434.1" << endl;
  fort35_file.close();

  // prepare fort.111 (forces) file if it exists (for training FORCES in trainset)
  //if (boost::filesystem::exists("forgeo")) {
  //   boost::filesystem::copy_file("forgeo", "fort.111",
  //       boost::filesystem::copy_option::overwrite_if_exists);
  };


#endif
};


vector <double> Swarm::get_com(Swarm newSwarm) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_com()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_com()" << endl;
};
#endif

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
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_disp()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_disp()" << endl;
};
#endif
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
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_worse()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_worse()" << endl;
};
#endif
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
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered get_best()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered get_best()" << endl;
};
#endif

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

void Swarm::MakeFilesfromPath(Swarm & newSwarm) {
#ifdef WITH_MPI
if (core == 0) {
  if (verbose == true) {
     cout << "CPU: " << core << " entered MakeFilesfromPath()" << endl;
  };
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered MakeFilesfromPath()" << endl;
};
#endif
if (core == 0) {
   if (preppath > 0 && ensembleave == false) {
      for (int p = 0; p < 1; p++) {
         Par NewPar;
         newSwarm.AddPar(NewPar);
      };

      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      boost::filesystem::ifstream pathinfo_file("path.info");
      boost::filesystem::ifstream vectordump_file("vector.dump.path");

      if (pathinfo_file.fail() || vectordump_file.fail()) {
        cout << "Error: 'path.info or vector.dump.path' files not found!" << endl;
        MPI_Abort(MPI_COMM_WORLD,8);
        pathinfo_file.close();
      }else{
      log << "\n";
      log << "Preparing training files from path.info and vector.dump.path files..." << endl;
      log << "Note: 6 smallest |frequencies| will be ignored. Assuming non-linear molecule." << endl;
      };

      int natoms=preppath;
      std::ifstream fin("path.info");
      std::ifstream fin2("vector.dump.path");



      // struct to store energies, symmetries, mass-weighted hessian eigenvalues
      // and coordinates as given in path.info file
      struct STpoint{
        double energy;
        string symmetry;
        vector < vector <double> > Mwhessevals;
        vector < vector <double> > coordinates;
      };

      vector <STpoint> pathinfo;
      pathinfo.clear();
      STpoint mystpoint;


      std::string line;
      int numlines = 0;
      // for each line of path.info file
      while (std::getline(fin, line)) {
        numlines++;
        // create a new row
        std::vector < double > lineData;
        string val;
        std::istringstream lineStream(line);
        // add 1st row in block to energy
        if (numlines == 1) {
           mystpoint.energy = stod(line);
        };
        // add 2nd row in block to symmetry
        if (numlines == 2) {
           mystpoint.symmetry = line;
        };
        // for each value in line
        while (lineStream >> val && numlines > 2) {
             lineData.push_back(stod(val));
        };
        // add next rows to Mwhessevals
        if (numlines > 2 && numlines < 3+natoms) {
           mystpoint.Mwhessevals.push_back(lineData);
        };
        // add next rows to coordinates
        if (numlines > 2+natoms && numlines < 3+2*natoms) {
           mystpoint.coordinates.push_back(lineData);
        };
        // finished with 1st block (stationary point data)
        // reset counter and start again in next block
        // after adding the struct for current block to array pathinfo
        if (numlines == 2*natoms+2) {
           numlines = 0;
           pathinfo.push_back(mystpoint);
           mystpoint = {};
        };
      };


      // struct to store freqs and normal mode vectors as given in vector.dump.path file
      struct frqandmodes{
        double freqfromdump;
        vector <double> normalmodes; // size [3*natoms]
      };

      vector <frqandmodes> vdumpinfo; // vector to contain all frequencies and associated modes of a given stationary point
      vector < vector <frqandmodes> > totvdumpinfo; // vector to contain all previous vectors (i.e. for all stationary points)
      vdumpinfo.clear();
      totvdumpinfo.clear();

      frqandmodes myvdump;

      numlines=0;
      int numstructs=0;
      // for each line of vector.dump.path file
      while (std::getline(fin2, line)) {
        // create a new row
        std::vector < double > lineData;
        string val;
        std::istringstream lineStream(line);

        // add frequency element into struct myvdump
        if (mod(numlines,natoms+1)==0) {
           myvdump.freqfromdump = stod(line);
           //cout << "freq: " << stod(line) << endl;
        };

        // fill up normal modes of one frequency
        //myvdump.normalmodes.clear();
        while ( lineStream >> val && (mod(numlines,natoms+1) != 0) ) {
             myvdump.normalmodes.push_back(stod(val));
             //cout << "nmodes pushed back: " << stod(val) << endl;
        };
        //cout << "numlines = " << numlines << endl;
        //cout << "number of modes in myvdump.normalmodes = " << myvdump.normalmodes.size() << endl; 
        numlines++;

        // fill up the nmodes of current stationary point
        if (myvdump.normalmodes.size() == 3*natoms) {
           //cout << "myvdump.normalmodes reached 3*natoms modes. filling vdumpinfo" << endl;
           vdumpinfo.push_back(myvdump);
           myvdump = {};
        };

        // fill up the nmodes of all stationary points
        if (vdumpinfo.size() == 3*natoms) {
           //vdumpinfo.push_back(myvdump);
           //numstructs++;
           //cout << "finished inserting freqs and nmodes of stationary point #" << numstructs << endl;
           totvdumpinfo.push_back(vdumpinfo);
           myvdump = {};
           vdumpinfo.clear();
        };

      };

      //cout << "============ debug nmodes ===============" << endl;
      //cout << "totvdumpinfo.size = " << totvdumpinfo.size() << endl;
      //for (int k=0; k < totvdumpinfo.size(); k++) {
      //    cout << "struct #" << k << endl;
      //    for (int v=0; v < totvdumpinfo.at(k).size(); v++) {
      //        cout << "freq #" << v << endl;
      //        for (int r=0; r < totvdumpinfo[k].at(v).normalmodes.size(); r++) {
      //            cout << totvdumpinfo[k].at(v).normalmodes.at(r) << endl;
      //        };
      //    };
      //};


      // Prepare chopped frequencies (sqrt of Mwhessevals) file
      // remove 6 smallest (in absolute value) Mwhessevals and store
      // in tempfreq. 
      // note: needs to be modified for linear molecules.
      boost::filesystem::ofstream freqfile;
      vector <double> ::reverse_iterator col;
      vector <vector <double>> ::reverse_iterator row;
      boost::filesystem::create_directory("frequencies");

      // define a struct to hold a frequency and it's ID in Mwhessevals
      struct freqandID {
         int freqID;
         double freq;
      };

      // define a "greater than" comparator for freqandID type elements to use in max heap
      struct freqandID_greater_than {
          bool operator()(freqandID const& a, freqandID const& b) const {
              return a.freq > b.freq;
          };
      };

      // define a "less than" comparator for freqandID type elements to use in sort
      struct freqandID_less_than {
          bool operator()(freqandID const& a, freqandID const& b) const {
              return a.freq < b.freq;
          };
      };

      freqandID mytempfreq; // a single pair of ID and frequency
      vector <freqandID> tempfreq; // all pairs of ID-frequency

      tempfreq.clear();
      int chunk;
      int endcol;
      for (int k=0; k < pathinfo.size(); k++) {
          chunk=0;
          endcol=6;
         for (int i=0; i < pathinfo[k].Mwhessevals.size(); i++) {
            for (int j=0; j < pathinfo[k].Mwhessevals.at(i).size(); j++) {
              //add the first 3*natoms-6 Mwhessevals to the vector
              if (tempfreq.size() < 3*natoms-6)
              {
                  // convert mass-weighted hessian eigenvalues to frequencies
                  // convention: negative eigenvalue --> negative frequency
                  if (pathinfo[k].Mwhessevals.at(i).at(j) < 0.0) {
                       //tempfreq.push_back(-1.0*pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5));
                       mytempfreq.freqID=i+j;
                       mytempfreq.freq=-1.0*pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5);
                  }else{
                       //tempfreq.push_back(pow(pathinfo[k].Mwhessevals.at(i).at(j),0.5));
                       mytempfreq.freqID=i+j;
                       mytempfreq.freq=pow(pathinfo[k].Mwhessevals.at(i).at(j),0.5);
                  };
                  // fill in the tempfreq with a new mytempfreq struct
                  tempfreq.push_back(mytempfreq);

                  if ( tempfreq.size() == 3*natoms-6 )
                      // make the max_heap
                      std::make_heap(tempfreq.begin(), tempfreq.end(),freqandID_greater_than());
                  continue;
              }
              //cout << "tempfreq.size = " << tempfreq.size() << endl;

              // now check if the next element is larger than the top of the heap
              //cout << "structure #" << k << endl;
              //cout << "new candidate: " << pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5) << endl;
              //cout << "front of heap: " << fabs(tempfreq.front().freq) << endl;
              if (fabs(tempfreq.front().freq) < pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5))
              {
                  //cout << "exchange!" << endl;

                  // remove the front of the heap by placing it at the end of the vector
                  std::pop_heap(tempfreq.begin(), tempfreq.end(),freqandID_greater_than());

                  // get rid of that item now 
                  tempfreq.pop_back();

                  // add the new item after converting mass-weighted hessian eigenvalues to frequencies
                  // convention: negative eigenvalue --> negative frequency
                  if (pathinfo[k].Mwhessevals.at(i).at(j) < 0.0) {
                       //tempfreq.push_back(-1.0*pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5));
                       mytempfreq.freqID=i+j;
                       mytempfreq.freq=-1.0*pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5);
                       //cout << "pushed back into tempfreq: " << -1.0*pow(fabs(pathinfo[k].Mwhessevals.at(i).at(j)),0.5) << endl;
                  }else{
                       //tempfreq.push_back(pow(pathinfo[k].Mwhessevals.at(i).at(j),0.5));
                       mytempfreq.freqID=i+j;
                       mytempfreq.freq=pow(pathinfo[k].Mwhessevals.at(i).at(j),0.5);
                       //cout << "pushed back into tempfreq: " << pow(pathinfo[k].Mwhessevals.at(i).at(j),0.5) << endl;
                  };
                  tempfreq.push_back(mytempfreq);
                  //cout << "front of heap after push_back but before push_heap: " << tempfreq.front() << endl;

                  // heapify
                  std::push_heap(tempfreq.begin(), tempfreq.end(),freqandID_greater_than());
                  //cout << "front of heap after push_heap with greater: " << tempfreq.front() << endl;
              };
            };
         };

         // sort the chopped (3*natoms - 6) frequencies in ascending order
         std::sort(tempfreq.begin(), tempfreq.end(),freqandID_less_than());

         // ========== PRINTING frequencies and normal modes in chunks of 6 columns =========== //
         // assemble all indices of chopped frequencies into vector `matches`
         vector <int> matches;
         matches.clear();
         for (int i=0; i < tempfreq.size(); i++){
            matches.push_back(tempfreq.at(i).freqID);
         };
         //cout << "total number of matches: " << matches.size() << endl;

         // print the chopped (3*natoms-6) sorted frequencies (tempfreq) to files
         freqfile.open( "frequencies/SP_" + std::to_string(k+1) + "_vib" , ofstream::app ) ;
         while (chunk < tempfreq.size()) {
              freqfile << boost::format("%11a%6x") %"frequencies"%"";
              for (int i=chunk; i<chunk+endcol; i++){
                  tempfreq.at(i).freq = tempfreq.at(i).freq*frqconv;
                  freqfile << boost::format("%8.2f%1x") %tempfreq.at(i).freq %"";
              };
              freqfile << endl;

              // print normal modes that belong to the frequency ID's
              for (int r=0; r<3*natoms; r++) {
                  //for (auto n : matches) {
                  freqfile << boost::format("%17x") %"";
                  for (int i=chunk; i<chunk+endcol; i++){
                      int n;
                      n=matches.at(i);
                      freqfile << boost::format("%8.5f%1x") %totvdumpinfo[k].at(n).normalmodes.at(r)%"";
                  };
                  freqfile << endl;
              };

              chunk=chunk+6;
              if (chunk+6 > tempfreq.size()) {
                 endcol = tempfreq.size() - chunk;
              };
         }; // done printing for stationary point
         tempfreq.clear();
         freqfile.close();
      };


      // Prepare geo file
      vector <double> ::iterator acol;
      vector <vector <double>> ::iterator arow;
      boost::filesystem::ofstream geofile;
      int l=0;
      for (int k=0; k < pathinfo.size(); k++) {
         geofile.open( "geo.path" , ofstream::app ) ;
         geofile << "BIOGRF 200\n";
         geofile << boost::format("%6s%8s") %"DESCRP" %("SP_"+std::to_string(k+1));
         geofile << "\n";
         for (arow = pathinfo.at(k).coordinates.begin(); arow != pathinfo.at(k).coordinates.end(); ++arow) {
             l++;
             geofile << boost::format("%6s %5d %-5s %3s %1s %5s") %"HETATM" %l %"C" %"   " %" " %"   ";
             for (acol = arow->begin(); acol != arow->end(); ++acol) {
                 geofile << boost::format("%15.10f") %*acol;
             }
             geofile << boost::format(" %5s%3d%2d %8.5f\n") %"C" %"0" %"0" %"0.00000";
         }
         l=0;
         geofile << "END\n" << endl;
         geofile.close();
      };


      // Prepare trainset.in file
      boost::filesystem::ofstream trainfile;
      trainfile.open( "trainset.in.path", ofstream::app);
      // energy section
      trainfile << "ENERGY" << endl;
      for (int k=0; k < pathinfo.size(); k++) {
          // make a baseline - subtract energy of first stationary point from all energies
          double ediff = pathinfo.at(k).energy - pathinfo.at(0).energy;
          ediff=ediff*econv;
          trainfile << " 1.0 + SP_" + std::to_string(k+1) + "/1 - SP_1/1" + "        " << boost::format("%5.4f") %ediff << endl;
      };
      trainfile << "ENDENERGY" << endl;
      // frequencies section
      trainfile << "FREQUENCIES" << endl;
      for (int k=0; k < pathinfo.size(); k++) {
          trainfile << " SP_" + std::to_string(k+1) + "  10.0  10.0 " + "SP_" + std::to_string(k+1) + "_vib" << endl;
      };
      trainfile << "ENDFREQUENCIES" << endl;
      trainfile.close();

      log << "Done!" << endl;
      log.close();
      pathinfo_file.close();

   };
};
};

void Swarm::MakeEnsembleFF(Swarm & newSwarm) {
#ifdef WITH_MPI
if (core == 0) {
  if (verbose == true) {
     cout << "CPU: " << core << " entered MakeEnsembleFF()" << endl;
  };
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered MakeEnsembleFF()" << endl;
};
#endif
#ifdef WITH_MPI
if (core == 0) {
   if (ensembleave == true && preppath == 0) {
      for (int p = 0; p < 1; p++) {
         Par NewPar;
         newSwarm.AddPar(NewPar);
      };

      boost::filesystem::ofstream log("log.flocky", ofstream::app);
      log << "\n";
      log << "Preparing an ensemble-averaged force field." << endl;

      const std::string target_path( "ffbags" );
      const std::regex my_filter ( "ffield.gbest.*" );

      vector <vector <string>> myffield;
      
      numffbags = 0;
      vector <string> all_matching_files;
      
      boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
      for( boost::filesystem::directory_iterator i( target_path ); i != end_itr; ++i )
      {
          // Skip if not a file
          if( !boost::filesystem::is_regular_file( i->status() ) ) continue;
          
          std::smatch what;
          // skip if no match
          if (!regex_search (i->path().filename().string(),what,my_filter)) continue;

          // File matches, store it
          all_matching_files.push_back( i->path().filename().string() );
          numffbags++;
      }
      log << "\nflocky found an ensemble of " << numffbags << " gbest force fields to average." << endl;

      // loop over files of gbest ffields and store them
      newSwarm.GetPar(0).ffensemble.clear();
      for (auto myfile : all_matching_files) {
         myffield.clear();
         std::ifstream fin("ffbags/"+myfile);
         
         std::string line;
         int numlines = 0;
    
         // for each line
         while (std::getline(fin, line)) {
           //cout << line << endl;
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
           myffield.push_back(lineData);
         };
         fin.close();
         // store gbest ffields in ffensemble
         newSwarm.GetPar(0).ffensemble.push_back(myffield);
      };
      
      if (lg_yn == true) {
         newSwarm.GetPar(0).write_aveffield_lg();
      }else{
         newSwarm.GetPar(0).write_aveffield();
      };

      log << "\n";
      log << "An ensemble-averaged force field is ready!\n" << endl;
      log.close();
   };
};
#endif

#ifndef WITH_MPI
#endif
};

void Swarm::Populate(Swarm & newSwarm, int cycle) {
#ifdef WITH_MPI
if (verbose == true) {
  cout << "CPU: " << core << " entered Populate()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
  cout << "entered Populate()" << endl;
};
#endif

#ifdef WITH_MPI
if (ensembleave == false && preppath == 0) {
    if (core == 0 && verbose == true) {
       cout << "CPU: " << core << " entered Populate!" << endl; 
    };
    
     // define a new sub-communicator for reaxffcores and swarmcores //
        MPI_Comm ACTIVESWARM;  // swarmcores
        MPI_Comm PASSIVESWARM; // reaxffcores
        MPI_Comm *newcomm;
        int color;
        if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
            color = 444;
            newcomm = &ACTIVESWARM;
        };
        if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
            color = 333;
            newcomm = &PASSIVESWARM;
        };
        MPI_Comm_split( MPI_COMM_WORLD, color, core, newcomm );
    
        if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
             MPI_Comm_rank( ACTIVESWARM, &mycore_swarmcore);
             MPI_Comm_size( ACTIVESWARM, &size_swarmcores);
        };
    
        if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
             MPI_Comm_rank( PASSIVESWARM, &mycore_reaxffcore);
             MPI_Comm_size( PASSIVESWARM, &size_reaxffcores);
        };
    
    
      boost::filesystem::path pwd(boost::filesystem::current_path());
      string str_core = std::to_string(core);
      if (core == 0) {
        boost::filesystem::ofstream log("log.flocky", ofstream::app);
        log << "\n";
        log << "Swarm generation started. Please wait." << endl;
        log.close();
      };
      funceval = 0; // clear counter for cycles
    
      // find inverse parameters and store them in inversep
      if (regular == 1 || regular == 2) {
         get_inversep();
      };
    
      // ---------------------------------------------- //
      //     POPULATE: MAIN LOOP OVER SWARM MEMBERS
      // ---------------------------------------------- //
      for (int p = 0; p < NumP; p++) {
        string parID = std::to_string(p);
        string str_cycle = std::to_string(cycle);
    
        Par NewPar;
        newSwarm.AddPar(NewPar);
        // initial gbfit is INF for all processes and all swarm members
        gbfit = numeric_limits < double > ::infinity();
        // initial gbpos is 0.0 for all processes and all swarm members
        gbpos.clear();
        for (int m=0; m < dim; m++){
          gbpos.push_back(0.0);
    
        };
    
        if (ptrainset > 1) {
         // generate trainset and geo subsets
         newSwarm.GetPar(p).write_trainset();
         newSwarm.GetPar(p).write_geo();
        };
    
        // cp geo to geo.parID so each particle works with its own geo file
        boost::filesystem::copy_file("CPU." + str_core + "/geo",
          "CPU." + str_core + "/geo." + str_cycle + "." + parID,
            boost::filesystem::copy_option::overwrite_if_exists);
        
    
        if (core == 0) {
          // If contff == y, then take force field's current values for the position of particle 0 (others are random)
          if (contff == true) {
            vector < double > ffpos;
            double tempstand;
            ffpos.clear();
            for (int m = 0; m < dim; m++) {
                ffpos.push_back(stod(newSwarm.GetPar(0).ffieldmat.at(newSwarm.GetPar(0).ffline.at(m)).at(newSwarm.GetPar(0).ffcol.at(m))));
            };
            newSwarm.GetPar(0).set_pos(ffpos);
            // reset parameters if outside domain before standardization possible
            newSwarm.GetPar(0).check_bounds_contff();
            // standardize positions
            for (int m = 0; m < dim; m++) {
                newSwarm.GetPar(0).pos.at(m) = ( newSwarm.GetPar(0).pos.at(m) - newSwarm.GetPar(0).mindomain.at(m) ) / ( newSwarm.GetPar(0).maxdomain.at(m) - newSwarm.GetPar(0).mindomain.at(m)  );
            };
          };
          contff = false;
        };
    
        // evaluate fitness and set bfit = curfit
        newSwarm.GetPar(p).state.cycle = cycle;
        iter = 0;
        newSwarm.GetPar(p).state.iter = iter;
        newSwarm.GetPar(p).state.parid = p;
    
        // dropout
        if (regular == 3) {
           newSwarm.GetPar(p).dropout(dropvalue);
        };
    
        // local minimization
        if (localmin == 1 || localmin == 2) {
           newSwarm.GetPar(p).iterate();
        }else{
           newSwarm.GetPar(p).set_fitness(newSwarm.GetPar(p).eval_fitness(newSwarm.GetPar(p).get_pos_vec(), this));
        };
    
        // if we parallelize the training set, update gbfit and gbpos only among swarmcores
        if (ptrainset > 1) {
           if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());
              if (newSwarm.GetPar(p).get_fitness() < gbfit) {
                gbfitfound = true;
                parid_gbfit = p;
                gbfit = newSwarm.GetPar(p).get_fitness();
                gbpos.clear();
                gbpos = newSwarm.GetPar(p).get_pos_vec();
                //write_ffield_gbest(core, cycle, iter, p);
              }else{
                 gbfitfound = false;
              };
              // cleaning ffield.tmp.* files after write_ffield_gbest already
              // copied the correct ffield.tmp.* file as the ffield.gbest.*.0.*
              //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp." + str_cycle+".0." + "." + parID);
           };
        }else {
          newSwarm.GetPar(p).set_fitness(newSwarm.GetPar(p).get_fitness());
          newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());
          if (newSwarm.GetPar(p).get_fitness() < gbfit) {
            gbfitfound = true;
            parid_gbfit = p;
            gbfit = newSwarm.GetPar(p).get_fitness();
            gbpos.clear();
            gbpos = newSwarm.GetPar(p).get_pos_vec();
            //write_ffield_gbest(core, cycle, iter, p);
          }else{
             gbfitfound = false;
          };
          // cleaning ffield.tmp.* files after write_ffield_gbest already
          // copied the correct ffield.tmp.* file as the ffield.gbest.*.0.*
          //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp." + str_cycle+".0." + "." + parID);
        };
    
      }; // done loop on members
    
      // if we parallelize the training set, do the following only among swarmcores
      if (ptrainset > 1) {
            // pair struct to hold the global best fitness across processes and its core rank
            struct {
              double tmp_fit;
              int tmp_cpu;
            }  min_vals_in[1], min_vals_out[1];
            // store current fit on each process
            min_vals_in[0].tmp_fit = gbfit;
            // store core id of that current process
            min_vals_in[0].tmp_cpu = mycore_swarmcore; 
            // get global best fitness *across processes* and corresponding core rank and store them in min_vals_out
            if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
               MPI_Allreduce( & min_vals_in, & min_vals_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, ACTIVESWARM);
            };
            // global best fitness across all processes
            gbfit = min_vals_out[0].tmp_fit;
            // core rank the above fitness came from
            cpuid_gbfit = min_vals_out[0].tmp_cpu;
            // broadcast contents of gbpos vector from rank cpuid_gbfit
            if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
               MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, ACTIVESWARM);
            };
    
            // pair struct to hold the global best fitness across processes and its parID
             // The parID is required in detection of overfitting to cp the correct ffield.gbest file
             struct {
               double tmp_fit2;
               int tmp_parid;
             } /*min_vals_in2, min_vals_out2;*/ min_vals_in2[1], min_vals_out2[1];
    
             // store current fit on each process
              min_vals_in2[0].tmp_fit2 = gbfit;
             // store par id of that current process
              min_vals_in2[0].tmp_parid = parid_gbfit;
             // get global best fitness *across processes* and corresponding parID and store them in min_vals_out
             if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
                MPI_Allreduce( & min_vals_in2, & min_vals_out2, 1, MPI_DOUBLE_INT, MPI_MINLOC, ACTIVESWARM);
             };
             // parid of the gbfit
             parid_gbfit = min_vals_out2[0].tmp_parid;
    
             // now convert cpuid_gbfit value from ACTIVESWARM to MPI_COMM_WORLD
             if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end() && gbfitfound == true) {
                 cpuid_gbfit = swarmcores.at(cpuid_gbfit);
             };
    
             // do write_ffield_gbest
             if (core == cpuid_gbfit && find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
                //cout << "CPU: " << core << " (cpuid_gbfit) writing ffield_gbest: ffield.tmp."+to_string(cycle)+"."+"0."+to_string(parid_gbfit) << endl;
                write_ffield_gbest(cpuid_gbfit, cycle, 0, parid_gbfit);
                boost::filesystem::remove("CPU." + std::to_string(core) + "/ffield.tmp."+std::to_string(cycle)+"."+"0."+std::to_string(parid_gbfit));
             };
    
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
              log << "Initial global best fit: " << boost::format("        %1.10e") %gbfit << endl;
              log << "flocky optimization started!" << endl;
              log.close();
           };
    
              //newSwarm.printdisp(newSwarm, 0, cycle, 1);
              newSwarm.printpos(newSwarm, 0, cycle, 1);
              if (uq == true) {
                newSwarm.printUQFF(newSwarm, 0, cycle, 1);
                newSwarm.printUQQoI(newSwarm, 0, cycle, 1);
              };
            
    
            // ---- free the sub-communicators ------ //
            for (const int& swarmcore : swarmcores) {
              if (core == swarmcore) {
                  newcomm = &ACTIVESWARM;
                  MPI_Comm_free(newcomm);
              };
            };
            for (const int& reaxffcore : reaxffcores) {
              if (core == reaxffcore) {
                  newcomm = &PASSIVESWARM;
                  MPI_Comm_free(newcomm);
              };
            };
            // -------------------------------------- //
      }else {
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
    
         if (core == cpuid_gbfit) {
            write_ffield_gbest(cpuid_gbfit, cycle, 0, parid_gbfit);
            boost::filesystem::remove("CPU." + std::to_string(core) + "/ffield.tmp."+std::to_string(cycle)+"."+"0."+std::to_string(parid_gbfit));
         };
    
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
           log << "Initial global best fit: " << boost::format("        %1.10e") %gbfit << endl;
           log << "flocky optimization started!" << endl;
           log.close();
         };
    
         //newSwarm.printdisp(newSwarm, 0, cycle, 1);
         newSwarm.printpos(newSwarm, 0, cycle, 1);
         if (uq == true) {
           newSwarm.printUQFF(newSwarm, 0, cycle, 1);
           newSwarm.printUQQoI(newSwarm, 0, cycle, 1);
         };
      };
    
    if (core == 0 && verbose == true) {
       cout << "core " << core << " left Populate!" << endl;
    };
};
#endif

#ifndef WITH_MPI
if (verbose == true) {
   cout << "Swarm entered Populate!" << endl;
};

  boost::filesystem::path pwd(boost::filesystem::current_path());
  boost::filesystem::ofstream log("log.flocky", ofstream::app);
  log << "\n";
  log << "Swarm generation started. Please wait." << endl;
  log.close();

  funceval = 0; // clear counter for cycles

  // find inverse parameters and store them in inversep
  if (regular == 1 || regular == 2) {
     get_inversep();
  };

  // ---------------------------------------------- //
  //     POPULATE: MAIN LOOP OVER SWARM MEMBERS
  // ---------------------------------------------- //

  for (int p = 0; p < NumP; p++) {
    string parID = std::to_string(p);
    string str_cycle = std::to_string(cycle);
    
    Par NewPar;
    newSwarm.AddPar(NewPar);

    // initial gbfit is INF for all processes and all swarm members
    gbfit = numeric_limits < double > ::infinity();
    // initial gbpos is 0.0 for all processes and all swarm members
    gbpos.clear();
    for (int m=0; m < dim; m++){
      gbpos.push_back(0.0);
    };

    // cp geo to geo.parID so each particle works with its own geo file
    boost::filesystem::copy_file("geo","geo." + str_cycle + "." + parID,
        boost::filesystem::copy_option::overwrite_if_exists);

    if (core == 0) {
      // If contff == y, then take force field's current values for the position of particle 0 (others are random)
      if (contff == true) {
          vector < double > ffpos;
          double tempstand;
          ffpos.clear();
          for (int m = 0; m < dim; m++) {
               ffpos.push_back(stod(newSwarm.GetPar(0).ffieldmat.at(newSwarm.GetPar(0).ffline.at(m)).at(newSwarm.GetPar(0).ffcol.at(m))));
          };
          newSwarm.GetPar(0).set_pos(ffpos);
          // reset parameters if outside domain before standardization possible
          newSwarm.GetPar(0).check_bounds_contff();
          // standardize positions
          for (int m = 0; m < dim; m++) {
              newSwarm.GetPar(0).pos.at(m) = ( newSwarm.GetPar(0).pos.at(m) - newSwarm.GetPar(0).mindomain.at(m) ) / ( newSwarm.GetPar(0).maxdomain.at(m) - newSwarm.GetPar(0).mindomain.at(m)  );
          };
      };
         contff = false;
    };
      
    
    // evaluate fitness and set bfit = curfit
    vector <double> numgrad;
    newSwarm.GetPar(p).state.cycle = cycle;
    iter = 0;
    newSwarm.GetPar(p).state.iter = iter;
    newSwarm.GetPar(p).state.parid = p;

    // evaluate fitness
    vector <double> numgrad;

    // dropout
    if (regular == 3) {
       newSwarm.GetPar(p).dropout(dropvalue);
    };

    // local minimization
    if (localmin == 1 || localmin == 2) {
       newSwarm.GetPar(p).iterate();
    };

    newSwarm.GetPar(p).set_fitness(newSwarm.GetPar(p).eval_fitness(newSwarm.GetPar(p).get_pos_vec(), this));
    newSwarm.GetPar(p).set_bfit(newSwarm.GetPar(p).get_fitness());

    if (newSwarm.GetPar(p).get_fitness() < gbfit) {
      gbfit = newSwarm.GetPar(p).get_fitness();
      gbpos.clear();
      gbpos = newSwarm.GetPar(p).get_pos_vec();
      write_ffield_gbest(0, cycle, 0, p);
    };
    // cleaning ffield.tmp.* files after write_ffield_gbest already
    // copied the correct ffield.tmp.* file as the ffield.gbest.*.0.*
    boost::filesystem::remove("ffield.tmp." + str_cycle + ".0." + "." + parID);
  }; // done loop over members

if (verbose == true) {
   cout << "Swarm done loop on members!" << endl;
};

  //initial_disp = newSwarm.get_disp(newSwarm);

  // clean up any old files belonging to previous job
  if ( boost::filesystem::exists( "opti_log.out." + std::to_string(cycle)) ){
    boost::filesystem::remove( "opti_log.out." + std::to_string(cycle) );
  };
  if ( boost::filesystem::exists( "disp_log.out." + std::to_string(cycle)) ){
    boost::filesystem::remove( "disp_log.out." + std::to_string(cycle) );
  };

  newSwarm.printopt(newSwarm, 0, cycle, 1);
  boost::filesystem::ofstream log2("log.flocky", ofstream::app);
  log2 << "\nSwarm generation completed." << endl;
  log2 << "Initial global best fit: " << boost::format("        %1.10e") %gbfit << endl;
  log2 << "flocky optimization started!" << endl;
  log2.close();

  //newSwarm.printdisp(newSwarm, 0, cycle, 1);
  newSwarm.printpos(newSwarm, 0, cycle, 1);
  if (uq == true){
    newSwarm.printUQFF(newSwarm, 0, cycle, 1);
    newSwarm.printUQQoI(newSwarm, 0, cycle, 1);
  };

if (verbose == true) {
   cout << "Swarm left Populate!" << endl;
};
#endif
};



void Swarm::Propagate(Swarm & newSwarm, int cycle) {
#ifdef WITH_MPI
if (verbose == true) {
   cout << "CPU: " << core << " entered Propagate()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered Propagate()" << endl;
};
#endif
#ifdef WITH_MPI
if (ensembleave == false && preppath == 0) {
   if (core == 0 && verbose == true) {
      cout << "CPU: " << core << " entered Propagate!" << endl;
   };
   
    // define a new sub-communicator for reaxffcores and swarmcores //
       MPI_Comm ACTIVESWARM;  // swarmcores
       MPI_Comm PASSIVESWARM; // reaxffcores
       MPI_Comm *newcomm;
       int color;
   
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
           color = 444;
           newcomm = &ACTIVESWARM;
       };
       if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
           color = 333;
           newcomm = &PASSIVESWARM;
       };
       MPI_Comm_split( MPI_COMM_WORLD, color, core, newcomm );
   
       if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
            MPI_Comm_rank( ACTIVESWARM, &mycore_swarmcore);
            MPI_Comm_size( ACTIVESWARM, &size_swarmcores);
       };
   
       if (find(reaxffcores.begin(), reaxffcores.end(), core) != reaxffcores.end()) {
            MPI_Comm_rank( PASSIVESWARM, &mycore_reaxffcore);
            MPI_Comm_size( PASSIVESWARM, &size_reaxffcores);
       };
    // ------------ end definition of new communicator ---------- //
   
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
          // if we parallelize the training set, do the following only among swarmcores
          if (ptrainset > 1) {
             if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
                 // Update velocities and positions
                 newSwarm.GetPar(p).update_vel(inertiafac, confac, gbpos, iter);
                 if (newSwarm.GetPar(p).fails > faili) {
                   newSwarm.GetPar(p).update_pos_levy(gbpos, iter, inertiafac);
                   newSwarm.GetPar(p).fails = 0;
                 } else {
                   newSwarm.GetPar(p).update_pos();
                 };
             };
             newSwarm.GetPar(p).state.cycle = cycle;
             newSwarm.GetPar(p).state.iter = iter;
             newSwarm.GetPar(p).state.parid = p;
   
             // evaluate fitness. if doing localmin with ptrainset > 1, the generation of trainset subsets
             // is performed inside eval_fitness. If not doing localmin, generation of trainset subsets is performed here.
   
             // dropout
             if (regular == 3) {
                newSwarm.GetPar(p).dropout(dropvalue);
             };
   
             if (localmin == 1 || localmin == 2) {
                newSwarm.GetPar(p).iterate();
             }else{
                newSwarm.GetPar(p).set_fitness(newSwarm.GetPar(p).eval_fitness(newSwarm.GetPar(p).get_pos_vec(), this));
             };
   
             if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
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
                    //write_ffield_gbest(core, cycle, iter, p);
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
             };
          }else {
                 // Update velocities and positions
                 newSwarm.GetPar(p).update_vel(inertiafac, confac, gbpos, iter);
                 if (newSwarm.GetPar(p).fails > faili) {
                   newSwarm.GetPar(p).update_pos_levy(gbpos, iter, inertiafac);
                   newSwarm.GetPar(p).fails = 0;
                 } else {
                   newSwarm.GetPar(p).update_pos();
                 };
   
                 newSwarm.GetPar(p).state.cycle = cycle;
                 newSwarm.GetPar(p).state.iter = iter;
                 newSwarm.GetPar(p).state.parid = p;
   
                // dropout
                if (regular == 3) {
                   newSwarm.GetPar(p).dropout(dropvalue);
                };
   
                 if (localmin == 1 || localmin == 2) {
                    newSwarm.GetPar(p).iterate();
                 } else {
                    vector <double> numgrad;
                    newSwarm.GetPar(p).set_fitness(newSwarm.GetPar(p).eval_fitness(newSwarm.GetPar(p).get_pos_vec(), this));
                 };
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
                     //write_ffield_gbest(core, cycle, iter, p);
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
          };
       }; // done loop over swarm members
   
       // if we parallelize the training set, do the following only among swarmcores
       if (ptrainset > 1) {
           //newSwarm.get_com(newSwarm);
   
           // pair struct to hold the global best fitness across processes and its core rank
           struct {
             double tmp_fit;
             int tmp_cpu;
           } min_vals_in[1], min_vals_out[1];
   
           // store current fit on each process
           min_vals_in[0].tmp_fit = gbfit;
           // store core id of that current process
           min_vals_in[0].tmp_cpu = mycore_swarmcore;
           // get global best fitness *across processes* and corresponding core rank and store them in min_vals_out
           if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              MPI_Allreduce( & min_vals_in, & min_vals_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, ACTIVESWARM);
           };
   
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
           if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              MPI_Allreduce( & min_vals_in2, & min_vals_out2, 1, MPI_DOUBLE_INT, MPI_MINLOC, ACTIVESWARM);
           };
           // parid of the gbfit
           parid_gbfit = min_vals_out2[0].tmp_parid;
   
           // global best fitness across all processes
           gbfit = min_vals_out[0].tmp_fit;
           // core rank the above fitness came from
           cpuid_gbfit = min_vals_out[0].tmp_cpu;
           // broadcast contents of gbpos vector from rank cpuid_gbfit
           if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              MPI_Bcast(gbpos.data(), gbpos.size(), MPI_DOUBLE, cpuid_gbfit, ACTIVESWARM);
           };
   
           // now convert cpuid_gbfit value from ACTIVESWARM to MPI_COMM_WORLD
           if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
              cpuid_gbfit = swarmcores.at(cpuid_gbfit);
           };
   
           // do write_ffield_gbest
           if (core == cpuid_gbfit && find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end() && gbfitfound == true) {
              //cout << "CPU: " << core << " (cpuid_gbfit) writes ffield.gbest: " << "ffield.tmp." + to_string(cycle)+"."+to_string(iter)+"."+to_string(parid_gbfit) << endl;
              write_ffield_gbest(cpuid_gbfit, cycle, iter, parid_gbfit);
              //boost::filesystem::remove(pwd.string() + "/CPU." + std::to_string(core) + "/ffield.tmp."+std::to_string(cycle)+"."+std::to_string(iter)+"."+std::to_string(parid_gbfit));
           };
           boost::filesystem::remove("CPU." + std::to_string(core) + "/ffield.tmp."+std::to_string(cycle)+"."+std::to_string(iter)+"."+std::to_string(parid_gbfit));
   
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
              if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
                 MPI_Bcast( &firstovfit, 1, MPI_C_BOOL, cpuid_gbfit, ACTIVESWARM);
              };
             //MPI_Bcast( &firstovfit, 1, MPI_C_BOOL, cpuid_gbfit, ACTIVESWARM);
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
   
           // reset swarm randomly if gbfit == INF for ninf consecutive iterations
           if (gbfit == numeric_limits < double > ::infinity()) {
              ninf = ninf + 1;
              if (ninf == 10) {
                  std::uniform_real_distribution < double > dist2(0.0,1.0);
                  for (int p = 0; p < NumP; p++) {
                     for (int i = 0; i < dim; i++) {
                        newSwarm.GetPar(p).set_posdim(i, dist2(generator));
                        newSwarm.GetPar(p).set_veldim(i, 0.5*dist2(generator) - newSwarm.GetPar(p).get_pos(i));
                        newSwarm.GetPar(p).bpos.at(i) = newSwarm.GetPar(p).pos.at(i);
                     };
                     ninf = 0;
                  };
              };
           };
   
       }else { // done if ptrainset > 1
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
   
           if (core == cpuid_gbfit && gbfitfound == true) {
              write_ffield_gbest(cpuid_gbfit, cycle, iter, parid_gbfit);
              //boost::filesystem::remove(pwd.string() + "/CPU." + std::to_string(core) + "/ffield.tmp."+std::to_string(cycle)+"."+std::to_string(iter)+"."+std::to_string(parid_gbfit));
           };
           boost::filesystem::remove("CPU." + std::to_string(core) + "/ffield.tmp."+std::to_string(cycle)+"."+std::to_string(iter)+"."+std::to_string(parid_gbfit));
   
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
   
           // reset swarm randomly if gbfit == INF for ninf consecutive iterations
           if (gbfit == numeric_limits < double > ::infinity()) {
              ninf = ninf + 1;
              if (ninf == 10) {
                  std::uniform_real_distribution < double > dist2(0.0,1.0);
                  for (int p = 0; p < NumP; p++) {
                     for (int i = 0; i < dim; i++) {
                        newSwarm.GetPar(p).set_posdim(i, dist2(generator));
                        newSwarm.GetPar(p).set_veldim(i, 0.5*dist2(generator) - newSwarm.GetPar(p).get_pos(i));
                        newSwarm.GetPar(p).bpos.at(i) = newSwarm.GetPar(p).pos.at(i);
                     };
                     ninf = 0;
                  };
              };
           };
   
       };
     }; // done loop over iterations
   
     if (ptrainset == 1) {
         MPI_Allreduce(& funceval, & funceval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
     }else {
         if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
             MPI_Allreduce(& funceval, & funceval, 1, MPI_INT, MPI_SUM, ACTIVESWARM);
         };
     };
     MPI_Barrier(MPI_COMM_WORLD);
   
     if (core == 0) {
       double temphys;
       boost::filesystem::ofstream log("log.flocky", ofstream::app);
       log << "\n\nTraining completed successfuly!\n";
       log << "Total ReaxFF calls: " << funceval << endl;
       log << "\nGlobal best ReaxFF fit: " << boost::format("        %1.10e") %newSwarm.get_gbfit() << endl;
       log << "Global best ReaxFF parameters:" << endl;
       log << "[ ";
       for (int m = 0; m < dim; m++) {
         temphys = newSwarm.get_gbpos().at(m)*(newSwarm.GetPar(0).maxdomain.at(m) - newSwarm.GetPar(0).mindomain.at(m)) + newSwarm.GetPar(0).mindomain.at(m);
         // physical gbest position
         //log << boost::format("%8.4f") %newSwarm.get_gbpos().at(m) << " ";
         // standardized gbest position
         log << boost::format("%8.4f") %temphys << " ";
       };
       log << "]\n" << endl;
       log.close();
     };
   
     // ---- free the sub-communicators ------ //
     for (const int& swarmcore : swarmcores) {
       if (core == swarmcore) {
           newcomm = &ACTIVESWARM;
           MPI_Comm_free(newcomm);
       };
     };
     for (const int& reaxffcore : reaxffcores) {
       if (core == reaxffcore) {
           newcomm = &PASSIVESWARM;
           MPI_Comm_free(newcomm);
       };
     };
   if (core == 0 && verbose == true) {
      cout << "core " << core << " left Populate!" << endl;
   };
};
#endif

#ifndef WITH_MPI
if (verbose == true) {
   cout << "Swarm entered Propagate!" << endl;
};

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

      newSwarm.GetPar(p).state.cycle = cycle;
      newSwarm.GetPar(p).state.iter = iter;
      newSwarm.GetPar(p).state.parid = p;

      // dropout
      if (regular == 3) {
         newSwarm.GetPar(p).dropout(dropvalue);
      };

      if (localmin == 1 || localmin == 2) {
         newSwarm.GetPar(p).iterate();
      } else {
          vector <double> numgrad;
          newSwarm.GetPar(p).set_fitness(newSwarm.GetPar(p).eval_fitness(newSwarm.GetPar(p).get_pos_vec(), this));
      };

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

  double temphys;
  boost::filesystem::ofstream log("log.flocky", ofstream::app);
  log << "\n\nTraining completed successfuly!\n";
  log << "Total ReaxFF calls: " << funceval << endl;
  log << "\nGlobal best ReaxFF fit: " << boost::format("        %1.10e") %newSwarm.get_gbfit() << endl;
  log << "Global best ReaxFF parameters:" << endl;
  log << "[ "; 
  for (int m = 0; m < dim; m++) {
    // physical gbest position
    //log << boost::format("%8.4f") %newSwarm.get_gbpos().at(m) << " ";
    // standardized gbest position
    temphys = newSwarm.get_gbpos().at(m)*(newSwarm.GetPar(0).maxdomain.at(m) - newSwarm.GetPar(0).mindomain.at(m)) + newSwarm.GetPar(0).mindomain.at(m);
    log << boost::format("%8.4f") %temphys << " ";
  };
  log << "]" << endl;
  log.close();

if (verbose == true) {
   cout << "Swarm left Propagate!" << endl;
};
#endif
};

void Swarm::write_ffield_gbest(int core, int cycle, int iter, int par) {
#ifdef WITH_MPI
if (verbose == true) {
   cout << "core " << core << " entered write_ffield_gbest()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered write_ffield_gbest()" << endl;
};
#endif

  // cp current ffield to be the global best and current analysis files to global best analysis files
#ifndef WITH_MPI

  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(par);

  boost::filesystem::ifstream test_ffieldtmp("ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID);
  if (!test_ffieldtmp.fail()) {
    boost::filesystem::copy_file("ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID,
    "ffield.gbest." + str_cycle + "." + str_iter + "." + str_parID, boost::filesystem::copy_option::overwrite_if_exists);
  };
  test_ffieldtmp.close();

  boost::filesystem::ifstream test_fort99("fort.99");
  if (!test_fort99.fail()) {
    boost::filesystem::copy_file("fort.99", "results.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  };
  test_fort99.close();

  // clean temp files
  //boost::filesystem::remove(pwd.string() + "/ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID);

#endif
#ifdef WITH_MPI

  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_cycle = std::to_string(cycle);
  string str_core = std::to_string(core);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(par);

  boost::filesystem::ifstream test_ffieldtmp("CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID);
  if (!test_ffieldtmp.fail()) {
      boost::filesystem::copy_file("CPU." + str_core + "/ffield.tmp."+str_cycle+"."+str_iter+"."+str_parID,
        "ffield.gbest." + str_cycle + "." + str_iter+"."+str_parID, boost::filesystem::copy_option::overwrite_if_exists);
  };
  test_ffieldtmp.close();

  boost::filesystem::ifstream test_fort99("CPU." + str_core + "/fort.99");
  if (!test_fort99.fail()) {
     boost::filesystem::copy_file("CPU." + str_core + "/fort.99",
       "results.out." + str_cycle, boost::filesystem::copy_option::overwrite_if_exists);
  };
  test_fort99.close();

  // clean temp files
  //boost::filesystem::remove(pwd.string() + "/CPU." + str_core + "/ffield.tmp." + str_cycle + "." + str_iter + "." + str_parID);

#endif
};

void Swarm::detovfit(Swarm &newSwarm, int cpuid_gbfit, int cycle, int iter, int parid_gbfit) {
#ifdef WITH_MPI
if (verbose == true) {
   cout << "core " << core << " entered detovfit()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered detovfit()" << endl;
};
#endif
#ifdef WITH_MPI
  boost::filesystem::path pwd(boost::filesystem::current_path());
  string str_core = std::to_string(cpuid_gbfit);
  string str_cycle = std::to_string(cycle);
  string str_iter = std::to_string(iter);
  string str_parID = std::to_string(parid_gbfit);
  // prepare mandatory files for reaxff in 'testovfit' directory
  if (firstovfit == true) {
      boost::filesystem::create_directory("testovfit");
      boost::filesystem::copy_file("ffield.initial." + str_cycle,
        "testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    boost::filesystem::copy_file("ffield.gbest." + str_cycle + "." + str_iter + "." + str_parID,
      "testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file("geo.val", "testovfit/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

  if (fixcharges == true){
    boost::filesystem::copy_file("charges_val","testovfit/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file("control", "testovfit/control",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::copy_file("valset.in", "testovfit/trainset.in",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::ofstream iopt_file("testovfit/fort.20");
  iopt_file << "0";
  iopt_file.close();
  ofstream outfile35("testovfit/fort.35");
  outfile35 << "23434.1" << endl;
  outfile35.close();

  // cd to testovfit directory
  string old_path = pwd.string();
  boost::filesystem::path p("testovfit");
  boost::filesystem::current_path(p);

  // execute tapreaxff
  boost::filesystem::copy_file("tapreaxff", "testovfit/tapreaxff",
  boost::filesystem::copy_option::overwrite_if_exists);
  //arguments for tapreaxff, will run: tapreaxff                                                                                                                                  
  char *args[3] = { "./tapreaxff", "", NULL} ;
  pid_t c_pid, pid;
  int status;

  /* create a child process */
  c_pid = fork();

  if (c_pid == 0){
    /* CHILD */
    char *filename  = "run.log";
    int outfile = open(filename, O_CREAT | O_WRONLY, S_IRWXU);
    if (outfile == -1){
          fprintf(stderr, "Error: failed to create file %s\n", filename);
    }else{
          /* redirect the standard output from this process to the file. */
          if(dup2(outfile, STDOUT_FILENO) != STDOUT_FILENO){
            fprintf(stderr, "Error: failed to redirect standard output\n");
          }
          /* redirect tapreaxff stdout to outfile*/
          dup2 (outfile, STDOUT_FILENO);
          /* redirect tapreaxff stderr to /dev/null */
          dup2(open("/dev/null", 0), 2);

        /* printf("Child: executing args\n"); */
        // execute args                                                                                                                                                               
        execvp( args[0], args);
        // only get here if exec failed                                                                                                                                             
        perror("execve failed");
        wait(&status);
    };
  }else if (c_pid > 0){
    /* PARENT */

    if( (pid = wait(&status)) < 0){
      perror("wait");
      _exit(1);
    };
    //printf("Parent: finished\n");

  }else{
    perror("fork failed");
    _exit(1);
  };

  if (WIFSIGNALED (status)) {
    cout << "Error: tapreaxff exited abnormaly on CPU:" << core << "\n";
    MPI_Abort(MPI_COMM_WORLD,4); 
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
  if (str.at(0) == '*' || !isdigit(str.at(0))) {
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
      boost::filesystem::copy_file("ffield.initial." + str_cycle,
        "testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    boost::filesystem::copy_file("ffield.gbest." + str_cycle + "." + str_iter + "." + str_parID,
      "testovfit/fort.4", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file("geo.val",
    "testovfit/fort.3", boost::filesystem::copy_option::overwrite_if_exists);

  if (fixcharges == true){
    boost::filesystem::copy_file("charges_val",
      "testovfit/fort.26", boost::filesystem::copy_option::overwrite_if_exists);
  };

  boost::filesystem::copy_file("control", "testovfit/control",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::copy_file("valset.in", "testovfit/trainset.in",
    boost::filesystem::copy_option::overwrite_if_exists);

  boost::filesystem::ofstream iopt_file("testovfit/fort.20");
  iopt_file << "0";
  iopt_file.close();
  ofstream outfile35("testovfit/fort.35");
  outfile35 << "23434.1" << endl;
  outfile35.close();

  // cd to testovfit directory
  string old_path = pwd.string();
  boost::filesystem::path p("testovfit");
  boost::filesystem::current_path(p);

  // execute tapreaxff
  boost::filesystem::copy_file("tapreaxff", "testovfit/tapreaxff",
  boost::filesystem::copy_option::overwrite_if_exists);
  //arguments for tapreaxff, will run: tapreaxff                                                                                                                                  
  char *args[3] = { "./tapreaxff", "", NULL} ;
  pid_t c_pid, pid;
  int status;

  /* create a child process */
  c_pid = fork();

  if (c_pid == 0){
    /* CHILD */
    char *filename  = "run.log";
    int outfile = open(filename, O_CREAT | O_WRONLY, S_IRWXU);
    if (outfile == -1){
          fprintf(stderr, "Error: failed to create file %s\n", filename);
    }else{
          /* redirect the standard output from this process to the file. */
          if(dup2(outfile, STDOUT_FILENO) != STDOUT_FILENO){
            fprintf(stderr, "Error: failed to redirect standard output\n");
          }
          /* redirect tapreaxff stdout to outfile*/
          dup2 (outfile, STDOUT_FILENO);
          /* redirect tapreaxff stderr to /dev/null */
          dup2(open("/dev/null", 0), 2);

        /* printf("Child: executing args\n"); */
        // execute args                                                                                                                                                               
        execvp( args[0], args);
        // only get here if exec failed                                                                                                                                             
        perror("execve failed");
        wait(&status);
    };
  }else if (c_pid > 0){
    /* PARENT */

    if( (pid = wait(&status)) < 0){
      perror("wait");
      _exit(1);
    };
    //printf("Parent: finished\n");

  }else{
    perror("fork failed");
    _exit(1);
  };

  if (WIFSIGNALED (status)) {
    cout << "Error: tapreaxff exited abnormaly on CPU:" << core << "\n";
    exit(EXIT_FAILURE);
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
  if (str.at(0) == '*' || !isdigit(str.at(0))) {
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
if (core == 0){
 ofstream outfiledisp("disp_log.out." + std::to_string(cycle), ofstream::app);
 if (mod(iter, fr) == 0.0) {
   stringstream ss;
   ss << boost::format("%5i %15.4f") %iter %(newSwarm.get_disp(newSwarm)/initial_disp);
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
    ss << boost::format("%5i         %1.10e") %iter %newSwarm.get_gbfit(); 
    outfileopt << ss.str();
    outfileopt << endl;
    
  };
  outfileopt.close();
};


void Swarm::printUQQoI(Swarm & newSwarm, int iter, int cycle, int fr) {
#ifdef WITH_MPI
if (verbose == true) {
   cout << "core " << core << " entered printUQQoI()" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered printUQQoI" << endl;
};
#endif
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
  std::regex e7 ("Force on atom:");

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
          if (regex_search (line,m,e6)){
            boost::trim(line);
            using boost::is_any_of;
            boost::split(values, line, is_any_of("\t "),boost::token_compress_on);
            // UQAFORC.* file will be generated for every cycle containing errors of all members
            ofstream outfileUQAFORC("UQAFORC." + std::to_string(cycle), ofstream::app);
            stringstream ss1;
            // determine position of error backwards
            int tmpcell = values.size() - 2;
            ss1 << boost::format("%25.4f") %sqrt(stod(values.at(tmpcell)));
            outfileUQAFORC << ss1.str() << endl;
            ss1.str("");
            values.clear();
            outfileUQAFORC.close();
            continue;
          };
        }; // done reading all lines
      };
      myfile.close();
    }; // done loop on members
  };
};

void Swarm::printUQFF(Swarm & newSwarm, int iter, int cycle, int fr) {
#ifdef WITH_MPI
if (verbose == true) {
   cout << "core " << core << " entered printUQFF" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered printUQFF" << endl;
};
#endif
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
#ifdef WITH_MPI
if (verbose == true) {
   cout << "core " << core << " entered printpos" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "entered printpos" << endl;
};
#endif
  boost::filesystem::path pwd(boost::filesystem::current_path());
  double phystempos;
#ifdef WITH_MPI
if (core == 0 && verbose == true) {
   cout << "core " << core << " entered printpos!" << endl;
};

  if (ptrainset > 1) {
     //if (find(swarmcores.begin(), swarmcores.end(), core) != swarmcores.end()) {
        string str_core = std::to_string(core);
        boost::filesystem::create_directory("CPU." + str_core + "/pos");
        const char * path = "/pos/pos_log.out.";
        ofstream outfilepos("CPU." + str_core + path + std::to_string(cycle), ofstream::app);
        if (mod(iter, fr) == 0.0) {
          outfilepos << "#iter: " << iter << endl;
          stringstream ss1;
          stringstream ss2;
          stringstream ss3;
          for (int p = 0; p < NumP; p++) {
            phystempos = 0.0;
            ss1 << boost::format("%3i") %p;
            outfilepos << ss1.str();
            ss1.str("");
            for (int m = 0; m < dim; m++) {
              // transform back to physical positions
              phystempos = newSwarm.GetPar(p).get_pos(m)*(newSwarm.GetPar(p).maxdomain.at(m) - newSwarm.GetPar(p).mindomain.at(m)) + newSwarm.GetPar(p).mindomain.at(m);
              // stream standardized position
              //ss2 << boost::format("%18.4f") %newSwarm.GetPar(p).get_pos(m); 
              // stream transformed to physical position
              ss2 << boost::format("%18.4f") %phystempos;
              outfilepos << ss2.str();
              ss2.str("");
            };
            ss3 << boost::format("        %1.10e") %newSwarm.GetPar(p).get_fitness();
            outfilepos << ss3.str();
            ss3.str("");
            outfilepos << endl;
          };
        };
        outfilepos.close();
     //};
  }else {
     string str_core = std::to_string(core);
     boost::filesystem::create_directory("CPU." + str_core + "/pos");
     const char * path = "/pos/pos_log.out.";
     ofstream outfilepos("CPU." + str_core + path + std::to_string(cycle), ofstream::app);
     if (mod(iter, fr) == 0.0) {
       outfilepos << "#iter: " << iter << endl;
       stringstream ss1;
       stringstream ss2;
       stringstream ss3;
       for (int p = 0; p < NumP; p++) {
         phystempos = 0.0;
         ss1 << boost::format("%3i") %p;
         outfilepos << ss1.str();
         ss1.str("");
         for (int m = 0; m < dim; m++) {
           // transform back to physical positions
           phystempos = newSwarm.GetPar(p).get_pos(m)*(newSwarm.GetPar(p).maxdomain.at(m) - newSwarm.GetPar(p).mindomain.at(m)) + newSwarm.GetPar(p).mindomain.at(m);
           // stream standardized position
           //ss2 << boost::format("%18.4f") %newSwarm.GetPar(p).get_pos(m); 
           // stream transformed to physical position
           ss2 << boost::format("%18.4f") %phystempos;
           outfilepos << ss2.str();
           ss2.str("");
         };
         ss3 << boost::format("        %1.10e") %newSwarm.GetPar(p).get_fitness();
         outfilepos << ss3.str();
         ss3.str("");
         outfilepos << endl;
       };
     };
     outfilepos.close();
  };

if (core == 0 && verbose == true) {
   cout << "core " << core << " entered printpos!" << endl;
};
#endif
#ifndef WITH_MPI
if (verbose == true) {
   cout << "Swarm entered printpos!" << endl;
};

  boost::filesystem::create_directory("pos");
  const char * path = "pos/pos_log.out.";
  ofstream outfilepos(path + std::to_string(cycle), ofstream::app);
  if (mod(iter, fr) == 0.0) {
    outfilepos << "#iter: " << iter << endl;
    stringstream ss1;
    stringstream ss2;
    stringstream ss3;
    for (int p = 0; p < NumP; p++) {
      phystempos = 0.0;
      ss1 << boost::format("%3i") %p;
      outfilepos << ss1.str();
      ss1.str("");
      for (int m = 0; m < dim; m++) {
        // transform back to physical positions
        phystempos = newSwarm.GetPar(p).get_pos(m)*(newSwarm.GetPar(p).maxdomain.at(m) - newSwarm.GetPar(p).mindomain.at(m)) + newSwarm.GetPar(p).mindomain.at(m);
        // stream standardized position
        //ss2 << boost::format("%18.4f") %newSwarm.GetPar(p).get_pos(m); 
        // stream transformed to physical position
        ss2 << boost::format("%18.4f") %phystempos;
        outfilepos << ss2.str();
        ss2.str("");
      };
      ss3 << boost::format("        %1.10e") %newSwarm.GetPar(p).get_fitness();
      outfilepos << ss3.str();
      ss3.str("");
      outfilepos << endl;
    };
  };
  outfilepos.close();
#endif
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
