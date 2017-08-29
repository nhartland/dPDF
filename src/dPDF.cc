// n.p.hartland@ed.ac.uk  03/12

#include "LHAPDF/LHAPDF.h"

// Configuration file
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/chisquared.h"

#include "fastaddchi2.h"
#include "deuteronset.h"
#include "ns_network.h"
#include "cmaes.h"
#include "filter.h"
#include "colour.h"
#include "proton.h"

using namespace Colour;
using namespace std;

void gsl_handler (const char * msg, const char * source, int line, int)
{ std::cerr << "gsl: " << source<<":"<<line<<": " <<msg <<std::endl;}

int main(int argc, char* argv[]) {

  // Set LHAPDF verbosity
  LHAPDF::setVerbosity(0);
  NNPDF::SetVerbosity(0);

  // Set GSL error handler
  gsl_set_error_handler (&gsl_handler);

  if (argc < 3)
  {
    cerr << "Usage: dPDF <configuration file> <replica number>"<<std::endl;
    exit(-1);
  }
  // Read configuration file
  libconfig::Config dPDFconfig;
  dPDFconfig.readFile(argv[1]);

  const int replica = atoi(argv[2]);
  libconfig::Setting& fitsetting = dPDFconfig.lookup("fit");
  libconfig::Setting & replica_setting = fitsetting.add ("replica", libconfig::Setting::TypeInt); replica_setting = replica;

  const std::string fitname = dPDFconfig.lookup("fit.name");
  const std::string base_path = "./res/"+fitname;
  if (replica == 0)
  {
    mkdir(base_path.c_str(), 0777);
    mkdir((base_path+"/par").c_str(), 0777);
    mkdir((base_path+"/erf").c_str(), 0777);
    mkdir((base_path+"/dat").c_str(), 0777);
    mkdir((base_path+"/dat/systypes").c_str(), 0777);
  }

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  // Initialise RNG
  const int method = dPDFconfig.lookup("rng.mode");
  const int globalseed   = dPDFconfig.lookup("rng.seed");
  NNPDF::RandomGenerator::InitRNG(method,globalseed);

    // Setting the correct replica based seed
  int fitseed = 0;
  for (int i = 0; i < replica; i++)
    fitseed = NNPDF::RandomGenerator::GetRNG()->GetRandomInt();
  NNPDF::RandomGenerator::InitRNG(method,fitseed);

  cout              << "           RNG Init - Mode:"<<method<<"  Seed:"<<globalseed<<"  ReplicaSeed: "<<fitseed<<  endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  // Initialise and filter datasets
  std::vector<NNPDF::Experiment> experimental_data;
  std::vector<NNPDF::Experiment> training_data;
  std::vector<NNPDF::Experiment> validation_data;
  ReadData(dPDFconfig, experimental_data);
  InitData(dPDFconfig, experimental_data, training_data, validation_data);

  double nData_trn = 0;
  for (size_t i=0; i<training_data.size(); i++)
    nData_trn += training_data[i].GetNData();

  double nData_val = 0;
  for (size_t i=0; i<validation_data.size(); i++)
    nData_val += validation_data[i].GetNData();

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  cout              << "                            deuteron fitter                            "<<             endl;
  cout              << "                     Fit: "<<fitname<<"     Replica: "<<replica<<"   "<<             endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  const int lambda = dPDFconfig.lookup("fit.lambda");
  const int nparam = NostateMLP::get_nparam(pdf_architecture);
  int ite_look_back = -1; double erf_look_back = std::numeric_limits<double>::infinity();

  // Initialise minimiser
  CMAESMinimizer min(nparam, lambda, dPDFconfig.lookup("fit.sigma"));
  gsl_vector* deuteron_search_centre = gsl_vector_alloc(nparam); // Centre of search distribution
  gsl_vector* deuteron_look_back     = gsl_vector_alloc(nparam); // Look-back best fit
  min.NormVect(deuteron_search_centre);

  // Initialise proton parametrisation and cost functions
  LHAPDFSet pPDF(dPDFconfig.lookup("fit.proton"), replica + 1, lambda);
  LHAPDFSet bpPDF(dPDFconfig.lookup("fit.proton"), replica + 1, 1);
  ErfComputer training_cost(pPDF, training_data);
  ErfComputer training_cost_centre(bpPDF, training_data);
  ErfComputer validation_cost_centre(bpPDF, validation_data);

  // Error function output
  std::stringstream erf_filename;
  erf_filename << base_path<< "/erf/replica_"<<replica<<".dat";
  ofstream erf_file; erf_file.open(erf_filename.str());

  const int ngen = dPDFconfig.lookup("fit.ngen");
  for (int i=0; i< ngen; i++)
  {
    std::cout << "Iteration: "<<i <<" / " <<ngen <<std::endl;
    min.Iterate(deuteron_search_centre, &training_cost);

    const double erf_training   = training_cost_centre(deuteron_search_centre)/nData_trn;
    const double erf_validation = validation_cost_centre(deuteron_search_centre)/nData_val;
    erf_file << i << "  " <<  erf_training << "  "<< erf_validation<<std::endl; 

    if (erf_validation < erf_look_back || erf_validation < erf_training )
    {
      ite_look_back = i; erf_look_back = erf_validation;
      gsl_vector_memcpy(deuteron_look_back, deuteron_search_centre);
    }
  }

  const double trnchi2 = training_cost_centre(deuteron_look_back)/nData_trn;
  const double valchi2 = validation_cost_centre(deuteron_look_back)/nData_val;
  const double glbchi2 = ErfComputer(bpPDF, experimental_data)(deuteron_look_back)/(nData_val + nData_trn);
  erf_file << ite_look_back << "  " <<  trnchi2 << "  "<< valchi2<<"  "<<glbchi2<<std::endl; 
  erf_file.close();

  // Export best-fit parameters
  FILE * f = fopen ( (base_path + "/par/parameters_" + to_string(replica)+".dat").c_str(), "w");
  gsl_vector_fwrite (f, deuteron_look_back); fclose (f);

  gsl_vector_free(deuteron_search_centre);
  gsl_vector_free(deuteron_look_back);

  exit(0);
}

