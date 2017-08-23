// n.p.hartland@ed.ac.uk  03/12

#include "LHAPDF/LHAPDF.h"

// Configuration file
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
#include <sys/stat.h>

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/chisquared.h"

#include "fastaddchi2.h"
#include "deuteronset.h"
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
    mkdir((base_path+"/prt").c_str(), 0777);
    mkdir((base_path+"/par").c_str(), 0777);
    mkdir((base_path+"/pdf").c_str(), 0777);
    mkdir((base_path+"/erf").c_str(), 0777);
    mkdir((base_path+"/thp").c_str(), 0777);
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

  // Initialise prototype parametrisation
  DeuteronSet deuteron_search_mutants(dPDFconfig, dPDFconfig.lookup("fit.lambda"));
  DeuteronSet deuteron_search_centre (dPDFconfig, 1); // Centre of search distribution
  DeuteronSet deuteron_look_back     (dPDFconfig, 1); // Best-fit PDFs
  const int lambda = deuteron_search_mutants.GetMembers();
  const int nparam = deuteron_search_mutants.GetNParameters();
  int ite_look_back = -1; double erf_look_back = std::numeric_limits<double>::infinity();

  // Initialise proton parametrisation
  LHAPDFSet pPDF(dPDFconfig.lookup("fit.proton"), replica + 1, lambda);
  LHAPDFSet bpPDF(dPDFconfig.lookup("fit.proton"), replica + 1, 1);

  // Initialise minimiser
  CMAESMinimizer min(nparam, lambda, dPDFconfig.lookup("fit.sigma"));
  min.NormVect(deuteron_search_centre.GetParameters(0));

  // Error function output
  std::stringstream erf_filename;
  erf_filename << base_path<< "/erf/replica_"<<replica<<".dat";
  ofstream erf_file; erf_file.open(erf_filename.str());

  const int ngen = dPDFconfig.lookup("fit.ngen");
  for (int i=0; i< ngen; i++)
  {
    std::cout << "Iteration: "<<i <<" / " <<ngen <<std::endl;
    min.Iterate(&pPDF, &deuteron_search_mutants, deuteron_search_centre.GetParameters(0), training_data);

    // Report chi2
    // gsl_vector_memcpy(deuteron_search_centre.GetParameters(0), deuteron_search_mutants.GetBestFit());
    deuteron_search_centre.InitPDFSet();

    const double trnchi2 = ComputeMemberChi2(&bpPDF, &deuteron_search_centre, 0, training_data)/nData_trn;
    const double valchi2 = ComputeMemberChi2(&bpPDF, &deuteron_search_centre, 0, validation_data)/nData_val;
    erf_file << i << "  " <<  trnchi2 << "  "<< valchi2<<std::endl; 

    if (valchi2 < erf_look_back)
    {
      ite_look_back = i;
      erf_look_back = valchi2;
      gsl_vector_memcpy(deuteron_look_back.GetParameters(0), deuteron_search_centre.GetParameters(0));
      deuteron_look_back.InitPDFSet();
    }
  }

  const double trnchi2 = ComputeMemberChi2(&bpPDF, &deuteron_look_back, 0, training_data)/nData_trn;
  const double valchi2 = ComputeMemberChi2(&bpPDF, &deuteron_look_back, 0, validation_data)/nData_val;
  erf_file << ite_look_back << "  " <<  trnchi2 << "  "<< valchi2<<std::endl; 
  erf_file.close();

  // Export PDF values
  std::stringstream filename;
  filename << base_path<< "/pdf/replica_"<<replica<<".dat";
  ofstream outfile; outfile.open(filename.str());
  deuteron_look_back.ExportPDF(0,outfile);

  // Export equivalent proton PDFs
  std::stringstream protonfilename;
  protonfilename << base_path<< "/prt/replica_"<<replica<<".dat";
  ofstream protonfile; protonfile.open(protonfilename.str());
  ExportProton(pPDF, dPDFconfig, protonfile);
  protonfile.close();

  // Export best-fit parameters
  std::stringstream parfilename;
  parfilename << base_path<< "/par/parameters_"<<replica<<".dat";
  ofstream parfile; parfile.open(parfilename.str());
  deuteron_look_back.ExportPars(0,parfile);

  // Export theoretical predictions
  for (auto exp : experimental_data)
    for (int i=0; i<exp.GetNSet(); i++)
    {
      NNPDF::DataSet const& set = exp.GetSet(i);
      NNPDF::real* predictions = new NNPDF::real[set.GetNData()];
      ComputePredictions(&bpPDF, &deuteron_look_back, &set, predictions);

      std::stringstream datafilename;
      datafilename << base_path<< "/thp/"+set.GetSetName()+"_replica_"<<replica<<".dat";
      ofstream datafile; datafile.open(datafilename.str());

      for (int j=0; j<set.GetNData(); j++)
        datafile << predictions[j]<<std::endl;
      datafile.close();

      // Export filtered data
      if (replica == 0)
        set.Export(base_path + "/dat");
    }


  exit(0);
}

