// n.p.hartland@ed.ac.uk  03/12

#include "LHAPDF/LHAPDF.h"

// Configuration file
#include <libconfig.h++>
#include <gsl/gsl_vector.h>

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

  const std::string fitname = dPDFconfig.lookup("fit.name");
  const int replica = atoi(argv[2]);
  libconfig::Setting& fitsetting = dPDFconfig.lookup("fit");
  libconfig::Setting & replica_setting = fitsetting.add ("replica", libconfig::Setting::TypeInt); replica_setting = replica;
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
  std::vector<NNPDF::Experiment> trainExp;
  std::vector<NNPDF::Experiment> validExp;
  InitData(dPDFconfig, trainExp, validExp);

  double nData = 0;
  for (size_t i=0; i<trainExp.size(); i++)
    nData += trainExp[i].GetNData();

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  cout              << "                            deuteron fitter                            "<<             endl;
  cout              << "                     Fit: "<<fitname<<"     Replica: "<<replica<<"   "<<             endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  // Initialise prototype parametrisation
  DeuteronSet dpdf(dPDFconfig);
  const int lambda = dpdf.GetMembers();
  const int nparam = dpdf.GetNParameters();
  CMAESMinimizer min(nparam, lambda, dPDFconfig.lookup("fit.sigma"));
  min.NormVect(dpdf.GetBestFit());

  for (int i=0; i< 1000; i++)
    min.Iterate(&dpdf, trainExp);

  NNPDF::real* chi2 = new NNPDF::real[lambda]();
  for (auto exp : trainExp)
    FastAddChi2(&dpdf, &exp, chi2);
  for (int i=0; i<lambda; i++)
    std::cout << i <<"  "<<chi2[i]/nData<<std::endl;
  delete[] chi2;

  std::stringstream filename;
  filename << "res/replica_"<<replica<<".dat";
  ofstream outfile; outfile.open(filename.str());
  dpdf.ExportBestFit(outfile);

  std::stringstream protonfilename;
  protonfilename << "prt/replica_"<<replica<<".dat";
  ofstream protonfile; protonfile.open(protonfilename.str());
  ExportProton(dPDFconfig, protonfile);


  exit(0);
}

