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

double ComputeBestChi2(DeuteronSet& dpdf, LHAPDFSet const& pPDF, vector<Experiment> const& exps)
{
  dpdf.UseBestFit();
  double global_chi2 = 0;
  for (auto exp : exps)
  {
    NNPDF::real* chi2 = new NNPDF::real[dpdf.GetMembers()]();
    FastAddChi2(&pPDF, &dpdf, &exp, chi2);
    global_chi2 += chi2[0];
    std::cout << "Chi^2 for set " << exp.GetExpName() <<"  "<<chi2[0]/((double)exp.GetNData())<<std::endl;
    delete[] chi2;
  }
  return global_chi2;
}

int main(int argc, char* argv[]) {

  // Set LHAPDF verbosity
  LHAPDF::setVerbosity(0);

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

  double nData_trn = 0;
  for (size_t i=0; i<trainExp.size(); i++)
    nData_trn += trainExp[i].GetNData();

  double nData_val = 0;
  for (size_t i=0; i<validExp.size(); i++)
    nData_val += validExp[i].GetNData();

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  cout              << "                            deuteron fitter                            "<<             endl;
  cout              << "                     Fit: "<<fitname<<"     Replica: "<<replica<<"   "<<             endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  // Initialise prototype parametrisation
  DeuteronSet dpdf(dPDFconfig);
  const int lambda = dpdf.GetMembers();
  const int nparam = dpdf.GetNParameters();

  // Initialise proton parametrisation
  LHAPDFSet pPDF(dPDFconfig.lookup("fit.proton"), replica + 1, lambda);

  // Initialise minimiser
  CMAESMinimizer min(nparam, lambda, dPDFconfig.lookup("fit.sigma"));
  min.NormVect(dpdf.GetBestFit());
  const int ngen = dPDFconfig.lookup("fit.ngen");
  for (int i=0; i< ngen; i++)
  {
    std::cout << "Iteration: "<<i <<" / " <<ngen <<std::endl;
    min.Iterate(&pPDF, &dpdf, trainExp);

    if (i % 20 == 0 )
    {
      // Compute final chi2
      const double trnchi2 = ComputeBestChi2(dpdf, pPDF, trainExp)/nData_trn;
      std::cout << "Training chi2: " << trnchi2 <<std::endl;  

      const double valchi2 = ComputeBestChi2(dpdf, pPDF, validExp)/nData_val;
      std::cout << "Validation chi2: " << valchi2 <<std::endl;  
    } 
  }


  // Compute final chi2
  const double bfchi2 = ComputeBestChi2(dpdf, pPDF, trainExp)/nData_trn;
  std::cout << "Final chi2: " << bfchi2 <<std::endl;

  std::stringstream filename;
  filename << "res/replica_"<<replica<<".dat";
  ofstream outfile; outfile.open(filename.str());
  dpdf.ExportBestFit(outfile);

  std::stringstream protonfilename;
  protonfilename << "prt/replica_"<<replica<<".dat";
  ofstream protonfile; protonfile.open(protonfilename.str());
  ExportProton(pPDF, dPDFconfig, protonfile);
  protonfile.close();

  exit(0);
}

