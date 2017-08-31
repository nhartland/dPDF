// Multi-beam test

#include "LHAPDF/LHAPDF.h"

// Configuration file
#include <libconfig.h++>

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

#include "fastaddchi2.h"
#include "deuteronset.h"
#include "cmaes.h"
#include "filter.h"
#include "colour.h"
#include "proton.h"

using namespace Colour;
using namespace std;

int main(int argc, char* argv[]) {
  LHAPDF::setVerbosity(0);
  NNPDF::SetVerbosity(0);

  if (argc < 2)
  {
    cerr << "Usage: compute_th <configuration file>"<<std::endl;
    exit(-1);
  }
  // Read configuration file
  libconfig::Config dPDFconfig;
  dPDFconfig.readFile(argv[1]);

  const std::string fitname = dPDFconfig.lookup("fit.name");
  const std::string base_path = "./res/"+fitname;
  const int n_replicas = dPDFconfig.lookup("fit.nrep");

  // Make directories
  mkdir((base_path+"/thr").c_str(), 0777);

  // Initialise and filter datasets
  std::vector<NNPDF::FKSet> plot_data;
  ReadPlots(dPDFconfig, plot_data);

  LHAPDFSet    proton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  IsoProtonSet isoproton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  DeuteronSet deuteron = DeuteronSet::ReadSet(fitname, n_replicas);

  for (auto set : plot_data)
  {
      NNPDF::real* predictions = new NNPDF::real[n_replicas*set.GetNDataFK()];
      ComputePredictions(&proton, &deuteron, &set, predictions);
      NNPDF::ThPredictions d_pred(&deuteron, &set, predictions);

      ComputePredictions(&proton, &isoproton, &set, predictions);
      NNPDF::ThPredictions p_pred(&isoproton, &set, predictions);
      delete[] predictions;

      NNPDF::ThPredictions c_pred = d_pred / p_pred;

      ofstream file;
      std::cout << "Writing to " << base_path + "/thr/TH_"+set.GetDataName()+".dat" <<std::endl;
      file.open (base_path + "/thr/TH_"+set.GetDataName()+".dat" );
      c_pred.Print(file); file.close();
    }

  exit(0);
}