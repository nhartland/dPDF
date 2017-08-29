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
  const int nparam = NostateMLP::get_nparam(pdf_architecture);

  // Make directories
  mkdir((base_path+"/thp").c_str(), 0777);
  mkdir((base_path+"/thd").c_str(), 0777);

  // Initialise and filter datasets
  std::vector<NNPDF::Experiment> experimental_data;
  ReadData(dPDFconfig, experimental_data);

  LHAPDFSet   proton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  IsoProtonSet isoproton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  DeuteronSet deuteron = DeuteronSet::ReadSet(fitname, n_replicas);

  std::array<NNPDF::PDFSet*,2> deuterons = {&deuteron, &isoproton};
  std::array<std::string,2>    labels    = {"thd", "thp"};

  for (auto exp : experimental_data)
    for (int i=0; i<exp.GetNSet(); i++)
    {
      NNPDF::DataSet const& set = exp.GetSet(i);
      for (int j=0; j<deuterons.size(); j++)
      {
        NNPDF::real* predictions = new NNPDF::real[n_replicas*set.GetNData()];
        ComputePredictions(&proton, deuterons[j], &set, predictions);
        NNPDF::ThPredictions pred(deuterons[j], &set, predictions);
        delete[] predictions;

        ofstream file;
        std::cout << "Writing to " << base_path + "/"+labels[j]+"/TH_"+set.GetSetName()+".dat" <<std::endl;
        file.open (base_path + "/"+labels[j]+"/TH_"+set.GetSetName()+".dat" );
        pred.Print(file); file.close();
      }
    }

  exit(0);
}

