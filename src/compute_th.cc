// Multi-beam test

#include "LHAPDF/LHAPDF.h"

// Configuration file
#include <libconfig.h++>
#include <gsl/gsl_integration.h>

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


void gsl_handler (const char * msg, const char * source, int line, int)
{ std::cerr << "gsl: " << source<<":"<<line<<": " <<msg <<std::endl;}

int main(int argc, char* argv[]) {
  LHAPDF::setVerbosity(0);
  NNPDF::SetVerbosity(0);
  gsl_set_error_handler (&gsl_handler);

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
  mkdir((base_path+"/dat").c_str(), 0777);
  mkdir((base_path+"/dat/systypes").c_str(), 0777);
  mkdir((base_path+"/thp").c_str(), 0777);
  mkdir((base_path+"/thd").c_str(), 0777);
  mkdir((base_path+"/pdf").c_str(), 0777);

  // Initialise and filter datasets
  std::vector<NNPDF::Experiment> experimental_data;
  ReadData(dPDFconfig, experimental_data);

  LHAPDFSet    proton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  IsoProtonSet isoproton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  DeuteronSet  deuteron = DeuteronSet::ReadSet(fitname, n_replicas);
  deuteron.Export(base_path+"/pdf");

  bool gslerror;
  gsl_integration_workspace *gsl = gsl_integration_workspace_alloc (10000);
  for (int i=0; i<n_replicas; i++)
  {
    const double xg = deuteron.IntegratePDF(i,LHAPDFSet::EVLN_GLU,1,PDFSet::XFX,gslerror,gsl, 0.0, 1.0);
    const double xs = deuteron.IntegratePDF(i,LHAPDFSet::EVLN_SNG,1,PDFSet::XFX,gslerror,gsl, 0.0, 1.0);
    std::cout << i 
              << "  xGLU: " << xg
              << "  xSNG: " << xs
              << "  xPDF: " << xg+xs <<std::endl;
  }
  
  std::array<NNPDF::PDFSet*,2> deuterons = {&deuteron, &isoproton};
  std::array<std::string,2>    labels    = {"thd", "thp"};

  for (auto exp : experimental_data)
    for (int i=0; i<exp.GetNSet(); i++)
    {
      NNPDF::DataSet const& set = exp.GetSet(i);
      set.Export(base_path + "/dat/");
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

