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
    cerr << "Usage: plotter <configuration file>"<<std::endl;
    exit(-1);
  }
  // Read configuration file
  libconfig::Config dPDFconfig;
  dPDFconfig.readFile(argv[1]);

  const std::string fitname = dPDFconfig.lookup("fit.name");
  const std::string base_path = "./res/"+fitname;
  const int n_replicas = dPDFconfig.lookup("fit.nrep");
  const int nparam = NostateMLP::get_nparam(pdf_architecture);

  // Initialise and filter datasets
  std::vector<NNPDF::Experiment> experimental_data;
  ReadData(dPDFconfig, experimental_data);

  double nData = 0;
  for (auto exp: experimental_data)
    nData += exp.GetNData();

  // Read fit parameters
  std::vector<gsl_vector*> fit_parameters;
  for (int i=0; i<n_replicas; i++)
  {
    gsl_vector* parameters = gsl_vector_alloc(nparam);
    FILE * f = fopen ( (base_path + "/par/parameters_" + to_string(i)+ ".dat").c_str(), "r");
    assert( gsl_vector_fread (f, parameters) == 0 );
    fit_parameters.push_back(parameters);
  }

  LHAPDFSet   proton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  DeuteronSet deuteron(fit_parameters, NNPDF::PDFSet::ER_MC);

  exit(0);
}

