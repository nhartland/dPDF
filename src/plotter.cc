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

  // Initialise and filter datasets
  std::vector<NNPDF::Experiment> experimental_data;
  ReadData(dPDFconfig, experimental_data);

  double nData = 0;
  for (auto exp: experimental_data)
    nData += exp.GetNData();

  exit(0);
}

