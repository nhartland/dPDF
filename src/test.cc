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

#include "ns_network.h"
#include "fastaddchi2.h"
#include "deuteronset.h"
#include "cmaes.h"
#include "filter.h"
#include "colour.h"

using namespace Colour;
using namespace std;

int main(int argc, char* argv[]) {
  NNPDF::RandomGenerator::InitRNG(0,0);

  NNPDF::LHAPDFSet protonset("NNPDF30_nlo_as_0118", NNPDF::PDFSet::ER_MC);

  exit(0);
}

