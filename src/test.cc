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


  NostateMLP network({2,3,5});
  gsl_vector* parameters = gsl_vector_calloc( network.GetNParameters() );
  CMAESMinimizer min(network.GetNParameters(), 1, 0.1);
  min.NormVect(parameters);
  for (int i=0; i<network.GetNParameters(); i++) std::cout << gsl_vector_get(parameters, i) <<std::endl;


  NNPDF::real* ovals = new NNPDF::real[5]();
  network.Compute(parameters, 0.1, ovals);
  for (int i=0; i<5; i++) std::cout << i <<"  "<<ovals[i]<<std::endl;
  // std::cout << std::endl;
  // min.NormVect(parameters);
  // network.Compute(parameters, xvals, ovals);
  // for (int i=0; i<5; i++) std::cout << i <<"  "<<ovals[i]<<std::endl;

  exit(0);
}

