// Multi-beam test

#include "LHAPDF/LHAPDF.h"
#include "APFEL/APFEL.h"

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

// External PDF set
PDFSet* APFELSet = 0;
extern "C" void externalsetapfel_(double const& x, double const& Q, double *pdf)
{
  NNPDF::real* EVLN = new NNPDF::real[14];
  NNPDF::real* LHA = new NNPDF::real[14];
  APFELSet->GetPDF(x, Q*Q, 0, EVLN);
  NNPDF::PDFSet::EVLN2LHA(EVLN, LHA);

  for (int i=0; i<13; i++)
    pdf[i] = LHA[i];

  delete[] EVLN;
  delete[] LHA;
}

int main(int argc, char* argv[]) {
  LHAPDF::setVerbosity(0);
  NNPDF::SetVerbosity(0);

  if (argc < 2)
  {
    cerr << "Usage: compute_f2 <configuration file>"<<std::endl;
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
  mkdir((base_path+"/f2r").c_str(), 0777);

  LHAPDFSet   proton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  IsoProtonSet isoproton(dPDFconfig.lookup("fit.proton"), NNPDF::PDFSet::ER_MC);
  DeuteronSet deuteron = DeuteronSet::ReadSet(fitname, n_replicas);

  // Initial and final scales
  const double Q0 = 1;
  const double Q  = 100;

  // Initialize APFEL
  APFELSet = &deuteron;
  APFEL::SetPDFSet("external");  // If you comment out this line you should get the same result.
  APFEL::InitializeAPFEL_DIS();
  APFEL::EvolveAPFEL(Q0,Q);

  exit(0);
}

