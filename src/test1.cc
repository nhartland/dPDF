// Multi-beam test

#include "LHAPDF/LHAPDF.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"
#include "NNPDF/randomgenerator.h"

#include "proton.h"

using namespace std;

/*
  This code tests the proton-deuteron rotation in APFEL and in
  IsoProtonSet. We compute the predictions for \sigma_pd by
    a) Using the \sigma_pd FK table from APFEL
    b) Using a \sigma_pp FK table and setting T3=V3=0
  Under the assumptions used in NNPDF these should be the same.
  NOTE: The proton-deuteron rotation in APFEL assumes data is given /per nucleon/
        hence the final factor of /2.0
*/

int main(int argc, char* argv[]) {
  NNPDF::RandomGenerator::InitRNG(0,0);

  NNPDF::LHAPDFSet proton  ("NNPDF30_nlo_as_0118", NNPDF::PDFSet::ER_MC);
  IsoProtonSet     deuteron("NNPDF30_nlo_as_0118", NNPDF::PDFSet::ER_MC);

  NNPDF::FKTable   p_table("./theory_65/FK_DYE886R_P.dat");
  NNPDF::FKTable   d_table("./theory_65/FK_DYE886R_D.dat");

  NNPDF::ThPredictions singlebeam(&proton, &d_table);
  NNPDF::ThPredictions multibeam(&proton, &deuteron, &p_table);

  NNPDF::ThPredictions ratio = multibeam/singlebeam;
  for (int i=0; i<ratio.GetNData(); i++)
    std::cout << ratio.GetObsCV(i)/2.0 <<std::endl;

  exit(0);
}

