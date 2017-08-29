// Multi-beam test

#include "LHAPDF/LHAPDF.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"
#include "NNPDF/randomgenerator.h"

#include "deuteronset.h"


using namespace std;

/*
  This code tests the DeuteronSet constructor and copy-constructor
  If all goes well it should output a list of (1 0) pairs
*/

int main(int argc, char* argv[]) {
  const int nreplica = 100;
  const int nparam   = NostateMLP::get_nparam(pdf_architecture);

  const std::string base_path = "./res/base_parameters_40_hidden";
  std::vector<gsl_vector*> test_parameters;
  for (int i=0; i<nreplica; i++)
  {
    gsl_vector* parameters = gsl_vector_alloc(nparam);
    FILE * f = fopen ( (base_path + "/par/parameters_" + to_string(i)+ ".dat").c_str(), "r");
    assert( gsl_vector_fread (f, parameters) == 0 );
    test_parameters.push_back(parameters);
  }

  DeuteronSet    test_set(test_parameters, NNPDF::PDFSet::ER_MC);
  DeuteronSet    copy_set(test_set);

  NNPDF::FKTable test_table("./theory_65/FK_DYE886R_D.dat");
  NNPDF::ThPredictions test_pred(&test_set, &test_table);
  NNPDF::ThPredictions copy_pred(&copy_set, &test_table);
  NNPDF::ThPredictions ratio_pred = test_pred / copy_pred;
  ratio_pred.Print(std::cout);

  exit(0);
}

