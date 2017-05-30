#pragma once
#include "NNPDF/pdfset.h"
#include "ns_network.h"
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
using NNPDF::real;
using std::vector;

  static const std::vector<int> activeFlavours = {1,2,10,11};
  static const int n_activeFlavours = static_cast<int>(activeFlavours.size());

  class DeuteronSet : public NNPDF::PDFSet
  {
  public:
    DeuteronSet(libconfig::Config const& s):
    PDFSet(s.lookup("fit.name"), s.lookup("fit.lambda"), ER_NONE),
    fParametrisation({2,10, n_activeFlavours}),
    fBestFit(gsl_vector_calloc( fParametrisation.GetNParameters() ))
    {
      for (int i=0; i<fMembers; i++)
        fParameters.push_back(gsl_vector_calloc( fParametrisation.GetNParameters() ) );
      std::cout << "Initialised DeuteronSet with " << fParametrisation.GetNParameters() << " parameters " <<std::endl;
    };
    ~DeuteronSet(){ for (auto i : fParameters) gsl_vector_free(i); }

    void InitPDFSet() {return;};
    NostateMLP const& GetParametrisation() const {return fParametrisation;};
    int const& GetNParameters() const {return fParametrisation.GetNParameters();}; 
    gsl_vector* GetParameters( int const& imem ) { return fParameters[imem]; };
    gsl_vector* GetBestFit()  { return fBestFit; };

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      NNPDF::real* fitbasis = new NNPDF::real[n_activeFlavours];
      fParametrisation.Compute(fParameters[n], x, fitbasis);
      for (int i =0; i<n_activeFlavours; i++ )
        pdf[activeFlavours[i]] = fitbasis[i];
      delete[] fitbasis; 
    	return;
    };

  private:
    NostateMLP fParametrisation;
    gsl_vector* fBestFit;
    vector<gsl_vector*> fParameters;
  };