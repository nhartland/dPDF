#pragma once
#include "NNPDF/pdfset.h"
#include "ns_network.h"
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
using NNPDF::real;
using std::vector;
  class DeuteronSet : public NNPDF::PDFSet
  {
  public:
    DeuteronSet(libconfig::Config const& s):
    PDFSet(s.lookup("fit.name"), s.lookup("fit.lambda"), ER_NONE),
    fParametrisation({2,5,3}),
    fBestFit(gsl_vector_calloc( fParametrisation.GetNParameters() ))
    {
      for (int i=0; i<fMembers; i++)
        fParameters.push_back(gsl_vector_calloc( fParametrisation.GetNParameters() ) );
      std::cout << "Initialised DeuteronSet with " << fParametrisation.GetNParameters() << " parameters " <<std::endl;
    };
    ~DeuteronSet(){ for (auto i : fParameters) gsl_vector_free(i); }

    void InitPDFSet() {return;};
    int const& GetNParameters() const {return fParametrisation.GetNParameters();}; 
    gsl_vector* GetParameters( int const& imem ) { return fParameters[imem]; };
    gsl_vector* GetBestFit()  { return fBestFit; };

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      real *xvals = new NNPDF::real[2]();
      xvals[0] = x; xvals[1] = static_cast<NNPDF::real>( log(x) );
      fParametrisation.Compute(fParameters[n], xvals, pdf);
      // for (int i=0; i<14; i++) std::cout << x <<"  "<<Q2<<"  "<<n<<"  "<<i <<"  "<<pdf[i] <<std::endl;
      delete[] xvals;
    	return;
    };

  private:
    NostateMLP fParametrisation;
    gsl_vector* fBestFit;
    vector<gsl_vector*> fParameters;
  };