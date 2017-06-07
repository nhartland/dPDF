#pragma once
#include "NNPDF/pdfset.h"
#include "ns_network.h"
#include "fk_xgrid.h"
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
#include <array>
using NNPDF::real;
using std::vector;

  static const std::vector<int> activeFlavours = {1,2}; // 10
  static const int n_activeFlavours = static_cast<int>(activeFlavours.size());

  class DeuteronSet : public NNPDF::PDFSet
  {
  public:
    DeuteronSet(libconfig::Config const& s):
    PDFSet(s.lookup("fit.name"), s.lookup("fit.lambda"), ER_NONE),
    fParametrisation({2,20, n_activeFlavours}),
    fBestFit(gsl_vector_calloc( fParametrisation.GetNParameters() )),
    nn_1(new double[fMembers*n_activeFlavours])
    {
      for (int i=0; i<fMembers; i++)
        fParameters.push_back(gsl_vector_calloc( fParametrisation.GetNParameters() ) );
      std::cout << "Initialised DeuteronSet with " << fParametrisation.GetNParameters() << " parameters " <<std::endl;
    };
    ~DeuteronSet(){ 
      for (auto i : fParameters) gsl_vector_free(i); 
      delete[] nn_1;
    }

    // Compute preprocessing
    void InitPDFSet() {
      std::fill(nn_1, nn_1 + fMembers*n_activeFlavours, 0);    
      for (int n=0; n<fMembers; n++)
      {
        std::array<NNPDF::real, 14> pdf;
        GetPDF(1,1,n, &pdf[0]);
        for (int ifl=0; ifl<n_activeFlavours; ifl++)
          nn_1[n_activeFlavours*n + ifl] = pdf[activeFlavours[ifl]];
      }
    };

    NostateMLP const& GetParametrisation() const {return fParametrisation;};
    int const& GetNParameters() const {return fParametrisation.GetNParameters();}; 
    gsl_vector* GetParameters( int const& imem ) { return fParameters[imem]; };
    gsl_vector* GetBestFit()  { return fBestFit; };

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      NNPDF::real* fitbasis = new NNPDF::real[n_activeFlavours];
      fParametrisation.Compute(fParameters[n], x, fitbasis);
      for (int i=0; i<14; i++) pdf[i] = 0;
      for (int i =0; i<n_activeFlavours; i++ )
        pdf[activeFlavours[i]] = std::abs(std::abs(fitbasis[i]) - nn_1[n_activeFlavours*n + i]);
      pdf[10] = pdf[1]; // T8 = Singlet
      delete[] fitbasis; 
    	return;
    };

    void ExportBestFit(std::ostream& os)
    {
      const double ymin = XGrid::appl_fy(1E-5); 
      const double ymax = XGrid::appl_fy(1.5);
      const int nx = 200;
 
      // Copy best-fit to fParameters[0]
      gsl_vector_memcpy (fParameters[0], fBestFit);
      InitPDFSet(); 
      for (int i=0; i<nx; i++)
      {
        const NNPDF::real x = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx - 1))*i);
        std::array<NNPDF::real, 14> pdf; GetPDF(x,1,0, &pdf[0]);
        os << x;
        for (int i =0; i<n_activeFlavours; i++ )
          os << "  "<< pdf[activeFlavours[i]];
        os <<std::endl;
      }
    }

  private:
    NostateMLP fParametrisation;
    gsl_vector* fBestFit;
    vector<gsl_vector*> fParameters;
    double* nn_1;
  };