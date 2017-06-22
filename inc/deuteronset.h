#pragma once
#include "NNPDF/pdfset.h"
#include "ns_network.h"
#include "fk_xgrid.h"
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <limits>
#include <array>
using NNPDF::real;
using std::vector;

  // static const std::vector<int> activeFlavours = {1,2,3,5,10};
  // static const std::vector<int> activeFlavours = {1,2,3};
  static const std::vector<int> activeFlavours = {1,2};
  static const int n_activeFlavours = static_cast<int>(activeFlavours.size());
  class DeuteronSet : public NNPDF::PDFSet
  {
  public:
    DeuteronSet(libconfig::Config const& s):
    PDFSet(s.lookup("fit.name"), s.lookup("fit.lambda"), ER_NONE),
    fParametrisation({2,30, n_activeFlavours}),
    fGSLWork( gsl_integration_workspace_alloc (10000) ),
    fBestFit(gsl_vector_calloc( fParametrisation.GetNParameters() )),
    nn_1(new double[fMembers*n_activeFlavours]),
    nn_norm(new double[fMembers*n_activeFlavours])
    {
      for (int i=0; i<fMembers; i++)
        fParameters.push_back(gsl_vector_calloc( fParametrisation.GetNParameters() ) );
      std::cout << "Initialised DeuteronSet with " << fParametrisation.GetNParameters() << " parameters " <<std::endl;
    };
    ~DeuteronSet(){ 
      for (auto i : fParameters) gsl_vector_free(i); 
      gsl_integration_workspace_free(fGSLWork);
      delete[] nn_1;
      delete[] nn_norm;
    }

    // Compute preprocessing
    void InitPDFSet() {
      std::fill(nn_1,    nn_1    + fMembers*n_activeFlavours, 0);   
      std::fill(nn_norm, nn_norm + fMembers*n_activeFlavours, 1);     
      for (int n=0; n<fMembers; n++)
      {
        // Compute large-x preprocessing
        std::array<NNPDF::real, 14> pdf;
        GetPDF(2.0,1,n, &pdf[0]); // Evaluate NN(x_max)
        for (int ifl=0; ifl<n_activeFlavours; ifl++)
          nn_1[n_activeFlavours*n + ifl] = pdf[activeFlavours[ifl]];

        // Compute MSR normalisation
        bool gslerror = false;
        const double xsng = IntegratePDF(n,1,1,PDFSet::XFX,gslerror,fGSLWork, 0, 2);
        const double xglu = IntegratePDF(n,2,1,PDFSet::XFX,gslerror,fGSLWork, 0, 2);
        nn_norm[n_activeFlavours*n +1] = (1.0-xsng)/xglu;
      }
    };

    NostateMLP const& GetParametrisation() const {return fParametrisation;};
    int const& GetNParameters() const {return fParametrisation.GetNParameters();}; 
    gsl_vector* GetParameters( int const& imem ) { return fParameters[imem]; };
    gsl_vector* GetBestFit()  { return fBestFit; };

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      NNPDF::real* fitbasis = new NNPDF::real[n_activeFlavours];
      fParametrisation.Compute(fParameters[n], x/2.0, fitbasis);
      for (int i=0; i<14; i++) pdf[i] = 0;
      for (int i =0; i<n_activeFlavours; i++ )
        // pdf[activeFlavours[i]] = nn_norm[n_activeFlavours*n + i]*std::abs(std::abs(fitbasis[i]) - nn_1[n_activeFlavours*n + i]);
        pdf[activeFlavours[i]] = nn_norm[n_activeFlavours*n + i]*(fitbasis[i] - nn_1[n_activeFlavours*n + i]);

      pdf[10] = pdf[1]; // T8 = Singlet
      pdf[5] = pdf[3]; // V8 = Valence
      delete[] fitbasis; 
    	return;
    };

    void UseBestFit() {
      for (int i=0; i<fMembers; i++)
        gsl_vector_memcpy (fParameters[i], fBestFit);
      InitPDFSet();
    }

    void ExportBestFit(std::ostream& os)
    {
      const double ymin = XGrid::appl_fy(1E-5); 
      const double ymax = XGrid::appl_fy(2.0);
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
    gsl_integration_workspace* fGSLWork; 
    gsl_vector* fBestFit;
    vector<gsl_vector*> fParameters;
    double* nn_1;     // NN(1) for large-x preprocessing
    double* nn_norm;  // Multiplicative normalisation for PDF 
  };