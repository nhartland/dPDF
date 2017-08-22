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
  static const std::vector<int> activeFlavours = {1,2,3,10};
  static const int n_activeFlavours = static_cast<int>(activeFlavours.size());
  class DeuteronSet : public NNPDF::PDFSet
  {
  public:
    DeuteronSet(libconfig::Config const& s, int const& lambda):
    PDFSet(s.lookup("fit.name"), lambda, ER_NONE),
    fParametrisation({2,30, n_activeFlavours}),
    fGSLWork( gsl_integration_workspace_alloc (10000) ),
    nn_2(new double[fMembers*n_activeFlavours]),
    nn_norm(new double[fMembers*n_activeFlavours])
    {
      for (int i=0; i<fMembers; i++)
        fParameters.push_back(gsl_vector_calloc( fParametrisation.GetNParameters() ) );
      std::cout << "Initialised DeuteronSet with " << fParametrisation.GetNParameters() << " parameters " <<std::endl;
    };
    ~DeuteronSet(){ 
      for (auto i : fParameters) gsl_vector_free(i); 
      gsl_integration_workspace_free(fGSLWork);
      delete[] nn_2;
      delete[] nn_norm;
    }

    // Compute preprocessing
    void InitPDFSet() {
      std::fill(nn_2,    nn_2    + fMembers*n_activeFlavours, 0);   
      std::fill(nn_norm, nn_norm + fMembers*n_activeFlavours, 1);     
      for (int n=0; n<fMembers; n++)
      {
        // Compute large-x preprocessing
        std::array<NNPDF::real, 14> pdf;
        GetPDF(2.0,1,n, &pdf[0]); // Evaluate NN(2)
        for (int ifl=0; ifl<n_activeFlavours; ifl++)
          nn_2[n_activeFlavours*n + ifl] = pdf[activeFlavours[ifl]];

        GetPDF(0,1,n, &pdf[0]); // Evaluate NN(0)

        // Compute sum rules
        bool gslerror = false;
        const double pval = IntegratePDF(n,3,1,PDFSet::FX,gslerror,fGSLWork, 0.0, 2.0);
        if (gslerror) std::cout << "Valence alpha: " << fabs(0.5+0.1*gsl_vector_get(fParameters[n], fParametrisation.GetNParameters()-1)) <<std::endl;
        const double xsng = IntegratePDF(n,1,1,PDFSet::XFX,gslerror,fGSLWork, 0.0, 2.0);
        const double xglu = IntegratePDF(n,2,1,PDFSet::XFX,gslerror,fGSLWork, 0.0, 2.0);
        nn_norm[n_activeFlavours*n + 0] = (2.0-xsng)/xglu;
        nn_norm[n_activeFlavours*n + 2] = 6.0/pval;
      }
    };

    NostateMLP const& GetParametrisation() const {return fParametrisation;};
    int const& GetNParameters() const {return fParametrisation.GetNParameters();}; 
    gsl_vector* GetParameters( int const& imem ) { return fParameters[imem]; };

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      NNPDF::real* fitbasis = new NNPDF::real[n_activeFlavours];
      fParametrisation.Compute(fParameters[n], x, fitbasis);

      for (int i=0; i<14; i++) pdf[i] = 0;
      for (int i =0; i<n_activeFlavours; i++ )
        // pdf[activeFlavours[i]] = nn_norm[n_activeFlavours*n + i]*std::abs(std::abs(fitbasis[i]) - nn_2[n_activeFlavours*n + i]);
        pdf[activeFlavours[i]] = nn_norm[n_activeFlavours*n + i]*(fitbasis[i] - nn_2[n_activeFlavours*n + i]);
     
      // Valence preprocessing
      pdf[3] *= pow(x/2.0, fabs(0.5+0.1*gsl_vector_get(fParameters[n], fParametrisation.GetNParameters()-1))); // Valence low-x sum rule
      // pdf[10] = pdf[1]; // T8 = Singlet
      pdf[5] = pdf[3]; // V8 = Valence
      delete[] fitbasis; 
    	return;
    };

    void ExportPDF(int const& imem, std::ostream& os)
    {
      const double ymin = XGrid::appl_fy(1E-5); 
      const double ymax = XGrid::appl_fy(2.0);
      const int nx = 200;
 
      InitPDFSet(); 
      for (int i=0; i<nx; i++)
      {
        const NNPDF::real x = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx - 1))*i);
        std::array<NNPDF::real, 14> pdf; GetPDF(x,1,imem, &pdf[0]);
        os << x;
        for (int i =0; i<n_activeFlavours; i++ )
          os << "  "<< pdf[activeFlavours[i]];
        os <<std::endl;
      }
    }

    void ExportPars(int const& imem, std::ostream& os)
    {
      os << std::scientific << std::setprecision(20);
      for (int i=0; i<fParametrisation.GetNParameters(); i++)
        os << gsl_vector_get(fParameters[imem], i) <<std::endl;
    }

  private:
    NostateMLP fParametrisation;
    gsl_integration_workspace* fGSLWork; 
    vector<gsl_vector*> fParameters;
    double* nn_2;     // NN(2) for large-x preprocessing
    double* nn_norm;  // Multiplicative normalisation for PDF 
  };