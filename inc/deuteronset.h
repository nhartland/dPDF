#pragma once
#include "NNPDF/pdfset.h"
#include "ns_network.h"
#include "fk_xgrid.h"
#include <libconfig.h++>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <limits>
#include <array>
#include <iomanip>

using NNPDF::real;
using std::vector;

  // static const std::vector<int> activeFlavours = {1,2,3,5,10};
  // static const int n_activeFlavours = static_cast<int>(activeFlavours.size());
  // static const std::vector<int> pdf_architecture = {2,40, n_activeFlavours};

  static const std::vector<int> activeFlavours = {1,2};
  static const int n_activeFlavours = static_cast<int>(activeFlavours.size());
  static const std::vector<int> pdf_architecture = {2,10, n_activeFlavours};


  class DeuteronSet : public NNPDF::PDFSet
  {
  public:
    DeuteronSet(std::vector<gsl_vector*> const& parameters, erType error = ER_NONE):
    PDFSet("Deuteron", parameters.size(), error),
    fParametrisation(pdf_architecture),
    fGSLWork( gsl_integration_workspace_alloc (10000) ),
    nn_1(fMembers*n_activeFlavours,0),
    nn_norm(fMembers*n_activeFlavours,1)
    {
      for (size_t n=0; n<fMembers; n++)
      {
        // Copy parameters
        gsl_vector* newpar = gsl_vector_calloc(fParametrisation.GetNParameters());
        gsl_vector_memcpy(newpar, parameters[n]); fParameters.push_back(newpar);

        // Compute large-x preprocessing
        std::array<NNPDF::real, 14> pdf;
        GetPDF(1.0,1,n, &pdf[0]); // Evaluate NN(1)
        for (int ifl=0; ifl<n_activeFlavours; ifl++)
          nn_1[n_activeFlavours*n + ifl] = pdf[activeFlavours[ifl]];
        GetPDF(0,1,n, &pdf[0]); // Evaluate NN(0)

        // Compute sum rules
        bool gslerror = false;
        // const double pval = IntegratePDF(n,3,1,PDFSet::FX,gslerror,fGSLWork, 0.0, 1.0);
        const double xsng = IntegratePDF(n,1,1,PDFSet::XFX,gslerror,fGSLWork, 0.0, 1.0);
        const double xglu = IntegratePDF(n,2,1,PDFSet::XFX,gslerror,fGSLWork, 0.0, 1.0);
        nn_norm[n_activeFlavours*n + 0] = (1.0-xsng)/xglu;
        // nn_norm[n_activeFlavours*n + 2] = 6.0/pval;
        if(gslerror) nn_norm[n_activeFlavours*n] = std::numeric_limits<double>::infinity(); // Integration error: set gluon norm to infty
      }
    };

    // Copy-constructor
    DeuteronSet(DeuteronSet const& other):
    PDFSet(other),
    fParametrisation(other.fParametrisation),
    fGSLWork( gsl_integration_workspace_alloc (10000) ),
    nn_1(other.nn_1),
    nn_norm(other.nn_norm)
    {
      for (size_t n=0; n<fMembers; n++)
      {
        gsl_vector* newpar = gsl_vector_calloc(fParametrisation.GetNParameters());
        gsl_vector_memcpy(newpar, other.fParameters[n]); fParameters.push_back(newpar);
      }
    };


    ~DeuteronSet(){ 
      for (auto i : fParameters) gsl_vector_free(i); 
      gsl_integration_workspace_free(fGSLWork);
    }

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      NNPDF::real* fitbasis = new NNPDF::real[n_activeFlavours];
      fParametrisation.Compute(fParameters[n], x, fitbasis);

      for (int i=0; i<14; i++) pdf[i] = 0;
      for (int i =0; i<n_activeFlavours; i++ )
        pdf[activeFlavours[i]] = nn_norm[n_activeFlavours*n + i]*(fitbasis[i] - nn_1[n_activeFlavours*n + i]);
     
      // Valence preprocessing
      // pdf[3] *= pow(x, fabs(0.5+0.1*gsl_vector_get(fParameters[n], fParametrisation.GetNParameters()-1))); // Valence low-x sum rule
      pdf[10] = pdf[1]; // T8 = Singlet
      // pdf[5] = pdf[3]; // V8 = Valence
      delete[] fitbasis; 
    	return;
    };

    void Export(std::string const& prefix)
    {
      const double ymin = XGrid::appl_fy(1E-5); 
      const double ymax = XGrid::appl_fy(1.0);
      const int nx = 200;
 
      for (int imem = 0; imem<fMembers; imem++)
      {
        std::ofstream os; os.open(prefix + "/replica_"+std::to_string(imem)+".dat");
        for (int i=0; i<nx; i++)
        {
          const NNPDF::real x = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx - 1))*i);
          std::array<NNPDF::real, 14> pdf; GetPDF(x,1,imem, &pdf[0]);
          os << x;
          for (int i =0; i<14; i++ ) os << "  "<< pdf[i];
          os <<std::endl;
        }
        os.close();
      }
    }

    // Read a set of deuteron fit parameters from file
    static DeuteronSet ReadSet(std::string const& setname, int const& n_replicas)
    {
        const int nparam = NostateMLP::get_nparam(pdf_architecture);
        const std::string base_path = "./res/"+setname;
        std::vector<gsl_vector*> fit_parameters;
        for (int i=0; i<n_replicas; i++)
        {
          gsl_vector* parameters = gsl_vector_alloc(nparam);
          FILE * f = fopen ( (base_path + "/par/parameters_" + std::to_string(i)+ ".dat").c_str(), "r");
          assert( gsl_vector_fread (f, parameters) == 0 );
          fit_parameters.push_back(parameters);
        }

        DeuteronSet set(fit_parameters, NNPDF::PDFSet::ER_MC);
        for (auto ipar : fit_parameters) gsl_vector_free(ipar);
        return set;
    }

  private:
    NostateMLP fParametrisation;
    gsl_integration_workspace* fGSLWork; 
    vector<gsl_vector*> fParameters;
    vector<double> nn_1;     // NN(2) for large-x preprocessing
    vector<double> nn_norm;  // Multiplicative normalisation for PDF 
  };