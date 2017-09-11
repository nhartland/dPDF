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

using std::vector;
using std::distance;
using std::find;
using namespace NNPDF;

static const vector<int> activeFlavours = {PDFSet::EVLN_SNG, PDFSet::EVLN_GLU};//, PDFSet::EVLN_T8};
static const int n_activeFlavours = static_cast<int>(activeFlavours.size());
static const vector<int> pdf_architecture = {2,10, n_activeFlavours};

static const ptrdiff_t FIT_SNG = distance(activeFlavours.begin(), find(activeFlavours.begin(), activeFlavours.end(), PDFSet::EVLN_SNG));
static const ptrdiff_t FIT_GLU = distance(activeFlavours.begin(), find(activeFlavours.begin(), activeFlavours.end(), PDFSet::EVLN_GLU));
static const ptrdiff_t FIT_VAL = distance(activeFlavours.begin(), find(activeFlavours.begin(), activeFlavours.end(), PDFSet::EVLN_VAL));
static const ptrdiff_t FIT_T8  = distance(activeFlavours.begin(), find(activeFlavours.begin(), activeFlavours.end(), PDFSet::EVLN_T8));
static const ptrdiff_t FIT_V8  = distance(activeFlavours.begin(), find(activeFlavours.begin(), activeFlavours.end(), PDFSet::EVLN_V8));
static bool fitting( ptrdiff_t pdf ) {return pdf < n_activeFlavours;};

class DeuteronSet : public PDFSet
{
public:
  DeuteronSet(std::vector<gsl_vector*> const& parameters, erType error = ER_NONE):
  PDFSet("Deuteron", parameters.size(), error),
  fParametrisation(pdf_architecture),
  fGSLWork( gsl_integration_workspace_alloc (10000) ),
  nn_1(fMembers*n_activeFlavours,0),
  nn_norm(fMembers*n_activeFlavours,1)
  {
    const double Q0dum = 0;
    for (size_t n=0; n<fMembers; n++)
    {
      // Copy parameters
      gsl_vector* newpar = gsl_vector_calloc(fParametrisation.GetNParameters());
      gsl_vector_memcpy(newpar, parameters[n]); fParameters.push_back(newpar);

      // Compute large-x preprocessing
      std::array<real, 14> pdf;
      GetPDF(1.0, Q0dum, n, &pdf[0]); // Evaluate NN(1)
      for (int ifl=0; ifl<n_activeFlavours; ifl++)
        nn_1[n_activeFlavours*n + ifl] = pdf[activeFlavours[ifl]];

      // Compute sum rules
      bool gslerror = false;
      if ( fitting(FIT_SNG) && fitting(FIT_GLU) ) // A bit superfluous this
      {  
        const double xsng = IntegratePDF(n,EVLN_SNG,Q0dum,PDFSet::XFX,gslerror,fGSLWork, 0.0, 1.0);
        const double xglu = IntegratePDF(n,EVLN_GLU,Q0dum,PDFSet::XFX,gslerror,fGSLWork, 0.0, 1.0);
        nn_norm[n_activeFlavours*n + FIT_SNG] = (1.0-xglu)/xsng;
      }
      if ( fitting(FIT_VAL) )
      {
        const double pval = IntegratePDF(n,EVLN_VAL,Q0dum,PDFSet::FX,gslerror,fGSLWork, 0.0, 1.0);
        nn_norm[n_activeFlavours*n + FIT_VAL] = 6.0/pval;
      }
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

  virtual void GetPDF(real const& x, real const& Q2, int const& n, real* pdf) const
  {
    real* fitbasis = new real[n_activeFlavours]();
    fParametrisation.Compute(fParameters[n], x, fitbasis);

    for (int i=0; i<14; i++) pdf[i] = 0;
    for (int i =0; i<n_activeFlavours; i++ )
      pdf[activeFlavours[i]] = nn_norm[n_activeFlavours*n + i]*(fitbasis[i] - nn_1[n_activeFlavours*n + i]);
    delete[] fitbasis; 

    // Valence preprocessing
    if ( fitting(FIT_VAL) ) pdf[EVLN_VAL] *= pow(x, fabs(0.5+0.1*gsl_vector_get(fParameters[n], fParametrisation.GetNParameters()-1)));
    if ( !fitting(FIT_T8) ) pdf[EVLN_T8] = pdf[EVLN_SNG];   // T8 = Singlet
    if ( !fitting(FIT_V8) ) pdf[EVLN_V8] = pdf[EVLN_VAL];   // V8 = Valence

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
        const real x = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx - 1))*i);
        std::array<real, 14> pdf; GetPDF(x,1,imem, &pdf[0]);
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
