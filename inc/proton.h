// proton.h
// Handling of proton 'c-factors'

#pragma once

#include <libconfig.h++>
#include "NNPDF/lhapdfset.h"
#include "deuteronset.h"

using NNPDF::LHAPDFSet;
using NNPDF::real;

void GenerateNormalisations(libconfig::Setting const& set, libconfig::Config const& settings);
void ExportProton(libconfig::Config const& settings, std::ostream& os);

class IsoProtonSet: public LHAPDFSet
{
public:
	IsoProtonSet(const std::string name, erType etype):
	LHAPDFSet(name, etype)
	{};

  IsoProtonSet(const std::string name, int member, erType etype):
  LHAPDFSet(name, member, etype)
  {};

 	virtual void GetPDF(real const& x, real const& Q2, int const& n, real* pdf) const
 	{
 		LHAPDFSet::GetPDF(x,Q2,n,pdf);
    for (int i=0; i<14; i++) pdf[i] *= 2.0;
 		pdf[EVLN_V3] = 0;
 		pdf[EVLN_T3] = 0;
 	};

};

void ExportProton(NNPDF::LHAPDFSet const& proton, libconfig::Config const& settings, std::ostream& os)
{
  const double Q0 = settings.lookup("qcd.Q0");
  const double ymin = XGrid::appl_fy(1E-5); 
  const double ymax = XGrid::appl_fy(2.0);
  const int nx = 200;

  for (int i=0; i<nx; i++)
  {
    const NNPDF::real x_d = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx ))*i);
    const NNPDF::real x_p = (x_d > 1) ? 1:x_d;
    std::array<NNPDF::real, 14> pdf; 
    proton.GetPDF(x_p, Q0*Q0 ,0, &pdf[0]);
    os << x_d;
    for (int i =0; i<n_activeFlavours; i++ )
      os << "  "<< 2.0*pdf[activeFlavours[i]];
    os <<std::endl;
  }
}