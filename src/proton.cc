// proton.h
// Handling of proton 'c-factors'
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>

#include <libconfig.h++>

#include "deuteronset.h"

void GenerateNormalisations(libconfig::Setting const& set, libconfig::Config const& settings)
{
  if (!set.exists("normby_proton")) return;
  const int replica = settings.lookup("fit.replica");
  const NNPDF::LHAPDFSet proton(settings.lookup("fit.proton"),  NNPDF::LHAPDFSet::ER_MC);

  const int theoryIndex = settings.lookup("qcd.THID");
  const char *tablename = set["normby_proton"];
  const char *setname   = set["name"];
  std::stringstream fkpath;
  fkpath <<  "theory_"<<theoryIndex<<"/FK_"<<tablename<<".dat";

  const NNPDF::FKTable tab(fkpath.str());
  const NNPDF::ThPredictions theory(&proton, &tab);

  std::stringstream opath;
  opath << "res/norm/CF_"<< setname << "_"<<replica<<".dat";
  std::ofstream os; os.open(opath.str());

  os  <<  "*******************************************************************************************" << std::endl
    <<  " SetName: "<< setname << std::endl
    <<  " Author:  dPDF" << std::endl
    <<  " Date:    Today" << std::endl
    <<  " CodesUsed: None " << std::endl
    <<  " TheoryInput: None" << std::endl
    <<  " PDFset: None" << std::endl
    <<  " Warnings: None" << std::endl
    <<  " ********************************************************************************************" << std::endl;
  for (int i=0; i< theory.GetNData(); i++)
    os << 1.0/theory.GetObs(i, replica) <<"  "<<0<<std::endl;
  os.close();
}

void ExportProton(libconfig::Config const& settings, std::ostream& os)
{
  const int replica = settings.lookup("fit.replica");
  const NNPDF::LHAPDFSet proton(settings.lookup("fit.proton"),  NNPDF::LHAPDFSet::ER_MC);

  const double ymin = XGrid::appl_fy(1E-5); 
  const double ymax = XGrid::appl_fy(1.5);
  const int nx = 200;

  for (int i=0; i<nx; i++)
  {
    const NNPDF::real x = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx - 1))*i);
    std::array<NNPDF::real, 14> pdf; proton.GetPDF(x,1.0,replica, &pdf[0]);
    os << x;
    for (int i =0; i<n_activeFlavours; i++ )
      os << "  "<< pdf[activeFlavours[i]];
    os <<std::endl;
  }
}

