// FK table generation

#include "LHAPDF/LHAPDF.h"
#include "APFEL/APFEL.h"
#include "APFEL/APFELdev.h"

// Configuration file
#include <libconfig.h++>

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"
#include "NNPDF/nnpdfdb.h"
#include "NNPDF/fkgenerator.h"

#include "fk_xgrid.h"
#include "fastaddchi2.h"
#include "deuteronset.h"
#include "cmaes.h"
#include "filter.h"
#include "colour.h"
#include "proton.h"

using namespace Colour;
using namespace std;

void gsl_handler (const char * msg, const char * source, int line, int)
{ std::cerr << "gsl: " << source<<":"<<line<<": " <<msg <<std::endl;}

// Initialise APFEL for evolution factors
// Note here nx is the number of x-points to be output to the Fk table
// Therefore it doesn't include x=1. This is added manually in this function
void initEvolgrid(int const& nx, double const& xmin)
{
  double* xg = new double[nx+1];
  std::cout << " Initialising  "<< nx << " points starting from x = " << xmin <<std::endl; 

  const double ymin = XGrid::appl_fy(0.99*xmin); 
  const double ymax = XGrid::appl_fy(1.0);
  
  // Populate grid
  for (int i=0; i<=nx; i++)
    xg[i] = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx))*i);
  xg[nx] = 1;
  
  // Set evolution operator parameters
  APFEL::SetFastEvolution(false);
  APFEL::LockGrids(false);
  APFEL::SetEpsilonTruncation(1E-1);
  APFEL::EnableEvolutionOperator(true); 
  APFEL::SetNumberOfGrids(1);           
  APFEL::SetExternalGrid(1,nx,3,xg);  

  delete[] xg;
  return;
}


double QC = -1;
double diskernel(std::string const& obs, double const& x, double const& Q, double const& y, int const& i, int const& beta)
{
  if (Q < APFEL::GetMuF0()) return 0;
  if (Q != QC)
  {
    APFEL::SetFKObservable(obs);
    APFEL::ComputeStructureFunctionsAPFEL(APFEL::GetMuF0(),Q);
    QC = Q;
  }
  return APFEL::FKSimulator(x,Q,y,i,beta);
}

double getQ2Max_CommonData(NNPDF::CommonData const& cd)
{
  double q2max = 0;
  for (int i=0; i<cd.GetNData(); i++)
   q2max = fmax(q2max, cd.GetKinematics(i,1));
  return q2max;
}

double getXmin_CommonData(NNPDF::CommonData const& cd)
{
  double xmin = 1;
  for (int i=0; i<cd.GetNData(); i++)
   xmin = fmin(xmin, cd.GetKinematics(i,0));
  return xmin;
}

int main(int argc, char* argv[]) {
  LHAPDF::setVerbosity(0);
  NNPDF::SetVerbosity(0);

  if (argc < 2)
  {
    cerr << "Usage: fkgen <configuration file>"<<std::endl;
    exit(-1);
  }
  // Read configuration file
  libconfig::Config dPDFconfig;
  dPDFconfig.readFile(argv[1]);

  const int thID = dPDFconfig.lookup("qcd.THID");
  std::map<std::string, std::string> thMap;
  NNPDF::IndexDB theoryDB("data/theory.db", "theoryIndex");
  theoryDB.ExtractMap(thID, APFEL::kValues, thMap);
  APFEL::SetParam(thMap);


  NNPDF::CommonData data = CommonData::ReadFile("data/commondata/DATA_SLACD.dat", "data/commondata/systypes/SYSTYPE_SLACD_DEFAULT.dat");
  APFEL::SetQLimits( 1.0, getQ2Max_CommonData(data) );
  initEvolgrid(20, getXmin_CommonData(data)/2.0);
  APFEL::InitializeAPFEL_DIS();
  APFEL::EvolveAPFEL(1,1);

  std::stringstream description;
  description << "***************************" << std::endl << "SLACD" << std::endl<< "***************************";

  NNPDF::FKHeader FKhead;
  FKhead.AddTag(FKHeader::BLOB, "GridDesc", description.str());
  FKhead.AddTag(FKHeader::GRIDINFO, "SETNAME", data.GetSetName());
  FKhead.AddTag(FKHeader::GRIDINFO, "NDATA", data.GetNData());
  FKhead.AddTag(FKHeader::GRIDINFO, "HADRONIC", false );
  FKhead.AddTag(FKHeader::VERSIONS, "APFEL", APFEL::GetVersion());
  FKhead.AddTag(FKHeader::VERSIONS, "libnnpdf", NNPDF::getVersion());
  FKhead.AddTag(FKHeader::GRIDINFO, "NX", APFEL::nIntervals());
  stringstream xGheader;
  for (int i=0; i<APFEL::nIntervals(); i++)
    xGheader << std::setprecision(16) << std::scientific << APFEL::xGrid(i) <<std::endl;
  FKhead.AddTag(FKHeader::BLOB, "xGrid", xGheader.str());
  for (auto imap = thMap.begin(); imap != thMap.end(); imap++)
    FKhead.AddTag(FKHeader::THEORYINFO, imap->first, imap->second);

  FKhead.ResetFlavourMap();
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );

  for (size_t d=0; d<data.GetNData(); d++)
  {
    std::cout << d <<"/" <<data.GetNData() <<std::endl;
    const double x = data.GetKinematics(d,0)/2.0;
    const double Q = sqrt(data.GetKinematics(d,1));
    const double y = data.GetKinematics(d,2);

    for(size_t ix=0; ix<FK->GetNx(); ix++) 
      for(int ifl=0; ifl<14; ifl++) 
        FK->Fill(d, ix, ifl, diskernel("DIS_F2P", x, Q, y, ifl, ix) );
  }

  FK->Finalise();
  FK->Print("FK_SLACD.dat");

  // // Make directories
  // mkdir((base_path+"/dat").c_str(), 0777);
  // mkdir((base_path+"/dat/systypes").c_str(), 0777);
  // mkdir((base_path+"/thr").c_str(), 0777);

  exit(0);
}