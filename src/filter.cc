// nPDF -nh 11/14

#include "filter.h"
#include "colour.h"
#include "fastaddchi2.h"
#include "deuteronset.h"

#include <NNPDF/fastkernel.h>
#include <NNPDF/pdfset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/randomgenerator.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>

using std::cout;
using std::endl;
using std::left;
using std::setw;

using NNPDF::sysError;

using namespace Colour;

bool passKinCuts(libconfig::Config const& settings, NNPDF::DataSet const& set, int const& idat)
{
  if (set.GetProc(idat).compare(0,3,std::string("DYP")) == 0)
  {
    const double DYinvMassMin = settings.lookup("cuts.DYMassMin");
    const double DYinvMassMax = settings.lookup("cuts.DYMassMax");
    const double invM  = sqrt(set.GetKinematics(idat,1));
    return ( invM > DYinvMassMin && invM < DYinvMassMax);
  }
  if (set.GetProc(idat).compare(0,3,std::string("DIS")) == 0)
  {
    const double x  = set.GetKinematics(idat,0);
    const double Q2 = set.GetKinematics(idat,1);
    const double W2 = Q2*(1-x)/x;
    double W2cut; settings.lookupValue("cuts.W2cut", W2cut);
    double Q2cut; settings.lookupValue("cuts.Q2cut", Q2cut);
    if (W2 <= W2cut) return false;
    if (Q2 <= Q2cut) return false;
  } 
  return true;
}

// Initialises all datasets
void ReadData(libconfig::Config const& settings, std::vector<NNPDF::Experiment>& outexp)
{
  cout<<endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  cout <<              "                        NUCLEAR NNPDF FILTER                         " <<endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  const bool vb = NNPDF::PDFSet::Verbose;
  NNPDF::PDFSet::Verbose = true;

  // Read T0 PDF set
  libconfig::Setting& fitsettings = settings.lookup("fit");
  NNPDF::LHAPDFSet t0set_proton(fitsettings["proton"],  NNPDF::LHAPDFSet::ER_MC);
  DeuteronSet      t0set_deuteron = DeuteronSet::ReadSet(fitsettings["t0set"], t0set_proton.GetMembers());

  NNPDF::PDFSet::Verbose = vb;

  cout <<              "                         T0 Set Initialised                          " <<endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  // Read datasets
  libconfig::Setting& sets = settings.lookup("data.sets");

  cout << setw(20) << left << "SETNAME" << "  "
       << setw(10) << left << "SYSTYPE" << "   "
       << setw(10) << left << "NDATA" << "  "
       << setw(10) << left << "NDATA_CUT" << "   "
       << endl;

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  std::map<std::string,NNPDF::DataSet> setMap;

  for (int i=0; i<sets.getLength(); i++)
  {
    std::string setname;
    int    systype = 0;

    libconfig::Setting& set = sets[i];
    set.lookupValue("name", setname);
    set.lookupValue("systype", systype);

    NNPDF::DataSet dset = LoadDataSet(set, settings);
    NNPDF::DataSet fset = FilterData(dset, settings);
    SetT0(fset, t0set_proton, t0set_deuteron);

    cout << setw(20) << left << setname << "  "
         << setw(10) << left << systype << "   "
         << setw(10) << left << dset.GetNData() << "  "
         << setw(10) << left << fset.GetNData() << "   "
         << endl;

    // Add set to vector
    setMap.insert(std::pair<std::string, NNPDF::DataSet>(setname,fset));
  }

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;  

  // Read datasets
  libconfig::Setting& exps = settings.lookup("data.experiments");

  cout << setw(20) << left << "EXPERIMENT" << "  "
       << setw(10) << left << "T/V SPLIT" << "   "
       << setw(10) << left << "NDATA" << "   "
       << endl;

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

  for (int i=0; i<exps.getLength(); i++)
  {
    libconfig::Setting& exp = exps[i];

    std::string expname;
    exp.lookupValue("name", expname);

    double tvsplit = 0;
    exp.lookupValue("tvsplit", tvsplit);

    // Fetch constituent sets
    std::vector<NNPDF::DataSet> subSets;
    for (int i=0; i<exp["sets"].getLength(); i++)
    {
      const char *setname = exp["sets"][i];
      subSets.push_back(setMap.find(setname)->second);
    }

    // Push experiment
    outexp.emplace_back(subSets,expname);

    cout << setw(20) << left << expname << "  "
         << setw(10) << left << tvsplit << "   "
         << setw(10) << left << outexp[outexp.size()-1].GetNData() << "  "
         << endl;
  }


  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  cout << FG_GREEN  << "                         FILTERING COMPLETE                          "<<FG_DEFAULT <<endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  
  cout <<endl;
  cout <<endl;
}

/* Initialises datasets and performs tr-valid split */
void InitData(libconfig::Config const& settings, std::vector<NNPDF::Experiment> experiments, std::vector<NNPDF::Experiment>& trainexp, std::vector<NNPDF::Experiment>& validexp)
{
  // Make artificial data replica
  for (size_t i=0; i<experiments.size(); i++)
    experiments[i].MakeReplica();

  cout<<endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
  cout <<              "                          TR-VAL SPLITTING                           " <<endl;
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;

    cout << setw(20) << left << "DATASET" << "  "
       << setw(10) << left << "N_TRAIN" << "   "
       << setw(10) << left << "N_VALID" << "   "
       << endl;

  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;


  libconfig::Setting& exps = settings.lookup("data.experiments");
  for (int i=0; i<exps.getLength(); i++)
  {
    std::string expname;
    double tvsplit = 0;

    libconfig::Setting& exp = exps[i];
    exp.lookupValue("name", expname);
    exp.lookupValue("tvsplit", tvsplit);

    if (expname.compare(experiments[i].GetExpName()) != 0)
    {
      std::cerr << "InitData Error: experiment names do not match: " <<expname <<" vs "<<experiments[i].GetExpName()<<std::endl;
      exit(-1);
    }

    cout << "EXPERIMENT: " << expname <<endl;
    cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;


    std::vector<NNPDF::DataSet> trainingSets;
    std::vector<NNPDF::DataSet> validationSets;

    for (int j=0; j<experiments[i].GetNSet(); j++)
    {
      std::vector<int> trainingMask;
      std::vector<int> validationMask;

      NNPDF::DataSet const& currentSet = experiments[i].GetSet(j);
      for (int k=0; k<currentSet.GetNData(); k++)
      {
        if (NNPDF::RandomGenerator::GetRNG()->GetRandomUniform() < tvsplit)
          trainingMask.push_back(k);
        else
          validationMask.push_back(k);
      }

      cout  << setw(20) << left << currentSet.GetSetName() << "  "
            << setw(10) << left << trainingMask.size() << "   "
            << setw(10) << left << validationMask.size() << "  "
            << endl;

      if (trainingMask.size() > 0)
        trainingSets.emplace_back(currentSet, trainingMask);
      if (validationMask.size() > 0)
        validationSets.emplace_back(currentSet, validationMask);
    }

    if (trainingSets.size() > 0 )
      trainexp.emplace_back(trainingSets, expname);
    if (validationSets.size() > 0 )
    validexp.emplace_back(validationSets, expname);

    const int nDat_TRN = (trainingSets.size() == 0) ? 0 : trainexp[trainexp.size()-1].GetNData();
    const int nDat_VAL = (validationSets.size() == 0) ? 0 : validexp[validexp.size()-1].GetNData();

    cout << setw(20) << left << "TOTAL:" << "  "
       << setw(10) << left << nDat_TRN << "   "
       << setw(10) << left << nDat_VAL << "   "
       << endl;

  }
  cout << FG_YELLOW << "---------------------------------------------------------------------"<<FG_DEFAULT <<endl;
}


/* Function to return a dataset */ 
NNPDF::DataSet LoadDataSet(libconfig::Setting const& set, libconfig::Config const& settings)
{
    std::string setname;
    set.lookupValue("name", setname);

	  // Filenames
    std::stringstream datname;
    std::stringstream sysname;

    const std::string cdata_path = "./data/commondata";
    datname << cdata_path<<"/DATA_"+setname+".dat";
    sysname << cdata_path<<"/systypes/SYSTYPE_"+setname+"_DEFAULT.dat";

    // Read commondata file
    NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(datname.str(), sysname.str());

    // Read C-factors
    std::vector<std::string> cfactors;

    // Read FK table operator
    std::string fkopstring = "NULL";
	  set.lookupValue("operator", fkopstring);
    NNPDF::SigmaOp FKOperator = NNPDF::FKSet::parseOperator(fkopstring);

    const int theoryIndex = settings.lookup("qcd.THID");
    std::vector<NNPDF::FKTable*> fkTables;
    for (int i=0; i<set["tables"].getLength(); i++)
    {
    	const char *tablename = set["tables"][i];
    	std::stringstream fkpath;
    	fkpath <<  "theory_"<<theoryIndex<<"/FK_"<<tablename<<".dat";
    	fkTables.push_back(new NNPDF::FKTable(fkpath.str(), cfactors));
	  }

	// Note FKSet owns its FKTables, no need to delete them
	NNPDF::FKSet fk_set(FKOperator, fkTables);
	return NNPDF::DataSet(cd, fk_set);
}

NNPDF::DataSet FilterData(NNPDF::DataSet const& set, libconfig::Config const& settings)
{
    std::vector<int> datamask;
    for (int i=0; i<set.GetNData(); i++)
        if (passKinCuts(settings,set,i)) datamask.push_back(i);
    return NNPDF::DataSet(set,datamask);
}

void SetT0(NNPDF::DataSet& set, NNPDF::PDFSet const& proton, NNPDF::PDFSet const& deuteron )
{
    NNPDF::real* theory = new NNPDF::real[set.GetNData()*deuteron.GetMembers()];
    ComputePredictions(&proton, &deuteron, &set, theory);
    NNPDF::ThPredictions t0pred(&deuteron, &set, theory); delete[] theory;
    set.SetT0(t0pred);
}


