// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "fastaddchi2.h"
#include <NNPDF/thpredictions.h>
#include <NNPDF/chisquared.h>
#include <NNPDF/exceptions.h>
#include <NNPDF/utils.h>
#include "deuteronset.h"

using NNPDF::ThPredictions;

enum process_type {F2, F2R, DYR};
const std::map<std::string, process_type> process_map = {{"BCDMSD", F2},
                                                         {"SLACD", F2}, 
                                                         {"NMCPD_D", F2R}, 
                                                         {"DYE886R", DYR},
                                                         {"F2R1", F2},
                                                         {"F2R10", F2},
                                                         {"F2R100", F2},
                                                         {"F2R1000", F2}
                                                       };
void ComputePredictions(const PDFSet* proton, const PDFSet* deuteron, const FKSet* fkset, real * theory)
{
  auto iproc = process_map.find(fkset->GetDataName());
  if (iproc == process_map.end()) throw NNPDF::RuntimeException("ComputePredictions", "Cannot find process type for set " + fkset->GetDataName());
  const int nPDF = deuteron->GetMembers();
  const int nDAT = fkset->GetNDataFK();
  // Note: All data is given as F2 per nucleon (hence all the /2.0 s)
  switch (iproc->second)
  {
    case DYR: {
      ThPredictions::Convolute( proton, deuteron, fkset->GetFK(0), theory);
      const ThPredictions SigPP(proton, proton,   fkset->GetFK(0));
      for (size_t idat=0; idat<nDAT; idat++)
        for (size_t imem=0; imem<nPDF; imem++)
          theory[idat*nPDF + imem] /= 2.0*SigPP.GetObs()[idat*nPDF + imem];
      break; }
    case F2: {
      ThPredictions::Convolute(deuteron, fkset->GetFK(0),theory);
      for (size_t i=0; i<nPDF*nDAT; i++) theory[i] /= 2.0;
      break; }

    case F2R: {
      ThPredictions::Convolute(deuteron, fkset->GetFK(0),theory);
      const ThPredictions F2p(proton,   fkset->GetFK(1));
      for (size_t idat=0; idat<nDAT; idat++)
        for (size_t imem=0; imem<nPDF; imem++)
          theory[idat*nPDF + imem] /= 2.0*F2p.GetObs()[idat*nPDF + imem];
      break; }
  }
}

void Convolute(const PDFSet* proton, const PDFSet* deuteron, const Experiment* exp, real * theory)
{
  int index = 0;
  for (int i = 0; i < exp->GetNSet(); i++)
    {
      ComputePredictions(proton, deuteron, &exp->GetSet(i), theory+index);
      index += deuteron->GetMembers()*exp->GetSet(i).GetNData();
    }
}

void FastAddChi2(const PDFSet* proton, const PDFSet* deuteron, const Experiment* exp, real* chi2)
{
  // Set up theory array
  const int nMem = deuteron->GetMembers();
  real *theory = new real[exp->GetNData()*nMem]();

  // Perform convolution and chi^2 calculation
  Convolute(proton, deuteron, exp, theory);
  NNPDF::ComputeChi2(exp,nMem,theory,chi2);

  delete[] theory;
}

vector<real> ErfComputer::ComputeErf(vector<gsl_vector*>const& parameters) const
{
  DeuteronSet mutants(parameters);
  std::vector<real> Chi2Mem(mutants.GetMembers(), 0);
  for (auto exp : exps)
    FastAddChi2(&proton, &mutants, &exp, Chi2Mem.data());
  for (real& Chi2 : Chi2Mem) // Check for anomalous chi2 values (probably due to integration failure)
    if (Chi2 >= 1E20 || std::isnan(Chi2) || std::isinf(Chi2))
      Chi2 = std::numeric_limits<real>::infinity();
  return Chi2Mem;
}

real ErfComputer::operator()(gsl_vector* par) const
{
  const vector<gsl_vector*> parameters = {par}; 
  return ComputeErf(parameters)[0];
}