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

using NNPDF::ThPredictions;

enum process_type {F2, F2R, DYR};
const std::map<std::string, process_type> process_map = {{"BCDMSD", F2}, {"SLACD", F2}, {"NMCPD", F2R}, {"DYE886R", DYR}};
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
      {
        const NNPDF::real norm = 2.0*SigPP.GetObsCV(idat);
        for (size_t imem=0; imem<nPDF; imem++)
          theory[idat*nPDF + imem] /= norm;
      }
      break; }
    case F2: {
      ThPredictions::Convolute(deuteron, fkset->GetFK(0),theory);
      for (size_t i=0; i<nPDF*nDAT; i++) theory[i] /= 2.0;
      break; }

    case F2R: {
      ThPredictions::Convolute(deuteron, fkset->GetFK(0),theory);
      const ThPredictions F2p(proton,   fkset->GetFK(1));
      for (size_t idat=0; idat<nDAT; idat++)
      {
        const NNPDF::real norm = 2.0*F2p.GetObsCV(idat);
        for (size_t imem=0; imem<nPDF; imem++)
          theory[idat*nPDF + imem] /= norm;
      }
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