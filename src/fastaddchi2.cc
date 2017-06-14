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

using NNPDF::ThPredictions;

enum process_type {F2, F2R, DYR};
const std::map<std::string, process_type> process_map = {{"BCDMSD", F2}, {"SLACD", F2}, {"NMCPD", F2R}, {"DYE886R", DYR}};
void ComputePredictions(const PDFSet* proton, const PDFSet* deuteron, const FKSet* fkset, real * theory)
{
  auto iproc = process_map.find(fkset->GetDataName());
  if (iproc == process_map.end()) throw NNPDF::RuntimeException("ComputePredictions", "Cannot find process type for set " + fkset->GetDataName());
  const int nPDF = deuteron->GetMembers();
  const int nDAT = fkset->GetNDataFK();
  switch (iproc->second)
  {
    case DYR: {
      const ThPredictions SigPD(proton, deuteron, fkset->GetFK(0));
      const ThPredictions SigPP(proton, proton,   fkset->GetFK(0));
      for (size_t idat=0; idat<nDAT; idat++)
        for (size_t imem=0; imem<nPDF; imem++)
          theory[idat*nPDF + imem] = SigPD.GetObs(idat, imem) / SigPP.GetObsCV(idat);
      break; }
    case F2: {
      ThPredictions::Convolute(deuteron, fkset->GetFK(0),theory);
      break; }

    case F2R: {
      const ThPredictions F2d(deuteron, fkset->GetFK(0));
      const ThPredictions F2p(proton,   fkset->GetFK(1));
      for (size_t idat=0; idat<nDAT; idat++)
        for (size_t imem=0; imem<nPDF; imem++)
          theory[idat*nPDF + imem] = F2d.GetObs(idat, imem) / F2p.GetObsCV(idat);
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