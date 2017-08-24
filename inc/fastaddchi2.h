// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <NNPDF/common.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include "cmaes.h"

using NNPDF::PDFSet;
using NNPDF::LHAPDFSet;
using NNPDF::DataSet;
using NNPDF::Experiment;
using NNPDF::real;
using NNPDF::FKSet;
using std::vector;

class DeuteronSet;

// Fast methods for the computation of chi2s.
void FastAddChi2(const PDFSet*, const PDFSet*, const Experiment*, real* chi2);
void Convolute(const PDFSet*, const PDFSet*, const Experiment*, real *);
void ComputePredictions(const PDFSet*, const PDFSet*, const FKSet*, real*);
double ComputeMemberChi2(const PDFSet*, const PDFSet*, const int, std::vector<Experiment> const&);

struct ErfComputer : public CostComputer
{ 
public:
	ErfComputer(NNPDF::LHAPDFSet const& _proton, vector<Experiment> const& _exps):
	proton(_proton),
	exps(_exps) {};

	vector<real> operator()(vector<gsl_vector*>const& p) const {return ComputeErf(p);}; 
  	real         operator()(gsl_vector*)         	   const; 
private:
	LHAPDFSet 		   const& proton;
	vector<Experiment> const& exps;
	vector<real> ComputeErf(vector<gsl_vector*>const& p) const;
};