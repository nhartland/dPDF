// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <NNPDF/common.h>
#include <NNPDF/pdfset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
using NNPDF::PDFSet;
using NNPDF::DataSet;
using NNPDF::Experiment;
using NNPDF::real;
using NNPDF::FKSet;

// Fast methods for the computation of chi2s.
void FastAddChi2(const PDFSet*, const PDFSet*, const Experiment*, real* chi2);
void Convolute(const PDFSet*, const PDFSet*, const Experiment*, real *);
void ComputePredictions(const PDFSet*, const PDFSet*, const FKSet*, real*);
double ComputeMemberChi2(const PDFSet*, const PDFSet*, const int, std::vector<Experiment> const&);
