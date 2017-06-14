// proton.h
// Handling of proton 'c-factors'

#pragma once

#include <libconfig.h++>
#include "NNPDF/lhapdfset.h"
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

 	virtual void GetPDF(real const& x, real const& Q2, int const& n, real* pdf) const
 	{
 		LHAPDFSet::GetPDF(x,Q2,n,pdf);
 		pdf[EVLN_V3] = 0;
 		pdf[EVLN_T3] = 0;
 	};

};