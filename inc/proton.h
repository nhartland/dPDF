// proton.h
// Handling of proton 'c-factors'

#pragma once

#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>

#include <libconfig.h++>

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

	os 	<<	"*******************************************************************************************" << std::endl
		<<	" SetName: "<< setname << std::endl
		<<	" Author:  dPDF" << std::endl
		<<	" Date:    Today" << std::endl
		<<	" CodesUsed: None " << std::endl
		<<	" TheoryInput: None" << std::endl
		<<	" PDFset: None" << std::endl
		<<	" Warnings: None" << std::endl
		<<	" ********************************************************************************************" << std::endl;
	for (int i=0; i< theory.GetNData(); i++)
		os << 1.0/theory.GetObs(i, replica) <<"  "<<0<<std::endl;
	os.close();
}