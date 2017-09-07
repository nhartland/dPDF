// nPDF -nh 11/14

#pragma once

#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include <NNPDF/pdfset.h>

#include <libconfig.h++>


// Initialise nPDF datasets
void ReadPlots(libconfig::Config const& settings, std::vector<NNPDF::DataSet>&);
void ReadData (libconfig::Config const& settings, std::vector<NNPDF::Experiment>&);
void InitData (libconfig::Config const& settings, std::vector<NNPDF::Experiment>, std::vector<NNPDF::Experiment>&, std::vector<NNPDF::Experiment>&);

// Kinematical cuts
bool passKinCuts(libconfig::Config const& settings, NNPDF::DataSet const& set, int const& idat);

/* Note there is a fair bit of (perhaps unneccesary) copying
   going on in these functions. I'm trying to avoid returning pointers
   to skip over ownership problems - nh */

/* Function to return an experimental dataset */ 
NNPDF::DataSet LoadDataSet(libconfig::Setting const& set, libconfig::Config const& settings);

/* Function to perform the filtering of an experimental dataset */
NNPDF::DataSet FilterData(NNPDF::DataSet const& set, libconfig::Config const& settings);

// Function to calculate and set the t0 predictions vector in a DataSet
void SetT0(NNPDF::DataSet& set, NNPDF::PDFSet const& proton, NNPDF::PDFSet const& deuteron);
