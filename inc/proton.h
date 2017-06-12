// proton.h
// Handling of proton 'c-factors'

#pragma once

#include <libconfig.h++>

void GenerateNormalisations(libconfig::Setting const& set, libconfig::Config const& settings);
void ExportProton(libconfig::Config const& settings, std::ostream& os);