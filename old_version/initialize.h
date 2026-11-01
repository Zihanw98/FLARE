#pragma once
#include <string>

void ReadParameters(const std::string& filename);
// Empty filename retains the original hardcoded defaults.
void TempInit(const std::string& parameter_file = "");
void ApplyBoundaryConditions();