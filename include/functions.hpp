#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <string>

#include <RTypes.h>
#include <TH1F.h>
#include <TF1.h>

////////////////////////////////////////////////////////////////////////////////
// CUT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

bool cut(const Double_t val, const Double_t low, const Double_t high);

bool cut_l(const Double_t val, const Double_t low);

////////////////////////////////////////////////////////////////////////////////
// CONVERSION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

std::string int_to_string(const Long64_t value, int width);

////////////////////////////////////////////////////////////////////////////////
// FIT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

// 2d Gaussian with rotation
Double_t fitf(Double_t *x_, Double_t *par);

// Fit function for data feast t0 timestamp
// (Also for MC when MC of this function is implemented in Falaise)
Double_t feast_t0_fitf(Double_t *x_, Double_t *par);

////////////////////////////////////////////////////////////////////////////////
// PRINT FIT FUNCTION OUTPUT PARAMETERS
////////////////////////////////////////////////////////////////////////////////

void fit_param_print(std::ostream& os, TF1* func, int n_param);

void fit_param_print(std::ostream& os, TF2* func, int n_param);

////////////////////////////////////////////////////////////////////////////////
// WAVEFORM OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void waveform_print(TH1F* histo, Long64_t &canvas_name_counter, std::string &output_file_name, const std::string &output_file_directory);

#endif