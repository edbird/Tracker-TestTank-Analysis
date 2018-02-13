#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <string>
#include <ostream>

//#include <RTypes.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>


////////////////////////////////////////////////////////////////////////////////
// HISTOGRAM DIFFERENTIATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

const Int_t differentiate_method_simple{0};
const Int_t differentiate_method_simple_with_smooth{1};
const Int_t differentiate_method_filter_simulation{2};
void histogram_differentiate(TH1F* const differential_histo, const TH1F* const histo, Int_t method);
void histogram_differentiate_simple(TH1F* const differential_histo, const TH1F* const histo);
void histogram_differentiate_simple_with_smooth(TH1F* const differential_histo, const TH1F* const histo);
void histogram_differentiate_filter_simulation(TH1F* const differential_histo, const TH1F* const histo);

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

// Fit function for cathode time
Double_t cathode_time_fitf(Double_t *x_, Double_t *par);

// Plasma propagation time
Double_t ppt_fitf(Double_t *x_, Double_t *par);

// Z Position Cathode Time
Double_t zpos_cathode_time_fitf(Double_t *x_, Double_t *par);

// Z Position Cathode Time Profile (non-fit)
Double_t zpos_cathode_time_profilef(Double_t x, Double_t mean, Double_t theta);

// Gaussian distribution
Double_t gaussian_fitf(Double_t *x_, Double_t *par);

// Double Gaussian distribution
Double_t double_gaussian_fitf(Double_t *x_, Double_t *par);


////////////////////////////////////////////////////////////////////////////////
// TCANVAS FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void canvas(TH1F* const histogram_, const std::string& canvas_name_, const std::string& draw_opt_ = "");

void canvas(TH2F* const histogram_, const std::string& canvas_name_);

void canvas_fit(TH1F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_);

void canvas_fit(TH2F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_);

void canvas_scale_fit(TH1F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_);

void canvas_scale_fit(TH2F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_);


////////////////////////////////////////////////////////////////////////////////
// PRINT FIT FUNCTION OUTPUT PARAMETERS
////////////////////////////////////////////////////////////////////////////////

void fit_param_print(std::ostream& os, TF1* func);

void fit_param_print(std::ostream& os, TF2* func);


////////////////////////////////////////////////////////////////////////////////
// WAVEFORM OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void waveform_print(TH1F* histo, Long64_t &canvas_name_counter, std::string &output_file_name, const std::string &output_file_directory);

#endif
