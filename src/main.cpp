
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <random>

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TGaxis.h>

////////////////////////////////////////////////////////////////////////////////
// CUT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

bool cut(const Double_t val, const Double_t low, const Double_t high)
{
    if(low <= val)
    {
        if(val <= high)
            return true;
    }
    return false;
}

bool cut_l(const Double_t val, const Double_t low)
{
    if(low <= val)
    {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////
// CONVERSION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

std::string int_to_string(const Long64_t value, int width)
{
    std::string value_string = std::to_string(value);
    int string_width = value_string.size();
    if(string_width < width)
    {
        std::string pad;
        for(int w{string_width}; w < width; ++ w)
        {
            pad.push_back('0');
        }
        return std::string(pad + value_string);
    }
    return value_string;
}

////////////////////////////////////////////////////////////////////////////////
// FIT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

Double_t fitf(Double_t *x_, Double_t *par)
{
    Double_t x{x_[0]};
    Double_t y{x_[1]};

    Double_t A{par[0]};
    Double_t x0{par[1]};
    Double_t y0{par[2]};
    Double_t sigma_x{par[3]};
    Double_t sigma_y{par[4]};
    Double_t theta{par[5]};
    
    Double_t c{std::cos(theta)};
    Double_t s{std::sin(theta)};
    
    Double_t s2{std::sin(2.0 * theta)};
    
    Double_t cc{std::pow(c, 2.0)};
    Double_t ss{std::pow(s, 2.0)};
    
    Double_t sigma_x_2{std::pow(sigma_x, 2.0)};
    Double_t sigma_y_2{std::pow(sigma_y, 2.0)};
    
    Double_t a_{cc / (2.0 * sigma_x_2) + ss / (2.0 * sigma_y_2)};
    Double_t b_{(s2 / 2.0) * (1.0 / sigma_y_2 - 1.0 / sigma_x_2)};
    Double_t c_{ss / (2.0 * sigma_x_2) + cc / (2.0 * sigma_y_2)};
    
    return A * std::exp(-( a_*pow(x-x0, 2.0) + b_*(x-x0)*(y-y0) + c_*pow(y-y0, 2.0) ));
}

// Fit function for data feast t0 timestamp
// (Also for MC when MC of this function is implemented in Falaise)
Double_t feast_t0_fitf(Double_t *x_, Double_t *par)
{
    // variables
    Double_t x{x_[0]};

    // parameters
    Double_t A{par[0]}; // amplitude / normalization
    Double_t a{par[1]}; // start ramp up
    Double_t b{par[2]}; // end ramp up
    Double_t c{par[3]}; // start ramp down
    Double_t d{par[4]}; // end ramp down

    if(x <= a)
    {
        return 0.0;
    }
    else if(a < x && x <= b)
    {
        return A * (x - a) / (b - a);
    }
    else if(b < x && x <= c)
    {
        return A;
    }
    else if(c < x && x <= d)
    {
        return A * (1 - (x - c) / (d - c));
    }
    else
    {
        return 0.0;
    }


    std::cerr << "Warning, should not reach this line" << std::endl;
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// WAVEFORM OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void waveform_print(TH1F* histo, Long64_t &canvas_name_counter, std::string &output_file_name, const std::string &output_file_directory)
{

    if(histo != nullptr)
    {
        // print out the waveform of the t0 events for inspection
        //const std::string canvas_name_base("cathode_event_");
        const std::string canvas_name_base(output_file_name);
        std::string canvas_name(canvas_name_base + int_to_string(canvas_name_counter, 6));
        TCanvas *c = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);
        histo->Draw();
        output_file_name = std::string("./") + output_file_directory + std::string("/") + canvas_name + std::string(".png");
        c->SaveAs(output_file_name.c_str());
        delete c;
        ++ canvas_name_counter;
    }

}

////////////////////////////////////////////////////////////////////////////////
// TEST TANK DATA STRUCTURE
////////////////////////////////////////////////////////////////////////////////

typedef struct TestTankStorage
{
    Float_t time = 0;
    Float_t delay = 0;
    Float_t delay_since_good_trigger = 0;
    Int_t duration = 0;
    Float_t plasma_propagation_time = 0;
    Bool_t good_trigger = 0;
    Bool_t prev_good_trigger = 0;
    Bool_t with_cathode = 0;
    Float_t anode_peak = 0;
    Float_t anode_time = 0;
    Float_t cathode_peak = 0;
    Float_t cathode_time = 0;
    Float_t position = 0;
    Float_t half_position = 0;
    Float_t stop1 = 0;
    Float_t stop1_peak = 0;
    Float_t stop1_type = 0;
    Float_t stop2 = 0;
    Float_t stop2_peak = 0;
    Float_t stop2_type = 0;
    Int_t stopA = 0;
    Float_t deriv_rms = 0;
    Float_t feast_t0 = 0;
    Float_t feast_t1 = 0;
    Float_t feast_t2 = 0;
    Float_t feast_t3 = 0;
    Float_t feast_t4 = 0;
    TH1F *anode_histo = (TH1F*)0;
    TH1F *deriv_histo = (TH1F*)0;
    TH1F *cathode_histo = (TH1F*)0;
};

////////////////////////////////////////////////////////////////////////////////
// PROGRAM ARGUMENTS HELP
////////////////////////////////////////////////////////////////////////////////

void print_general_help(std::ostream& os, char* argv[])
{
    os  << argv[0] << " program arguments:\n"\
        << "\n"\
        << "Usage " << argv[0] << " --falaise-mc-input FILENAME\n"\
        << "\n"\
        << "Mandatory arguments:\n"\
        << "    --falaise-mc-input FILENAME     Set input file to FILENAME.\n"\
        << "        This is used when switching from analysis of default\n"\
        << "        datafile (cell8.root, obtained from testtank, real data)\n"\
        << "        to monte-carlo file generated using falaise.\n"\
        << "    --input FILENAME                Alias for --falaise-mc-input\n"\
        << "\n";
}

////////////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

    bool falaise_mc_ = false; // data is falaise mc type if set to true

    std::string filename = "cell8.root"; // set to default for testtank data
                                         // (real data)
    for(int ix{1}; ix < argc; ++ ix)
    {
        std::string arg{std::string(argv[ix])};

        if(ix < argc)
        {
            if(arg == std::string("-h") || arg == std::string("--help"))
            {
                print_general_help(std::cout, argv);
                std::cout.flush();
            }
            else if(arg == std::string("--input") || arg == std::string("--falaise-mc-input"))
            {
                if(++ ix < argc)
                {
                    filename = std::string(argv[ix]);
                    falaise_mc_ = true;
                }
                else
                {
                    std::cout << "[ WARNING ] : " << arg << " FILENAME" << std::endl;
                    std::cout << "[ WARNING ] : expected string FILENAME following argument " << arg << std::endl;
                }
            }
            else
            {
                std::cerr << "[ WARNING ] : unrecognized argument " << arg << " at index " << ix << std::endl;
            }
        }
    }

    // TODO: bug here - if arguments fail then program runs on default data anyway

    TGaxis::SetMaxDigits(2);
    
    ////////////////////////////////////////////////////////////////////////////
    // DATA LOAD AND SAVE
    ////////////////////////////////////////////////////////////////////////////
    
    TFile *f = new TFile(filename.c_str()); //("cell8.root");
    f->cd();
    TTree *t = (TTree*)f->Get("histo");
    t->SetDirectory(f);

    // Local storage, testtank format (origional C++ analysis code - this code)
    TestTankStorage store;
    
    // Default values
    store.time = 0;
    store.delay = 0;
    store.delay_since_good_trigger = 0;
    store.duration = 0;
    store.plasma_propagation_time = 0;
    store.good_trigger = 0;
    store.prev_good_trigger = 0;
    store.with_cathode = 0;
    store.anode_peak = 0;
    store.anode_time = 0;
    store.cathode_peak = 0;
    store.cathode_time = 0;
    store.position = 0;
    store.half_position = 0;
    store.stop1 = 0;
    store.stop1_peak = 0;
    store.stop1_type = 0;
    store.stop2 = 0;
    store.stop2_peak = 0;
    store.stop2_type = 0;
    store.stopA = 0;
    store.deriv_rms = 0;
    store.feast_t0 = 0;
    store.feast_t1 = 0;
    store.feast_t2 = 0;
    store.feast_t3 = 0;
    store.feast_t4 = 0;
    store.anode_histo = (TH1F*)0;
    store.deriv_histo = (TH1F*)0;
    store.cathode_histo = (TH1F*)0;
    
    
    t->SetBranchAddress("time", &store.time);
    t->SetBranchAddress("delay", &store.delay);
    t->SetBranchAddress("delay_since_good_trigger", &store.delay_since_good_trigger);
    t->SetBranchAddress("duration", &store.duration);
    t->SetBranchAddress("plasma_propagation_time", &store.plasma_propagation_time);
    t->SetBranchAddress("good_trigger", &store.good_trigger);
    t->SetBranchAddress("prev_good_trigger", &store.prev_good_trigger);
    t->SetBranchAddress("with_cathode", &store.with_cathode);
    t->SetBranchAddress("anode_peak", &store.anode_peak);
    t->SetBranchAddress("anode_time", &store.anode_time);
    t->SetBranchAddress("cathode_peak", &store.cathode_peak);
    t->SetBranchAddress("cathode_time", &store.cathode_time);
    t->SetBranchAddress("position", &store.position);
    t->SetBranchAddress("half_position", &store.half_position);
    t->SetBranchAddress("stop1", &store.stop1);
    t->SetBranchAddress("stop1_peak", &store.stop1_peak);
    t->SetBranchAddress("stop1_type", &store.stop1_type);
    t->SetBranchAddress("stop2", &store.stop2);
    t->SetBranchAddress("stop2_peak", &store.stop2_peak);
    t->SetBranchAddress("stop2_type", &store.stop2_type);
    t->SetBranchAddress("stopA", &store.stopA);
    t->SetBranchAddress("deriv_rms", &store.deriv_rms);
    t->SetBranchAddress("feast_t0", &store.feast_t0);
    t->SetBranchAddress("feast_t1", &store.feast_t1);
    t->SetBranchAddress("feast_t2", &store.feast_t2);
    t->SetBranchAddress("feast_t3", &store.feast_t3);
    t->SetBranchAddress("feast_t4", &store.feast_t4);
    //t->SetBranchAddress("anode_histo", &anode_histo);
    //t->SetBranchAddress("deriv_histo", &deriv_histo);
    //t->SetBranchAddress("cathode_histo", &cathode_histo);
    t->SetBranchAddress("anode_histo", &store.anode_histo);
    t->SetBranchAddress("deriv_histo", &store.deriv_histo);
    t->SetBranchAddress("cathode_histo", &store.cathode_histo);
    
    std::string filename_out{filename.substr(0, filename.rfind(".")) + std::string("_out") + filename.substr(filename.rfind("."))};
    TFile *f2 = new TFile(filename_out.c_str(), "RECREATE"); //("cell8_out.root", "RECREATE");
    f2->cd();
    TTree *t2 = new TTree("histo", "");
    t2->SetDirectory(f2);

    t2->Branch("time", &store.time);
    t2->Branch("delay", &store.delay);
    t2->Branch("delay_since_good_trigger", &store.delay_since_good_trigger);
    t2->Branch("duration", &store.duration);
    t2->Branch("plasma_propagation_time", &store.plasma_propagation_time);
    t2->Branch("good_trigger", &store.good_trigger);
    t2->Branch("prev_good_trigger", &store.prev_good_trigger);
    t2->Branch("with_cathode", &store.with_cathode);
    t2->Branch("anode_peak", &store.anode_peak);
    t2->Branch("anode_time", &store.anode_time);
    t2->Branch("cathode_peak", &store.cathode_peak);
    t2->Branch("cathode_time", &store.cathode_time);
    t2->Branch("position", &store.position);
    t2->Branch("half_position", &store.half_position);
    t2->Branch("stop1", &store.stop1);
    t2->Branch("stop1_peak", &store.stop1_peak);
    t2->Branch("stop1_type", &store.stop1_type);
    t2->Branch("stop2", &store.stop2);
    t2->Branch("stop2_peak", &store.stop2_peak);
    t2->Branch("stop2_type", &store.stop2_type);
    t2->Branch("stopA", &store.stopA);
    t2->Branch("deriv_rms", &store.deriv_rms);
    t2->Branch("feast_t0", &store.feast_t0);
    t2->Branch("feast_t1", &store.feast_t1);
    t2->Branch("feast_t2", &store.feast_t2);
    t2->Branch("feast_t3", &store.feast_t3);
    t2->Branch("feast_t4", &store.feast_t4);
    t2->Branch("anode_histo", &store.anode_histo);
    t2->Branch("deriv_histo", &store.deriv_histo);
    t2->Branch("cathode_histo", &store.cathode_histo);
    
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS GENERIC
    ////////////////////////////////////////////////////////////////////////////

    TH1F *h_plasma_propagation_time = new TH1F("h_plasma_propagation_time", "h_plasma_propagation_time", 50, 0.0, 500.0); //35.0, 45.0);
    h_plasma_propagation_time->SetStats(0);
    
    TH1F *h_cathode_time = new TH1F("h_cathode_time", "h_cathode_time", 50, -100.0, 100.0); //-50, 50.0);
    h_cathode_time->SetStats(0);
    
    //TH1F *h_feast_t0 = new TH1F("h_feast_t0", "h_feast_t0", 50, -1.0, 7.0); //4.76, 4.86);
    TH1F *h_feast_t0 = new TH1F("h_feast_t0", "h_feast_t0", 50, 4.76, 4.86);
    TH1F *h_feast_t1 = new TH1F("h_feast_t1", "h_feast_t1", 50, -100.0, 100.0); //0.0, 50.0);
    TH1F *h_feast_t2 = new TH1F("h_feast_t2", "h_feast_t2", 50, -100.0, 100.0); //0.0, 50.0);
    TH1F *h_feast_t3 = new TH1F("h_feast_t3", "h_feast_t3", 50, -100.0, 100.0); //0.0, 50.0);
    TH1F *h_feast_t4 = new TH1F("h_feast_t4", "h_feast_t4", 50, -100.0, 100.0); //0.0, 50.0); 
    
    h_feast_t0->SetStats(0);
    h_feast_t1->SetStats(0);
    h_feast_t2->SetStats(0);
    h_feast_t3->SetStats(0);
    h_feast_t4->SetStats(0);
    
    TH1F *h_feast_t0_diff = new TH1F("h_feast_t0_diff", "h_feast_t0_diff", 50, -50.0, 50.0);
    TH1F *h_feast_t1_diff = new TH1F("h_feast_t1_diff", "h_feast_t1_diff", 50, -50.0, 50.0);
    TH1F *h_feast_t2_diff = new TH1F("h_feast_t2_diff", "h_feast_t2_diff", 50, -50.0, 50.0);
    TH1F *h_feast_t3_diff = new TH1F("h_feast_t3_diff", "h_feast_t3_diff", 50, -50.0, 50.0);
    TH1F *h_feast_t4_diff = new TH1F("h_feast_t4_diff", "h_feast_t4_diff", 50, -50.0, 50.0); 
    
    h_feast_t0_diff->SetStats(0);
    h_feast_t1_diff->SetStats(0);
    h_feast_t2_diff->SetStats(0);
    h_feast_t3_diff->SetStats(0);
    h_feast_t4_diff->SetStats(0);
   
    ////////////////////////////////////////////////////////////////////////////
    // FIT FUNCTION FOR FEAST T0
    ////////////////////////////////////////////////////////////////////////////

    TF1* f_feast_t0 = new TF1("f_feast_t0", feast_t0_fitf, 4.76, 4.86, 5);
    f_feast_t0->SetParameter(0, 250.0);
    f_feast_t0->SetParameter(1, 4.780);
    f_feast_t0->SetParameter(2, 4.790);
    f_feast_t0->SetParameter(3, 4.829);
    f_feast_t0->SetParameter(4, 4.836);

    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS AND FIT FUNCTIONS (METHOD 1)
    ////////////////////////////////////////////////////////////////////////////
    
    TH1F *h_t0_smallest = new TH1F("h_t0_smallest", "h_t0_smallest", 50, -5.0, 3.0);
    TH1F *h_t_smallest = new TH1F("h_t_smallest", "h_t_smallest", 50, -0.4, 0.0); //4.3, 4.9);
    TH1F *h_t_next_smallest = new TH1F("h_t_next_smallest", "h_t_next_smallest", 50, 0.0, 0.4); //4.8, 5.2);
    
    TH1F *h_t_smallest_residual = new TH1F("h_t_smallest_residual", "h_t_smallest_residual", 50, -0.4, 0.0);
    TH1F *h_t_next_smallest_residual = new TH1F("h_t_next_smallest_residual", "h_t_next_smallest_residual", 50, 0.0, 0.4);
   
    h_t0_smallest->SetStats(0);
    h_t_smallest->SetStats(0);
    h_t_next_smallest->SetStats(0);

    h_t_smallest_residual->SetStats(0);
    h_t_next_smallest_residual->SetStats(0);

    TF1 *f_t0_smallest = new TF1("f_t0_smallest", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", -5.0, 3.0);
    f_t0_smallest->SetParameter(0, 20.0);
    f_t0_smallest->SetParameter(1, 0.0);
    f_t0_smallest->SetParameter(2, -1.0);
    f_t0_smallest->SetParameter(3, 0.0);
    f_t0_smallest->SetParameter(4, 0.0);
    
    TF1 *f_t_smallest = new TF1("f_t_smallest", "[0]*exp(-pow((x-[1])/[2], 2.0))", 4.3, 4.9);
    f_t_smallest->SetParameter(0, 300.0);
    f_t_smallest->SetParameter(1, -0.2);
    f_t_smallest->SetParameter(2, 0.05);
    
    TF1 *f_t_next_smallest = new TF1("f_t_next_smallest", "[0]*exp(-pow((x-[1])/[2], 2.0))", 4.8, 5.2);
    f_t_next_smallest->SetParameter(0, 350.0);
    f_t_next_smallest->SetParameter(1, 0.2);
    f_t_next_smallest->SetParameter(2, 0.05);
    
    TH2F *h_t_correlation = new TH2F("h_t_correlation", "h_t_correlation", 50, -0.4, 0.0, 50, 0.0, 0.4);
    TH2F *h_t_correlation_residual = new TH2F("h_t_correlation_residual", "h_t_correlation_residual", 50, -0.4, 0.0, 50, 0.0, 0.4);
    
    TH2F *h_t_correlation_mc = new TH2F("h_t_correlation_mc", "h_t_correlation_mc", 50, -0.4, 0.0, 50, 0.0, 0.4);
    
    h_t_correlation->SetStats(0);
    h_t_correlation_residual->SetStats(0);
    h_t_correlation_mc->SetStats(0);
    
    /*
    TF2 *f_t_correlation = new TF2("f_t_correlation", "[0]*exp(-([3]*pow(x-[1], 2.0) + 2.0*[4]*(x-[1])*(y-[2]) + [5]*pow(y-[2], 2.0)))", -0.4, 0.0, 0.0, 0.4);
    f_t_correlation->SetParameter(0, 1.0);
    f_t_correlation->SetParameter(1, -0.2);
    f_t_correlation->SetParameter(2, 0.2);
    f_t_correlation->SetParameter(3, 0.05);
    f_t_correlation->SetParameter(4, 0.0);
    f_t_correlation->SetParameter(5, 0.05);
    f_t_correlation->SetNpy(1000);
    f_t_correlation->SetContour(10);
    */
    
    TF2 *f_t_correlation = new TF2("f_t_correlation", fitf, -0.4, 0.0, 0.0, 0.4, 6, 2);
    
    f_t_correlation->SetParameter(0, 30.0);
    f_t_correlation->SetParameter(1, -0.2);
    f_t_correlation->SetParameter(2, 0.2);
    f_t_correlation->SetParameter(3, 0.05);
    f_t_correlation->SetParameter(4, 0.02);
    f_t_correlation->SetParameter(5, -0.3);
    
    f_t_correlation->SetParLimits(0, 0.0, 1000.0);
    f_t_correlation->SetParLimits(1, -1.0, 1.0);
    f_t_correlation->SetParLimits(2, -1.0, 1.0);
    f_t_correlation->SetParLimits(3, 0.0, 1000.0);
    f_t_correlation->SetParLimits(4, 0.0, 1000.0);
    f_t_correlation->SetParLimits(5, -1.6, 1.6);
    
    f_t_correlation->SetNpy(1000);
    f_t_correlation->SetContour(10);
    
    
    TF2 *f_t_correlation_mc = new TF2("f_t_correlation_mc", fitf, -0.4, 0.0, 0.0, 0.4, 6, 2);
    
    f_t_correlation_mc->SetParameter(0, 30.0);
    f_t_correlation_mc->SetParameter(1, -0.2);
    f_t_correlation_mc->SetParameter(2, 0.2);
    f_t_correlation_mc->SetParameter(3, 0.05);
    f_t_correlation_mc->SetParameter(4, 0.02);
    f_t_correlation_mc->SetParameter(5, -0.3);
    
    f_t_correlation_mc->SetParLimits(0, 0.0, 1000.0);
    f_t_correlation_mc->SetParLimits(1, -1.0, 1.0);
    f_t_correlation_mc->SetParLimits(2, -1.0, 1.0);
    f_t_correlation_mc->SetParLimits(3, 0.0, 1000.0);
    f_t_correlation_mc->SetParLimits(4, 0.0, 1000.0);
    f_t_correlation_mc->SetParLimits(5, -1.6, 1.6);
    
    f_t_correlation_mc->SetNpy(1000);
    f_t_correlation_mc->SetContour(10);
    
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS AND FIT FUNCTIONS (METHOD 2)
    ////////////////////////////////////////////////////////////////////////////
    
    TH1F *h_t_neg = new TH1F("h_t_neg", "h_t_neg", 50, -2.0, 2.0); //-0.4, 0.0);
    TH1F *h_t_pos = new TH1F("h_t_pos", "h_t_pos", 50, -2.0, 2.0); //0.0, 0.4);
    
    TH1F *h_t_neg_residual = new TH1F("h_t_neg_residual", "h_t_neg_residual", 50, -0.4, 0.0);
    TH1F *h_t_pos_residual = new TH1F("h_t_pos_residual", "h_t_pos_residual", 50, 0.0, 0.4);
   
    h_t_neg->SetStats(0);
    h_t_pos->SetStats(0);
    
    h_t_neg_residual->SetStats(0);
    h_t_pos_residual->SetStats(0);
    
    TF1 *f_t_neg = new TF1("f_t_neg", "[0]*exp(-pow((x-[1])/[2], 2.0))", 4.3, 4.9);
    f_t_neg->SetParameter(0, 300.0);
    f_t_neg->SetParameter(1, -0.2);
    f_t_neg->SetParameter(2, 0.05);
    
    TF1 *f_t_pos = new TF1("f_t_pos", "[0]*exp(-pow((x-[1])/[2], 2.0))", 4.8, 5.2);
    f_t_pos->SetParameter(0, 350.0);
    f_t_pos->SetParameter(1, 0.2);
    f_t_pos->SetParameter(2, 0.05);
    
    TH2F *h_t_cor = new TH2F("h_t_cor", "h_t_cor", 50, -0.4, 0.0, 50, 0.0, 0.4);
    TH2F *h_t_cor_residual = new TH2F("h_t_cor_residual", "h_t_cor_residual", 50, -0.4, 0.0, 50, 0.0, 0.4);
    
    TH2F *h_t_cor_mc = new TH2F("h_t_cor_mc", "h_t_cor_mc", 50, -0.4, 0.0, 50, 0.0, 0.4);
    
    h_t_cor->SetStats(0);
    h_t_cor_residual->SetStats(0);
    h_t_cor_mc->SetStats(0);
    
    /*
    TF2 *f_t_correlation = new TF2("f_t_correlation", "[0]*exp(-([3]*pow(x-[1], 2.0) + 2.0*[4]*(x-[1])*(y-[2]) + [5]*pow(y-[2], 2.0)))", -0.4, 0.0, 0.0, 0.4);
    f_t_correlation->SetParameter(0, 1.0);
    f_t_correlation->SetParameter(1, -0.2);
    f_t_correlation->SetParameter(2, 0.2);
    f_t_correlation->SetParameter(3, 0.05);
    f_t_correlation->SetParameter(4, 0.0);
    f_t_correlation->SetParameter(5, 0.05);
    f_t_correlation->SetNpy(1000);
    f_t_correlation->SetContour(10);
    */
    
    
    TF2 *f_t_cor = new TF2("f_t_cor", fitf, -0.4, 0.0, 0.0, 0.4, 6, 2);
    
    f_t_cor->SetParameter(0, 30.0);
    f_t_cor->SetParameter(1, -0.2);
    f_t_cor->SetParameter(2, 0.2);
    f_t_cor->SetParameter(3, 0.05);
    f_t_cor->SetParameter(4, 0.02);
    f_t_cor->SetParameter(5, -0.3);
    
    f_t_cor->SetParLimits(0, 0.0, 1000.0);
    f_t_cor->SetParLimits(1, -1.0, 1.0);
    f_t_cor->SetParLimits(2, -1.0, 1.0);
    f_t_cor->SetParLimits(3, 0.0, 1000.0);
    f_t_cor->SetParLimits(4, 0.0, 1000.0);
    f_t_cor->SetParLimits(5, -1.6, 1.6);
    
    f_t_cor->SetNpy(1000);
    f_t_cor->SetContour(10);
    
    
    TF2 *f_t_cor_mc = new TF2("f_t_cor_mc", fitf, -0.4, 0.0, 0.0, 0.4, 6, 2);
    
    f_t_cor_mc->SetParameter(0, 30.0);
    f_t_cor_mc->SetParameter(1, -0.2);
    f_t_cor_mc->SetParameter(2, 0.2);
    f_t_cor_mc->SetParameter(3, 0.05);
    f_t_cor_mc->SetParameter(4, 0.02);
    f_t_cor_mc->SetParameter(5, -0.3);
    
    f_t_cor_mc->SetParLimits(0, 0.0, 1000.0);
    f_t_cor_mc->SetParLimits(1, -1.0, 1.0);
    f_t_cor_mc->SetParLimits(2, -1.0, 1.0);
    f_t_cor_mc->SetParLimits(3, 0.0, 1000.0);
    f_t_cor_mc->SetParLimits(4, 0.0, 1000.0);
    f_t_cor_mc->SetParLimits(5, -1.6, 1.6);
    
    f_t_cor_mc->SetNpy(1000);
    f_t_cor_mc->SetContour(10);
    
    
    ////////////////////////////////////////////////////////////////////////////
    // DATA LOOP
    ////////////////////////////////////////////////////////////////////////////
   
    #define COUT_TIMESTAMP_GOOD 0
    #define COUT_TIMESTAMP_FAIL 0
    #define COUT_TIMESTAMP_WAIT 0
    #define WAVEFORM_PRINT_GOOD 0
    #define WAVEFORM_PRINT_FAIL 0

    Long64_t canvas_name_counter{0};
    
    std::vector<Float_t> width_vector;

    // 2018-01-18: Compute the mean anode time for the test-tank data
    // Required for falaise simulation
    // TODO: compute this using a histogram
    double mean_anode_time = 0.0;

    Long64_t count_accept{0};
    Long64_t count_reject{0};
    Long64_t count_accept_timestamp{0};
    Long64_t count_reject_timestamp{0};

    Long64_t byte_count = 0;
    Long64_t byte_count_entry = 0;
    Long64_t num_entry = t->GetEntries();
    //std::cout << "num_entry=" << num_entry << std::endl;
    for(Long64_t ix = 0; ix < num_entry; ++ ix)
    {    
        byte_count_entry = t->GetEntry(ix);
        if(byte_count_entry <= 0)
            continue;
        
        //std::cout << plasma_propagation_time << std::endl;
        //std::cin.get();
        h_plasma_propagation_time->Fill(store.plasma_propagation_time);
        
        // check anode and cathode peaks - only for real data
        bool event_check_flag{true};
        if(falaise_mc_ == false)
        {
            event_check_flag = cut_l(store.anode_peak, 250.0) && cut_l(store.cathode_peak, 350.0);
        }

        if
        (
            // if using real data, apply "real data" cuts
            (
                (falaise_mc_ == false) &&
                cut(store.plasma_propagation_time, 38.5, 42.0) &&
                cut(store.position, -0.95, 0.95) &&
                event_check_flag
            )
        ||
            // if using MC, apply the "MC cut" to PPT
            (
                (falaise_mc_ == true) &&
                cut(store.plasma_propagation_time, 57.0, 59.0) && // TODO: these are temp. test values
                cut(store.position, -0.95, 0.95) &&
                event_check_flag
            )
        )
        {
            ++ count_accept;

            mean_anode_time += store.feast_t0;
            
            //std::cout << "Good event" << std::endl;
            //std::cout << "Event index is: " << ix << std::endl;

            // TODO: shift all the timestamp variables so that they make more logical sense
            // this is currently done in the falaise module, however should be done here
            // instead

            #if COUT_TIMESTAMP_GOOD
                std::cout << "cathode               : " << store.cathode_time << "\n"\
                          << "t0                    : " << store.feast_t0 << "\n"\
                          << "t1                    : " << store.feast_t1 << "\n"\
                          << "t3                    : " << store.feast_t3 << "\n"\
                          << "t2                    : " << store.feast_t2 << "\n"\
                          << "t4                    : " << store.feast_t4 << "\n"\
                          << "cathode + t0          : " << store.cathode_time + store.feast_t0 << "\n"\
                          << "(t1 - t0) - cathode   : " << store.feast_t1 - (store.cathode_time + store.feast_t0) << "\n"\
                          << "t3 - (cathode + t0)   : " << store.feast_t3 - (store.cathode_time + store.feast_t0) << "\n"\
                          << "t2 - (cathode + t0)   : " << store.feast_t2 - (store.cathode_time + store.feast_t0) << "\n"\
                          << "t4 - (cathode + t0)   : " << store.feast_t4 - (store.cathode_time + store.feast_t0) << "\n";
                #if COUT_TIMESTAMP_WAIT
                    std::cin.get();
                #endif
            #endif
            
            //std::cout << "fill " << ix << std::endl;
            t2->Fill();
            
            h_cathode_time->Fill(store.cathode_time);
            
            h_feast_t0->Fill(store.feast_t0);
            h_feast_t1->Fill(store.feast_t1);
            h_feast_t2->Fill(store.feast_t2);
            h_feast_t3->Fill(store.feast_t3);
            h_feast_t4->Fill(store.feast_t4);
            
            h_feast_t0_diff->Fill(store.feast_t0 - store.cathode_time);
            h_feast_t1_diff->Fill(store.feast_t1 - store.cathode_time - store.feast_t0);
            h_feast_t2_diff->Fill(store.feast_t2 - store.cathode_time - store.feast_t0);
            h_feast_t3_diff->Fill(store.feast_t3 - store.cathode_time - store.feast_t0);
            h_feast_t4_diff->Fill(store.feast_t4 - store.cathode_time - store.feast_t0);
            
            // find a single, or pair of timestamps,
            // which have smallest abs difference to the
            // cathode timestamp
            // a single time is for t0,
            // a pair of times is for t other, where other
            // is a pair of 1,2,3,4
            
            // NOTE: Cathode times must be relative to feast_t0, however the
            // other timestamps are NOT relative to feast_t0
            // Why would anyone design it to be like this?
            // Seems like: Absolute time of event is triggered by anode rising
            // and other timestamp points on anode are triggered and recorded
            // using the same "absolute time" counter
            // BUT: that the triggering of the anode time (t0 event) starts a
            // new counter which is used for the cathode times
            // Is this correct?
            
            // Signed differences
            Double_t t0_diff = store.feast_t0 - (store.cathode_time /*+ feast_t0*/); // TODO: figure out what this is
            Double_t t1_diff = store.feast_t1 - (store.cathode_time + store.feast_t0);
            Double_t t2_diff = store.feast_t2 - (store.cathode_time + store.feast_t0);
            Double_t t3_diff = store.feast_t3 - (store.cathode_time + store.feast_t0);
            Double_t t4_diff = store.feast_t4 - (store.cathode_time + store.feast_t0); // TODO: maybe these are missleading
            
            // Absolute value differences
            Double_t t0_abs_diff = std::abs(store.feast_t0 - store.cathode_time); // TODO: maybe these are not missleading?
            Double_t t1_abs_diff = std::abs(store.feast_t1 - store.cathode_time);
            Double_t t2_abs_diff = std::abs(store.feast_t2 - store.cathode_time);
            Double_t t3_abs_diff = std::abs(store.feast_t3 - store.cathode_time);
            Double_t t4_abs_diff = std::abs(store.feast_t4 - store.cathode_time);
            
            // find smallest abs diff of all t 1/2/3/4
            std::vector<Double_t> sort_me;
            sort_me.push_back(t1_abs_diff);
            sort_me.push_back(t2_abs_diff);
            sort_me.push_back(t3_abs_diff);
            sort_me.push_back(t4_abs_diff);
            std::sort(sort_me.begin(), sort_me.end());
            Double_t small_abs_diff_1 = sort_me.at(0);
            Double_t small_abs_diff_2 = sort_me.at(1);
            
            // find nearest positive / negative value combination
            std::vector<Double_t> sort_me_neg;
            std::vector<Double_t> sort_me_pos;
            sort_me_neg.clear();
            sort_me_pos.clear();
            
            if(t1_diff < 0.0)
                sort_me_neg.push_back(t1_diff);
            else if(t1_diff > 0.0)
                sort_me_pos.push_back(t1_diff);
            
            if(t2_diff < 0.0)
                sort_me_neg.push_back(t2_diff);
            else if(t2_diff > 0.0)
                sort_me_pos.push_back(t2_diff);
            
            if(t3_diff < 0.0)
                sort_me_neg.push_back(t3_diff);
            else if(t3_diff > 0.0)
                sort_me_pos.push_back(t3_diff);
            
            if(t4_diff < 0.0)
                sort_me_neg.push_back(t4_diff);
            else if(t4_diff > 0.0)
                sort_me_pos.push_back(t4_diff);
            
            std::sort(sort_me_pos.begin(), sort_me_pos.end());
            std::sort(sort_me_neg.begin(), sort_me_neg.end());
            
            if((sort_me_pos.size() > 0) && (sort_me_neg.size() > 0))
            {
                //std::cout << *sort_me_pos.begin() << ", " << *sort_me_neg.rbegin() << std::endl;
                
                ++ count_accept_timestamp;

                h_t_neg->Fill(*sort_me_neg.rbegin());
                h_t_pos->Fill(*sort_me_pos.begin());
                
                h_t_cor->Fill(*sort_me_neg.rbegin(), *sort_me_pos.begin());
                
                std::string output_file_name("anode_");
                //canvas_name_counter = ix;
                Long64_t ix_copy{ix};
                #if WAVEFORM_PRINT_GOOD
                    waveform_print(store.anode_histo, ix_copy /*canvas_name_counter*/, output_file_name, "anode_histo_good");
                #endif
            }
            else
            {
                ++ count_reject_timestamp; 

                std::string output_file_name("anode_");
                //canvas_name_counter = ix;
                Long64_t ix_copy{ix};
                //waveform_print(ix_copy /*canvas_name_counter*/, cathode_histo, output_file_name);
                #if WAVEFORM_PRINT_FAIL
                    waveform_print(store.anode_histo, ix_copy /*canvas_name_counter*/, output_file_name, "anode_histo_fail");
                #endif
        
                // NOTE: This isn't the fail part of the first if statement
                std::cout << "Rejecting event with no timestamp either side" << std::endl;
                std::cout << "Event index is: " << ix << std::endl;
                std::cout << "Printing event to file: " << output_file_name << std::endl;
                std::cout << "sort_me_neg.size() = " << sort_me_neg.size() << std::endl;
                std::cout << "sort_me_pos.size() = " << sort_me_pos.size() << std::endl;
                
                #if COUT_TIMESTAMP_FAIL
                std::cout << "cathode               : " << store.cathode_time << "\n"\
                          << "t0                    : " << store.feast_t0 << "\n"\
                          << "t1                    : " << store.feast_t1 << "\n"\
                          << "t3                    : " << store.feast_t3 << "\n"\
                          << "t2                    : " << store.feast_t2 << "\n"\
                          << "t4                    : " << store.feast_t4 << "\n"\
                          << "cathode + t0          : " << store.cathode_time + store.feast_t0 << "\n"\
                          << "t1 - (cathode + t0)   : " << store.feast_t1 - (store.cathode_time + store.feast_t0) << "\n"\
                          << "t3 - (cathode + t0)   : " << store.feast_t3 - (store.cathode_time + store.feast_t0) << "\n"\
                          << "t2 - (cathode + t0)   : " << store.feast_t2 - (store.cathode_time + store.feast_t0) << "\n"\
                          << "t4 - (cathode + t0)   : " << store.feast_t4 - (store.cathode_time + store.feast_t0) << "\n";

                    #if COUT_TIMESTAMP_WAIT
                        std::cin.get();
                    #endif
                #endif
                
            }
            
            // is t0 diff smaller than all?
            if(t0_abs_diff < small_abs_diff_1)
            {
                // use only time t0
                h_t0_smallest->Fill(store.feast_t0 - store.cathode_time);
            }
            else if(std::abs(sort_me.at(2) - small_abs_diff_2) <= 3.0 * 0.065)
            {
                std::cout << "Rejecting event with close times" << std::endl;
            }
            else
            {
                // t0 may still be smaller than small_abs_diff_2
                // however assume it is not
                h_t_smallest->Fill(small_abs_diff_1 - store.feast_t0);
                h_t_next_smallest->Fill(small_abs_diff_2 - store.feast_t0);
                
                h_t_correlation->Fill(small_abs_diff_1 - store.feast_t0, small_abs_diff_2 - store.feast_t0);
                
                
                /*
                // print out the waveform of the t0 events for inspection
                const std::string canvas_name_base("cathode_event_");
                std::string canvas_name(canvas_name_base + int_to_string(canvas_name_counter, 6));
                TCanvas *c = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);
                cathode_histo->Draw();
                c->SaveAs((std::string("./t0_cathode_histo/") + canvas_name + std::string(".png")).c_str());
                delete c;
                ++ canvas_name_counter;
                
                Float_t maximum = cathode_histo->GetMaximum();
                Float_t maximum_pos = cathode_histo->GetBinCenter(cathode_histo->GetMaximumBin());
                
                Float_t half_maximum = 0.5 * maximum;
                Float_t half_maximum_pos_1, half_maximum_pos_2;
                int check_flag{0};
                for(int bin{cathode_histo->GetMaximumBin()}; bin >= 1; -- bin)
                {
                    if(cathode_histo->GetBinContent(bin) <= half_maximum)
                    {
                        half_maximum_pos_1 = cathode_histo->GetBinCenter(bin);
                        ++ check_flag;
                    }
                }
                for(int bin{cathode_histo->GetMaximumBin()}; bin <= cathode_histo->GetNbinsX(); ++ bin)
                {
                    if(cathode_histo->GetBinContent(bin) <= half_maximum)
                    {
                        half_maximum_pos_2 = cathode_histo->GetBinCenter(bin);
                        ++ check_flag;
                    }
                }
                if(check_flag == 2)
                {
                    Float_t width = half_maximum_pos_2 - half_maximum_pos_1;
                    width_vector.push_back(width);
                }
                */
            }
            
        }
        else
        {
            ++ count_reject;
        }
    }

    std::streamsize ss = std::cout.precision();
    std::cout.precision(10);
    mean_anode_time /= (double)count_accept;
    std::cout << "mean_anode_time = " << mean_anode_time << std::endl;
    std::cout.precision(ss);
    
    //std::cout << "I am done" << std::endl;
    
    t2->Write();
    //f2->Close();
    
    #define TIMESTAMP_CANVAS_ENABLE 1 // TODO:move
    #if TIMESTAMP_CANVAS_ENABLE
        TCanvas *c_plasma_propagation_time = new TCanvas("c_plasma_propagation_time", "c_plasma_propagation_time", 800, 600);
        h_plasma_propagation_time->Draw();
        c_plasma_propagation_time->SaveAs("c_plasma_propagation_time.C");
        c_plasma_propagation_time->SaveAs("c_plasma_propagation_time.png");
        c_plasma_propagation_time->SaveAs("c_plasma_propagation_time.pdf");
        h_plasma_propagation_time->Write();
        delete h_plasma_propagation_time;
        delete c_plasma_propagation_time;
        
        TCanvas *c_feast_t0 = new TCanvas("c_feast_t0", "c_feast_t0", 800, 600);
        h_feast_t0->Fit("f_feast_t0");
        h_feast_t0->Draw("E");
        c_feast_t0->SaveAs("c_feast_t0.C");
        c_feast_t0->SaveAs("c_feast_t0.png");
        c_feast_t0->SaveAs("c_feast_t0.pdf");
        h_feast_t0->Write(); 
        delete h_feast_t0;
        delete c_feast_t0;
        
        TCanvas *c_feast_t1 = new TCanvas("c_feast_t1", "c_feast_t1", 800, 600);
        h_feast_t1->Draw();
        c_feast_t1->SaveAs("c_feast_t1.C");
        c_feast_t1->SaveAs("c_feast_t1.png");
        c_feast_t1->SaveAs("c_feast_t1.pdf");
        h_feast_t1->Write();
        delete h_feast_t1;
        delete c_feast_t1;
        
        TCanvas *c_feast_t2 = new TCanvas("c_feast_t2", "c_feast_t2", 800, 600);
        h_feast_t2->Draw();
        c_feast_t2->SaveAs("c_feast_t2.C");
        c_feast_t2->SaveAs("c_feast_t2.png");
        c_feast_t2->SaveAs("c_feast_t2.pdf");
        h_feast_t2->Write();
        delete h_feast_t2;
        delete c_feast_t2;
        
        TCanvas *c_feast_t3 = new TCanvas("c_feast_t3", "c_feast_t3", 800, 600);
        h_feast_t3->Draw();
        c_feast_t3->SaveAs("c_feast_t3.C");
        c_feast_t3->SaveAs("c_feast_t3.png");
        c_feast_t3->SaveAs("c_feast_t3.pdf");
        h_feast_t3->Write();
        delete h_feast_t3;
        delete c_feast_t3;
            
        TCanvas *c_feast_t4 = new TCanvas("c_feast_t4", "c_feast_t4", 800, 600);
        h_feast_t4->Draw();
        c_feast_t4->SaveAs("c_feast_t4.C");
        c_feast_t4->SaveAs("c_feast_t4.png");
        c_feast_t4->SaveAs("c_feast_t4.pdf");
        h_feast_t4->Write();
        delete h_feast_t4;
        delete c_feast_t4;
    
        std::cout << "chisq = " << f_feast_t0->GetChisquare() / f_feast_t0->GetNDF() << std::endl;
        std::cout << "p0 = " << f_feast_t0->GetParameter(0) << " +- " << f_feast_t0->GetParError(0) << std::endl;
        std::cout << "p1 = " << f_feast_t0->GetParameter(1) << " +- " << f_feast_t0->GetParError(1) << std::endl;
        std::cout << "p2 = " << f_feast_t0->GetParameter(2) << " +- " << f_feast_t0->GetParError(2) << std::endl;
        std::cout << "p3 = " << f_feast_t0->GetParameter(3) << " +- " << f_feast_t0->GetParError(3) << std::endl;
        std::cout << "p4 = " << f_feast_t0->GetParameter(4) << " +- " << f_feast_t0->GetParError(4) << std::endl;

    #endif
    
    /*
    TCanvas *c_feast_t0_diff = new TCanvas("c_feast_t0_diff", "c_feast_t0_diff", 800, 600);
    h_feast_t0_diff->Draw();
    c_feast_t0_diff->SaveAs("c_feast_t0_diff.C");
    c_feast_t0_diff->SaveAs("c_feast_t0_diff.png");
    c_feast_t0_diff->SaveAs("c_feast_t0_diff.pdf");
        h_feast_t0_diff->Write();
    delete h_feast_t0_diff;
    delete c_feast_t0_diff;
    
    TCanvas *c_feast_t1_diff = new TCanvas("c_feast_t1_diff", "c_feast_t1_diff", 800, 600);
    h_feast_t1_diff->Draw();
    c_feast_t1_diff->SaveAs("c_feast_t1_diff.C");
    c_feast_t1_diff->SaveAs("c_feast_t1_diff.png");
    c_feast_t1_diff->SaveAs("c_feast_t1_diff.pdf");
        h_feast_t1_diff->Write();
    delete h_feast_t1_diff;
    delete c_feast_t1_diff;
    
    TCanvas *c_feast_t2_diff = new TCanvas("c_feast_t2_diff", "c_feast_t2_diff", 800, 600);
    h_feast_t2_diff->Draw();
    c_feast_t2_diff->SaveAs("c_feast_t2_diff.C");
    c_feast_t2_diff->SaveAs("c_feast_t2_diff.png");
    c_feast_t2_diff->SaveAs("c_feast_t2_diff.pdf");
        h_feast_t2_diff->Write();
    delete h_feast_t2_diff;
    delete c_feast_t2_diff;
    
    TCanvas *c_feast_t3_diff = new TCanvas("c_feast_t3_diff", "c_feast_t3_diff", 800, 600);
    h_feast_t3_diff->Draw();
    c_feast_t3_diff->SaveAs("c_feast_t3_diff.C");
    c_feast_t3_diff->SaveAs("c_feast_t3_diff.png");
    c_feast_t3_diff->SaveAs("c_feast_t3_diff.pdf");
        h_feast_t3_diff->Write();
    delete h_feast_t3_diff;
    delete c_feast_t3_diff;
        
    TCanvas *c_feast_t4_diff = new TCanvas("c_feast_t4_diff", "c_feast_t4_diff", 800, 600);
    h_feast_t4_diff->Draw();
    c_feast_t4_diff->SaveAs("c_feast_t4_diff.C");
    c_feast_t4_diff->SaveAs("c_feast_t4_diff.png");
    c_feast_t4_diff->SaveAs("c_feast_t4_diff.pdf");
        h_feast_t4_diff->Write();
    delete h_feast_t4_diff;
    delete c_feast_t4_diff;
    */
    
    ////////////////////////////////////////////////////////////////////////////
    // CANVAS OUTPUT (METHOD 1)
    ////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n>>> f_t0_smallest" << std::endl;
    TCanvas *c_t0_smallest = new TCanvas("c_t0_smallest", "c_t0_smallest", 800, 600);
    h_t0_smallest->Fit("f_t0_smallest");
    h_t0_smallest->Draw("E");
    c_t0_smallest->SaveAs("c_t0_smallest.C");
    c_t0_smallest->SaveAs("c_t0_smallest.png");
    c_t0_smallest->SaveAs("c_t0_smallest.pdf");
        h_t0_smallest->Write();
    delete h_t0_smallest;
    delete c_t0_smallest;
    
    std::cout << "\n>>> f_t_smallest" << std::endl;
    TCanvas *c_t_smallest = new TCanvas("c_t_smallest", "c_t_smallest", 800, 600);
    h_t_smallest->Fit("f_t_smallest");
    h_t_smallest->Draw("E");
    c_t_smallest->SaveAs("c_t_smallest.C");
    c_t_smallest->SaveAs("c_t_smallest.png");
    c_t_smallest->SaveAs("c_t_smallest.pdf");
        h_t_smallest->Write();
    //delete h_t_smallest;
    delete c_t_smallest;
    
    std::cout << "\n>>> f_t_next_smallest" << std::endl;
    TCanvas *c_t_next_smallest = new TCanvas("c_t_next_smallest", "c_t_next_smallest", 800, 600);
    h_t_next_smallest->Fit("f_t_next_smallest");
    h_t_next_smallest->Draw("E");
    c_t_next_smallest->SaveAs("c_t_next_smallest.C");
    c_t_next_smallest->SaveAs("c_t_next_smallest.png");
    c_t_next_smallest->SaveAs("c_t_next_smallest.pdf");
        h_t_next_smallest->Write();
    //delete h_t_next_smallest;
    delete c_t_next_smallest;
    
    std::cout << "chisq = " << f_t0_smallest->GetChisquare() / f_t0_smallest->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t0_smallest->GetParameter(0) << " +- " << f_t0_smallest->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t0_smallest->GetParameter(1) << " +- " << f_t0_smallest->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t0_smallest->GetParameter(2) << " +- " << f_t0_smallest->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t0_smallest->GetParameter(3) << " +- " << f_t0_smallest->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t0_smallest->GetParameter(4) << " +- " << f_t0_smallest->GetParError(4) << std::endl;
    
    std::cout << "chisq = " << f_t_smallest->GetChisquare() / f_t_smallest->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_smallest->GetParameter(0) << " +- " << f_t_smallest->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_smallest->GetParameter(1) << " +- " << f_t_smallest->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_smallest->GetParameter(2) << " +- " << f_t_smallest->GetParError(2) << std::endl;
    
    std::cout << "chisq=" << f_t_next_smallest->GetChisquare() / f_t_next_smallest->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_next_smallest->GetParameter(0) << " +- " << f_t_next_smallest->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_next_smallest->GetParameter(1) << " +- " << f_t_next_smallest->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_next_smallest->GetParameter(2) << " +- " << f_t_next_smallest->GetParError(2) << std::endl;
    
    // create residuals plot
    for(Int_t i = 1; i <= h_t_smallest->GetNbinsX(); ++ i)
    {
        Double_t func_x = h_t_smallest_residual->GetXaxis()->GetBinCenter(i);
        Double_t func_eval = f_t_smallest->Eval(func_x);
        
        h_t_smallest_residual->SetBinContent(i, h_t_smallest->GetBinContent(i) - func_eval);
        h_t_smallest_residual->SetBinError(i, h_t_smallest->GetBinError(i));
    }
    
    for(Int_t i = 1; i <= h_t_next_smallest->GetNbinsX(); ++ i)
    {
        Double_t func_x = h_t_next_smallest_residual->GetXaxis()->GetBinCenter(i);
        Double_t func_eval = f_t_next_smallest->Eval(func_x);
        
        h_t_next_smallest_residual->SetBinContent(i, h_t_next_smallest->GetBinContent(i) - func_eval);
        h_t_next_smallest_residual->SetBinError(i, h_t_next_smallest->GetBinError(i));
    }
    
    
    TCanvas *c_t_smallest_residual = new TCanvas("c_t_smallest_residual", "c_t_smallest_residual", 800, 600);
    h_t_smallest_residual->Draw("E");
    c_t_smallest_residual->SaveAs("c_t_smallest_residual.C");
    c_t_smallest_residual->SaveAs("c_t_smallest_residual.png");
    c_t_smallest_residual->SaveAs("c_t_smallest_residual.pdf");
        h_t_smallest_residual->Write();
    delete h_t_smallest_residual;
    delete c_t_smallest_residual;
    
    TCanvas *c_t_next_smallest_residual = new TCanvas("c_t_next_smallest_residual", "c_t_next_smallest_residual", 800, 600);
    h_t_next_smallest_residual->Draw("E");
    c_t_next_smallest_residual->SaveAs("c_t_next_smallest_residual.C");
    c_t_next_smallest_residual->SaveAs("c_t_next_smallest_residual.png");
    c_t_next_smallest_residual->SaveAs("c_t_next_smallest_residual.pdf");
        h_t_next_smallest_residual->Write();
    delete h_t_next_smallest_residual;
    delete c_t_next_smallest_residual;
    
    
    delete h_t_smallest;
    delete h_t_next_smallest;
    
    delete f_t0_smallest;
    delete f_t_smallest;
    delete f_t_next_smallest;
    
    std::cout << "\n>>> f_t_correlation" << std::endl;
    TCanvas *c_t_correlation = new TCanvas("c_t_correlation", "c_t_correlation", 800, 600);
    h_t_correlation->Fit("f_t_correlation");
    h_t_correlation->Draw("colz");
    c_t_correlation->SaveAs("c_t_correlation.C");
    c_t_correlation->SaveAs("c_t_correlation.png");
    c_t_correlation->SaveAs("c_t_correlation.pdf");
        h_t_correlation->Write();
    //delete h_t_correlation;
    delete c_t_correlation;
    
    // create residuals plot
    for(Int_t j = 1; j <= h_t_correlation->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_t_correlation->GetNbinsX(); ++ i)
        {
            Double_t func_x = h_t_correlation_residual->GetXaxis()->GetBinCenter(i);
            Double_t func_y = h_t_correlation_residual->GetYaxis()->GetBinCenter(j);
            Double_t func_eval = f_t_correlation->Eval(func_x, func_y);
            
            h_t_correlation_residual->SetBinContent(i, j, h_t_correlation->GetBinContent(i, j) - func_eval);
            h_t_correlation_residual->SetBinError(i, j, h_t_correlation->GetBinError(i, j));
        }
    }
    
    TCanvas *c_t_correlation_residual = new TCanvas("c_t_correlation_residual", "c_t_correlation_residual", 800, 600);
    h_t_correlation_residual->Draw("colz");
    c_t_correlation_residual->SaveAs("c_t_correlation_residual.C");
    c_t_correlation_residual->SaveAs("c_t_correlation_residual.png");
    c_t_correlation_residual->SaveAs("c_t_correlation_residual.pdf");
        h_t_correlation_residual->Write();
    delete h_t_correlation_residual;
    delete c_t_correlation_residual;
    
    delete h_t_correlation;
    
    
    /*
    std::cout << "chisq=" << f_t_correlation->GetChisquare() / f_t_correlation->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_correlation->GetParameter(0) << " +- " << f_t_correlation->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_correlation->GetParameter(1) << " +- " << f_t_correlation->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_correlation->GetParameter(2) << " +- " << f_t_correlation->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t_correlation->GetParameter(3) << " +- " << f_t_correlation->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t_correlation->GetParameter(4) << " +- " << f_t_correlation->GetParError(4) << std::endl;
    std::cout << "p5 = " << f_t_correlation->GetParameter(5) << " +- " << f_t_correlation->GetParError(5) << std::endl;

    Double_t A{f_t_correlation->GetParameter(0)};
    Double_t x0{f_t_correlation->GetParameter(1)};
    Double_t y0{f_t_correlation->GetParameter(2)};
    Double_t a{f_t_correlation->GetParameter(3)};
    Double_t b{f_t_correlation->GetParameter(4)};
    Double_t c{f_t_correlation->GetParameter(5)};
    */
    
    std::cout << "chisq=" << f_t_correlation->GetChisquare() / f_t_correlation->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_correlation->GetParameter(0) << " +- " << f_t_correlation->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_correlation->GetParameter(1) << " +- " << f_t_correlation->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_correlation->GetParameter(2) << " +- " << f_t_correlation->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t_correlation->GetParameter(3) << " +- " << f_t_correlation->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t_correlation->GetParameter(4) << " +- " << f_t_correlation->GetParError(4) << std::endl;
    std::cout << "p5 = " << f_t_correlation->GetParameter(5) << " +- " << f_t_correlation->GetParError(5) << std::endl;

    Double_t A{f_t_correlation->GetParameter(0)};
    Double_t x0{f_t_correlation->GetParameter(1)};
    Double_t y0{f_t_correlation->GetParameter(2)};
    Double_t sigma_x{f_t_correlation->GetParameter(3)};
    Double_t sigma_y{f_t_correlation->GetParameter(4)};
    Double_t theta{f_t_correlation->GetParameter(5)};
    
    delete f_t_correlation;
    
    // Do MC
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> d1(0.0, sigma_x);
    std::normal_distribution<double> d2(0.0, sigma_y); // rotate FIRST then translate to x0, y0
     
    for(Long64_t evt{0}; evt < 3957; ++ evt)
    {
        double u = d1(gen);
        double v = d2(gen);
        
        // inverse rotation - rotate back
        double x = x0 + u * std::cos(-theta) - v * std::sin(-theta);
        double y = y0 + u * std::sin(-theta) + v * std::cos(-theta);
        
        h_t_correlation_mc->Fill(x, y);
    }
    
    std::cout << "\n>>> f_t_correlation_mc" << std::endl;
    TCanvas *c_t_correlation_mc = new TCanvas("c_t_correlation_mc", "c_t_correlation_mc", 800, 600);
    h_t_correlation_mc->Fit("f_t_correlation_mc");
    h_t_correlation_mc->Draw("colz");
    c_t_correlation_mc->SaveAs("c_t_correlation_mc.C");
    c_t_correlation_mc->SaveAs("c_t_correlation_mc.png");
    c_t_correlation_mc->SaveAs("c_t_correlation_mc.pdf");
        h_t_correlation_mc->Write();
    delete h_t_correlation_mc;
    delete c_t_correlation_mc;
    
    std::cout << "chisq=" << f_t_correlation_mc->GetChisquare() / f_t_correlation_mc->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_correlation_mc->GetParameter(0) << " +- " << f_t_correlation_mc->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_correlation_mc->GetParameter(1) << " +- " << f_t_correlation_mc->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_correlation_mc->GetParameter(2) << " +- " << f_t_correlation_mc->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t_correlation_mc->GetParameter(3) << " +- " << f_t_correlation_mc->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t_correlation_mc->GetParameter(4) << " +- " << f_t_correlation_mc->GetParError(4) << std::endl;
    std::cout << "p5 = " << f_t_correlation_mc->GetParameter(5) << " +- " << f_t_correlation_mc->GetParError(5) << std::endl;
    
    delete f_t_correlation_mc;
    
    ////////////////////////////////////////////////////////////////////////////
    // CANVAS OUTPUT (METHOD 2)
    ////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n>>> f_t_pos" << std::endl;
    TCanvas *c_t_pos = new TCanvas("c_t_pos", "c_t_pos", 800, 600);
    h_t_pos->Fit("f_t_pos");
    h_t_pos->Draw("E");
    c_t_pos->SaveAs("c_t_pos.C");
    c_t_pos->SaveAs("c_t_pos.png");
    c_t_pos->SaveAs("c_t_pos.pdf");
        h_t_pos->Write();
    //delete h_t_pos;
    delete c_t_pos;
    
    std::cout << "\n>>> f_t_neg" << std::endl;
    TCanvas *c_t_neg = new TCanvas("c_t_neg", "c_t_neg", 800, 600);
    h_t_neg->Fit("f_t_neg");
    h_t_neg->Draw("E");
    c_t_neg->SaveAs("c_t_neg.C");
    c_t_neg->SaveAs("c_t_neg.png");
    c_t_neg->SaveAs("c_t_neg.pdf");
        h_t_neg->Write();
    //delete h_t_neg;
    delete c_t_neg;
    
    std::cout << "chisq = " << f_t_pos->GetChisquare() / f_t_pos->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_pos->GetParameter(0) << " +- " << f_t_pos->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_pos->GetParameter(1) << " +- " << f_t_pos->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_pos->GetParameter(2) << " +- " << f_t_pos->GetParError(2) << std::endl;
    
    std::cout << "chisq=" << f_t_neg->GetChisquare() / f_t_neg->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_neg->GetParameter(0) << " +- " << f_t_neg->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_neg->GetParameter(1) << " +- " << f_t_neg->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_neg->GetParameter(2) << " +- " << f_t_neg->GetParError(2) << std::endl;
    
    // create residuals plot
    for(Int_t i = 1; i <= h_t_pos->GetNbinsX(); ++ i)
    {
        Double_t func_x = h_t_pos_residual->GetXaxis()->GetBinCenter(i);
        Double_t func_eval = f_t_pos->Eval(func_x);
        
        h_t_pos_residual->SetBinContent(i, h_t_pos->GetBinContent(i) - func_eval);
        h_t_pos_residual->SetBinError(i, h_t_pos->GetBinError(i));
    }
    
    for(Int_t i = 1; i <= h_t_neg->GetNbinsX(); ++ i)
    {
        Double_t func_x = h_t_neg_residual->GetXaxis()->GetBinCenter(i);
        Double_t func_eval = f_t_neg->Eval(func_x);
        
        h_t_neg_residual->SetBinContent(i, h_t_neg->GetBinContent(i) - func_eval);
        h_t_neg_residual->SetBinError(i, h_t_neg->GetBinError(i));
    }
    
    
    TCanvas *c_t_pos_residual = new TCanvas("c_t_pos_residual", "c_t_pos_residual", 800, 600);
    h_t_pos_residual->Draw("E");
    c_t_pos_residual->SaveAs("c_t_pos_residual.C");
    c_t_pos_residual->SaveAs("c_t_pos_residual.png");
    c_t_pos_residual->SaveAs("c_t_pos_residual.pdf");
        h_t_pos_residual->Write();
    delete h_t_pos_residual;
    delete c_t_pos_residual;
    
    TCanvas *c_t_neg_residual = new TCanvas("c_t_neg_residual", "c_t_neg_residual", 800, 600);
    h_t_neg_residual->Draw("E");
    c_t_neg_residual->SaveAs("c_t_neg_residual.C");
    c_t_neg_residual->SaveAs("c_t_neg_residual.png");
    c_t_neg_residual->SaveAs("c_t_neg_residual.pdf");
        h_t_neg_residual->Write();
    delete h_t_neg_residual;
    delete c_t_neg_residual;
    
    
    delete h_t_pos;
    delete h_t_neg;
    
    delete f_t_pos;
    delete f_t_neg;
    
    std::cout << "\n>>> f_t_cor" << std::endl;
    TCanvas *c_t_cor = new TCanvas("c_t_cor", "c_t_cor", 800, 600);
    h_t_cor->Fit("f_t_cor");
    h_t_cor->Draw("colz");
    c_t_cor->SaveAs("c_t_cor.C");
    c_t_cor->SaveAs("c_t_cor.png");
    c_t_cor->SaveAs("c_t_cor.pdf");
        h_t_cor->Write();
    //delete h_t_cor;
    delete c_t_cor;
    
    // create residuals plot
    for(Int_t j = 1; j <= h_t_cor->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_t_cor->GetNbinsX(); ++ i)
        {
            Double_t func_x = h_t_cor_residual->GetXaxis()->GetBinCenter(i);
            Double_t func_y = h_t_cor_residual->GetYaxis()->GetBinCenter(j);
            Double_t func_eval = f_t_cor->Eval(func_x, func_y);
            
            h_t_cor_residual->SetBinContent(i, j, h_t_cor->GetBinContent(i, j) - func_eval);
            h_t_cor_residual->SetBinError(i, j, h_t_cor->GetBinError(i, j));
        }
    }
    
    TCanvas *c_t_cor_residual = new TCanvas("c_t_cor_residual", "c_t_cor_residual", 800, 600);
    h_t_cor_residual->Draw("colz");
    c_t_cor_residual->SaveAs("c_t_cor_residual.C");
    c_t_cor_residual->SaveAs("c_t_cor_residual.png");
    c_t_cor_residual->SaveAs("c_t_cor_residual.pdf");
        h_t_cor_residual->Write();
    delete h_t_cor_residual;
    delete c_t_cor_residual;
    
    delete h_t_cor;
    
    
    /*
    std::cout << "chisq=" << f_t_cor->GetChisquare() / f_t_cor->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_cor->GetParameter(0) << " +- " << f_t_cor->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_cor->GetParameter(1) << " +- " << f_t_cor->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_cor->GetParameter(2) << " +- " << f_t_cor->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t_cor->GetParameter(3) << " +- " << f_t_cor->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t_cor->GetParameter(4) << " +- " << f_t_cor->GetParError(4) << std::endl;
    std::cout << "p5 = " << f_t_cor->GetParameter(5) << " +- " << f_t_cor->GetParError(5) << std::endl;

    Double_t A{f_t_cor->GetParameter(0)};
    Double_t x0{f_t_cor->GetParameter(1)};
    Double_t y0{f_t_cor->GetParameter(2)};
    Double_t a{f_t_cor->GetParameter(3)};
    Double_t b{f_t_cor->GetParameter(4)};
    Double_t c{f_t_cor->GetParameter(5)};
    */
    
    std::cout << "chisq=" << f_t_cor->GetChisquare() / f_t_cor->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_cor->GetParameter(0) << " +- " << f_t_cor->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_cor->GetParameter(1) << " +- " << f_t_cor->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_cor->GetParameter(2) << " +- " << f_t_cor->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t_cor->GetParameter(3) << " +- " << f_t_cor->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t_cor->GetParameter(4) << " +- " << f_t_cor->GetParError(4) << std::endl;
    std::cout << "p5 = " << f_t_cor->GetParameter(5) << " +- " << f_t_cor->GetParError(5) << std::endl;

    Double_t A_cor{f_t_cor->GetParameter(0)};
    Double_t x0_cor{f_t_cor->GetParameter(1)};
    Double_t y0_cor{f_t_cor->GetParameter(2)};
    Double_t sigma_x_cor{f_t_cor->GetParameter(3)};
    Double_t sigma_y_cor{f_t_cor->GetParameter(4)};
    Double_t theta_cor{f_t_cor->GetParameter(5)};
    
    delete f_t_cor;
    
    // Do MC
    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::normal_distribution<double> d1_cor(0.0, sigma_x_cor);
    std::normal_distribution<double> d2_cor(0.0, sigma_y_cor); // rotate FIRST then translate to x0, y0
     
    for(Long64_t evt{0}; evt < 3957; ++ evt) // TODO: number
    {
        double u = d1_cor(gen);
        double v = d2_cor(gen);
        
        // inverse rotation - rotate back
        double x = x0_cor + u * std::cos(-theta_cor) - v * std::sin(-theta_cor);
        double y = y0_cor + u * std::sin(-theta_cor) + v * std::cos(-theta_cor);
        
        h_t_cor_mc->Fill(x, y);
    }
    
    std::cout << "\n>>> f_t_cor_mc" << std::endl;
    TCanvas *c_t_cor_mc = new TCanvas("c_t_cor_mc", "c_t_cor_mc", 800, 600);
    h_t_cor_mc->Fit("f_t_cor_mc");
    h_t_cor_mc->Draw("colz");
    c_t_cor_mc->SaveAs("c_t_cor_mc.C");
    c_t_cor_mc->SaveAs("c_t_cor_mc.png");
    c_t_cor_mc->SaveAs("c_t_cor_mc.pdf");
        h_t_cor_mc->Write();
    delete h_t_cor_mc;
    delete c_t_cor_mc;
    
    std::cout << "chisq=" << f_t_cor_mc->GetChisquare() / f_t_cor_mc->GetNDF() << std::endl;
    std::cout << "p0 = " << f_t_cor_mc->GetParameter(0) << " +- " << f_t_cor_mc->GetParError(0) << std::endl;
    std::cout << "p1 = " << f_t_cor_mc->GetParameter(1) << " +- " << f_t_cor_mc->GetParError(1) << std::endl;
    std::cout << "p2 = " << f_t_cor_mc->GetParameter(2) << " +- " << f_t_cor_mc->GetParError(2) << std::endl;
    std::cout << "p3 = " << f_t_cor_mc->GetParameter(3) << " +- " << f_t_cor_mc->GetParError(3) << std::endl;
    std::cout << "p4 = " << f_t_cor_mc->GetParameter(4) << " +- " << f_t_cor_mc->GetParError(4) << std::endl;
    std::cout << "p5 = " << f_t_cor_mc->GetParameter(5) << " +- " << f_t_cor_mc->GetParError(5) << std::endl;
    
    delete f_t_cor_mc;
    
    
    
    
    
    
    //delete t2;
    f2->Close();
    delete f2;
    
    f->Close();
    delete f;

    std::cout << "Stats from Preselection Cuts:" << std::endl;
    std::cout << "Number of accepted events: " << count_accept << "\nNumber of rejected events: " << count_reject << std::endl;
    std::cout << "Ratio accepted: " << (Double_t)count_accept / (Double_t)(count_accept + count_reject) << std::endl;

    std::cout << "Stats from Timestamp Selection:" << std::endl;
    std::cout << "Number of accepted events: " << count_accept_timestamp << "\nNumber of rejected events: " << count_reject_timestamp << std::endl;
    std::cout << "Ratio accepted: " << (Double_t)count_accept_timestamp / (Double_t)(count_accept_timestamp + count_reject_timestamp) << std::endl;

    return 0;
}
