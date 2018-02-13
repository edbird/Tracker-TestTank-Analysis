
/*
 *  Everything was essentially trivial. The only reason that I did
 *  not finish sooner was that I did not posess the required
 *  informations at the start.
 *
 */

#include "functions.hpp"
#include "datatypes.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <random>

//#include <RTypes.h>
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

    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

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

    TGaxis::SetMaxDigits(3);
    
    ////////////////////////////////////////////////////////////////////////////
    // DATA LOAD AND SAVE
    ////////////////////////////////////////////////////////////////////////////
    
    // Input file
    TFile *f = new TFile(filename.c_str()); //("cell8.root");
    f->cd();
    TTree *t = (TTree*)f->Get("histo");
    t->SetDirectory(f);

    // Local storage, testtank format (origional C++ analysis code - this code)
    TestTankStorage store;
    
    // Set default values
    TestTankStorage_init(&store);
    
    // SetBranchAddress for TTree t
    TestTankStorage_setbranchaddress(&store, t, falaise_mc_);
    
    // Output file
    std::string filename_out{filename.substr(0, filename.rfind(".")) + std::string("_out") + filename.substr(filename.rfind("."))};
    TFile *f2 = new TFile(filename_out.c_str(), "RECREATE"); //("cell8_out.root", "RECREATE");
    f2->cd();
    TTree *t2 = new TTree("histo", "");
    t2->SetDirectory(f2);

    // Branch for TTree t2
    TestTankStorage_branch(&store, t2, falaise_mc_);
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS GENERIC
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // PLASMA PROPAGATION TIME HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    TH1F *h_plasma_propagation_time = new TH1F("h_plasma_propagation_time", "h_plasma_propagation_time", 50, 0.0, 70.0); //35.0, 45.0);
    h_plasma_propagation_time->SetStats(0);

    TH2F *h_cathode_time_plasma_propagation_time = new TH2F("h_cathode_time_plasma_propagation_time", "h_cathode_time_plasma_propagation_time", 50, -20.0, 80.0, 50, 0.0, 70.0);
    h_cathode_time_plasma_propagation_time->SetStats(0);

    //TH2F *h_position_

    ////////////////////////////////////////////////////////////////////////////
    // FIT FUNCTION FOR PLASMA PROPAGATION TIME
    ////////////////////////////////////////////////////////////////////////////
    
    TF1* f_plasma_propagation_time = new TF1("f_plasma_propagation_time", ppt_fitf, 35.0, 65.0, 5);
    f_plasma_propagation_time->SetNpx(10000);
    f_plasma_propagation_time->SetParameter(0, 0.3); //4000.0);
    //f_plasma_propagation_time->SetParameter(1, 41.0);
    //f_plasma_propagation_time->SetParameter(2, 0.3);
    //f_plasma_propagation_time->SetParameter(3, 41.0);
    //f_plasma_propagation_time->FixParameter(3, 38.0);
    f_plasma_propagation_time->SetParameter(1, 39.0);
    f_plasma_propagation_time->SetParameter(2, 41.0);
    f_plasma_propagation_time->SetParameter(3, 0.3);
    f_plasma_propagation_time->SetParameter(4, 0.0);

    ////////////////////////////////////////////////////////////////////////////
    // CATHODE HISTOGRAM AND FIT FUNCTION
    ////////////////////////////////////////////////////////////////////////////
    
    TH1F *h_cathode_time = new TH1F("h_cathode_time", "h_cathode_time", 50, -20.0, 80.0); //-50, 50.0); // was 40.0 for data
    h_cathode_time->SetStats(0);
    
    TF1 *f_cathode_time = new TF1("f_cathode_time", cathode_time_fitf, -20.0, 80.0, 5);
    f_cathode_time->SetNpx(10000);
    f_cathode_time->SetParameter(0, 70.0e-3); // TODO: need to normalize
    f_cathode_time->SetParameter(1, 0.0);
    f_cathode_time->SetParLimits(1, 0.0, 1.0e1);
    f_cathode_time->SetParameter(2, 10.0);
    f_cathode_time->SetParameter(3, 30.0);
    f_cathode_time->SetParameter(4, 45.0);

    ////////////////////////////////////////////////////////////////////////////
    // FEAST TIMESTAMP HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////
    
    //TH1F *h_feast_t0 = new TH1F("h_feast_t0", "h_feast_t0", 50, -1.0, 7.0); //4.76, 4.86);
    TH1F *h_feast_t0 = new TH1F("h_feast_t0", "h_feast_t0", 50, 4.76, 4.86);
    TH1F *h_feast_t1 = new TH1F("h_feast_t1", "h_feast_t1", 50, 0.0, 80.0); // TODO: change ranges here
    TH1F *h_feast_t2 = new TH1F("h_feast_t2", "h_feast_t2", 50, 0.0, 80.0); // 0.0, 50.0
    TH1F *h_feast_t3 = new TH1F("h_feast_t3", "h_feast_t3", 50, 0.0, 80.0);
    TH1F *h_feast_t4 = new TH1F("h_feast_t4", "h_feast_t4", 50, 0.0, 80.0); 
    
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
    f_feast_t0->SetNpx(10000);
    f_feast_t0->SetParameter(0, 45.0e-3); // 250.0
    f_feast_t0->SetParameter(1, 4.780);
    f_feast_t0->SetParameter(2, 4.790);
    f_feast_t0->SetParameter(3, 4.829);
    f_feast_t0->SetParameter(4, 4.836);
    
    ////////////////////////////////////////////////////////////////////////////
    // FEAST T0 TO FEAST TX CORRELATIONS
    ////////////////////////////////////////////////////////////////////////////

    TH2F* h_feast_t0_t1_cor = new TH2F("h_feast_t0_t1_cor", "h_feast_t0_t1_cor", 50, 4.76, 4.86, 50, -10.0, 70.0); // TODO: change range here
    TH2F* h_feast_t0_t2_cor = new TH2F("h_feast_t0_t2_cor", "h_feast_t0_t2_cor", 50, 4.76, 4.86, 50, -10.0, 70.0);
    TH2F* h_feast_t0_t3_cor = new TH2F("h_feast_t0_t3_cor", "h_feast_t0_t3_cor", 50, 4.76, 4.86, 50, -10.0, 70.0);
    TH2F* h_feast_t0_t4_cor = new TH2F("h_feast_t0_t4_cor", "h_feast_t0_t4_cor", 50, 4.76, 4.86, 50, -10.0, 70.0);
    
    h_feast_t0_t1_cor->SetStats(0);
    h_feast_t0_t2_cor->SetStats(0);
    h_feast_t0_t3_cor->SetStats(0);
    h_feast_t0_t4_cor->SetStats(0);
    
    //TF2* f_feast_t0_t1_cor = new TF2("f_feast_t0_t1_cor", NULL, 0.0, 0.0);
    
    ////////////////////////////////////////////////////////////////////////////
    // FEAST T0 TO CATHODE CORRELATION AND FEAST T1 TO CATHODE CORRELATION
    ////////////////////////////////////////////////////////////////////////////
    
    TH2F* h_feast_t0_cathode_time_cor = new TH2F("h_feast_t0_cathode_time_cor", "h_feast_t0_cathode_time_cor", 50, 4.76, 4.86, 50, -10.0, 70.0); // TODO: change range here
    h_feast_t0_cathode_time_cor->SetStats(0);
    
    /*
    TH2F* h_feast_t1_cathode_time_cor = new TH2F("h_feast_t1_cathode_time_cor", "h_feast_t1_cathode_time_cor", 50, 0.0, 70.0, 50, -10.0, 70.0); // TODO: change range here
    h_feast_t1_cathode_time_cor->SetStats(0);
    */
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS AND FIT FUNCTIONS (METHOD 1)
    ////////////////////////////////////////////////////////////////////////////
    
    TH1F *h_t0_smallest = new TH1F("h_t0_smallest", "h_t0_smallest", 50, -5.0, 3.0); // TODO: changed limits below, do same here?
    TH1F *h_t_smallest = new TH1F("h_t_smallest", "h_t_smallest", 50, -0.4, 0.0); //4.3, 4.9);
    TH1F *h_t_next_smallest = new TH1F("h_t_next_smallest", "h_t_next_smallest", 50, 0.0, 0.4); //4.8, 5.2); // 2018-01-17 changed
    
    TH1F *h_t_smallest_residual = new TH1F("h_t_smallest_residual", "h_t_smallest_residual", 50, -0.4, 0.0); //4.3, 4.9);
    TH1F *h_t_next_smallest_residual = new TH1F("h_t_next_smallest_residual", "h_t_next_smallest_residual", 50, 0.0, 0.4); //4.8, 5.2);
   
    h_t0_smallest->SetStats(0);
    h_t_smallest->SetStats(0);
    h_t_next_smallest->SetStats(0);

    h_t_smallest_residual->SetStats(0);
    h_t_next_smallest_residual->SetStats(0);

    TF1 *f_t0_smallest = new TF1("f_t0_smallest", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", -5.0, 3.0);
    f_t0_smallest->SetNpx(10000); 
    f_t0_smallest->SetParameter(0, 20.0);
    f_t0_smallest->SetParameter(1, 0.0);
    f_t0_smallest->SetParameter(2, -1.0);
    f_t0_smallest->SetParameter(3, 0.0);
    f_t0_smallest->SetParameter(4, 0.0);
    
    TF1 *f_t_smallest = new TF1("f_t_smallest", "[0]*exp(-pow((x-[1])/[2], 2.0))", -0.4, 0.0); //4.3, 4.9);
    f_t_smallest->SetNpx(10000); 
    f_t_smallest->SetParameter(0, 1000.0); // 300.0
    f_t_smallest->SetParameter(1, -0.2);
    f_t_smallest->SetParameter(2, 0.05);
    
    TF1 *f_t_next_smallest = new TF1("f_t_next_smallest", "[0]*exp(-pow((x-[1])/[2], 2.0))", 0.0, 0.4); //4.8, 5.2);
    f_t_next_smallest->SetNpx(10000); 
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
    
    f_t_correlation->SetNpy(10000);
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
    
    f_t_correlation_mc->SetNpy(10000);
    f_t_correlation_mc->SetContour(10);
    
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS AND FIT FUNCTIONS (METHOD 2)
    ////////////////////////////////////////////////////////////////////////////
    
    TH1F *h_t_neg = new TH1F("h_t_neg", "h_t_neg", 50, -0.4, 0.0);
    TH1F *h_t_pos = new TH1F("h_t_pos", "h_t_pos", 50, 0.0, 0.4);
    
    TH1F *h_t_neg_residual = new TH1F("h_t_neg_residual", "h_t_neg_residual", 50, -0.4, 0.0);
    TH1F *h_t_pos_residual = new TH1F("h_t_pos_residual", "h_t_pos_residual", 50, 0.0, 0.4);
   
    h_t_neg->SetStats(0);
    h_t_pos->SetStats(0);
    
    h_t_neg_residual->SetStats(0);
    h_t_pos_residual->SetStats(0);
    
    TF1 *f_t_neg = new TF1("f_t_neg", "[0]*exp(-pow((x-[1])/[2], 2.0))", 4.3, 4.9);
    f_t_neg->SetNpx(10000); 
    f_t_neg->SetParameter(0, 300.0);
    f_t_neg->SetParameter(1, -0.2);
    f_t_neg->SetParameter(2, 0.05);
    
    TF1 *f_t_pos = new TF1("f_t_pos", "[0]*exp(-pow((x-[1])/[2], 2.0))", 4.8, 5.2);
    f_t_pos->SetNpx(10000); 
    f_t_pos->SetParameter(0, 350.0);
    f_t_pos->SetParameter(1, 0.2);
    f_t_pos->SetParameter(2, 0.05);
    
    TH2F *h_t_cor = new TH2F("h_t_cor", "h_t_cor", 50, -0.4, 0.0, 50, 0.0, 0.4);
    TH2F *h_t_cor_residual = new TH2F("h_t_cor_residual", "h_t_cor_residual", 50, -0.4, 0.0, 50, 0.0, 0.4);
    
    TH2F *h_t_cor_mc = new TH2F("h_t_cor_mc", "h_t_cor_mc", 50, -0.4, 0.0, 50, 0.0, 0.4);
    
    h_t_cor->SetStats(0);
    h_t_cor_residual->SetStats(0);
    h_t_cor_mc->SetStats(0);
    
    ////////////////////////////////////////////////////////////////////////////
    // Raw correlation histograms
    ////////////////////////////////////////////////////////////////////////////
    TH2F *h_t1_t3 = new TH2F("h_t1_t3", "h_t1_t3", 50, 0.0, 80.0, 50, 0.0, 80.0);
    TH2F *h_t2_t4 = new TH2F("h_t2_t4", "h_t2_t4", 50, 0.0, 80.0, 50, 0.0, 80.0);

    h_t1_t3 -> SetStats(0);
    h_t2_t4 -> SetStats(0);

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
    
    f_t_cor->SetNpy(10000);
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
    
    f_t_cor_mc->SetNpy(10000);
    f_t_cor_mc->SetContour(10);
    

    ////////////////////////////////////////////////////////////////////////////
    // Z POSITION CATHODE TIME CORRELATION HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    TH1F *h_position = new TH1F("h_position", "h_position", 50, -1.0, 1.0);
    TH1F *h_truth_position = new TH1F("h_truth_position", "h_truth_position", 50, -1.5, 1.5);
    TH1F *h_half_position = new TH1F("h_half_position", "h_half_position", 50, 0.0, 1.0);

    h_position->SetStats(0);
    h_truth_position->SetStats(0);
    h_half_position->SetStats(0);

    TH2F *h_zpos_cathode_time = new TH2F("h_zpos_cathode_time", "h_zpos_cathode_time", 20, -1.0, 1.0, 200, 0.0, 45.0);
    TH2F *h_zpos_cathode_time_residual = new TH2F("h_zpos_cathode_time_residual", "h_zpos_cathode_time_residual", 20, -1.0, 1.0, 200, 0.0, 45.0);
    
    h_zpos_cathode_time->SetStats(0);
    h_zpos_cathode_time_residual->SetStats(0);
    
    TF2 *f_zpos_cathode_time = new TF2("f_zpos_cathode_time", zpos_cathode_time_fitf, -1.0, 1.0, 0.0, 45.0, 4, 2);

    f_zpos_cathode_time->SetParameter(0, 285.0);
    f_zpos_cathode_time->SetParameter(1, 1.0); // 1.18
    f_zpos_cathode_time->SetParameter(2, 0.05);
    f_zpos_cathode_time->SetParameter(3, -0.051);

    f_zpos_cathode_time->SetParLimits(2, 0.0, 1000.0);
 
    f_zpos_cathode_time->SetNpy(10000);
    f_zpos_cathode_time->SetContour(10);

    ////////////////////////////////////////////////////////////////////////////
    // "2 dimensional shifted profile histogram"
    // 2 dimensional histogram where y values are shifted towards the mean
    // of the 2-d fit above
    ////////////////////////////////////////////////////////////////////////////
    
    Double_t actual_x_low = 0.0;
    Double_t actual_x_high = 1.0;
    Int_t actual_bins = 1;
    {
        // compute number of bins and limits
        Double_t x_low = 0.0; Double_t x_high = 45.0;
        Int_t n_bins = 200;
        Double_t bin_width = (x_high - x_low) / (Double_t) n_bins;
        Double_t target_range = (4.0) - (-4.0);
        Int_t target_bins = (Int_t)std::ceil(target_range / bin_width);
        Int_t actual_range = (Double_t)target_bins * bin_width;
        actual_x_low = 0.0 - 0.5 * actual_range;
        actual_x_high = 0.0 + 0.5 * actual_range;
        actual_bins = target_bins;
        std::cout << "actual_bins=" << actual_bins << std::endl;
    }
    TH2F *h_zpos_cathode_time_shift = new TH2F("h_zpos_cathode_time_shift", "h_zpos_cathode_time_shift", 20, -1.0, 1.0, actual_bins, actual_x_low, actual_x_high);
    
    h_zpos_cathode_time_shift->SetStats(0);

    //TH2F *h_zpos_cathode_time_shift_residual = new TH2F("h_zpos_cathode_time_shift_residual", "h_zpos_cathode_time_shift_residual", 20, -1.0, 1.0, actual_bins, actual_x_low, actual_x_high);

    //h_zpos_cathode_time_shift->SetStats(0);
    
    ////////////////////////////////////////////////////////////////////////////
    // "profile" histogram
    // 1 dimension collapsed (in y direction) histogram
    ////////////////////////////////////////////////////////////////////////////

    TH1F *h_zpos_cathode_time_profile = new TH1F("h_zpos_cathode_time_profile", "h_zpos_cathode_time_profile", 30, -3.0, 3.0); // TODO

    h_zpos_cathode_time_profile->SetStats(0);
    h_zpos_cathode_time_profile->Sumw2(kTRUE);

    TF1 *f_zpos_cathode_time_profile_1 = new TF1("f_zpos_cathode_time_profile_1", gaussian_fitf, -3.0, 3.0, 3);

    f_zpos_cathode_time_profile_1->SetParameter(0, 650.0);
    f_zpos_cathode_time_profile_1->SetParameter(1, 0.0);
    f_zpos_cathode_time_profile_1->SetParameter(2, 0.5);

    f_zpos_cathode_time_profile_1->SetParLimits(2, 0.0, 999.9);

    f_zpos_cathode_time_profile_1->SetNpx(1000);
    
    TF1 *f_zpos_cathode_time_profile_2 = new TF1("f_zpos_cathode_time_profile_2", gaussian_fitf, -3.0, 3.0, 3);

    f_zpos_cathode_time_profile_2->SetParameter(0, 100.0);
    f_zpos_cathode_time_profile_2->SetParameter(1, 0.0);
    f_zpos_cathode_time_profile_2->SetParameter(2, 0.5);

    f_zpos_cathode_time_profile_2->SetParLimits(2, 0.0, 999.9);

    f_zpos_cathode_time_profile_2->SetNpx(1000);

    TF1 *f_zpos_cathode_time_profile = new TF1("f_zpos_cathode_time_profile", double_gaussian_fitf, -3.0, 3.0, 6);

    f_zpos_cathode_time_profile->SetParameter(0, 400.0);
    f_zpos_cathode_time_profile->SetParameter(1, -1.0);
    f_zpos_cathode_time_profile->SetParameter(2, 0.5);
    f_zpos_cathode_time_profile->FixParameter(3, 0.0);
    f_zpos_cathode_time_profile->SetParameter(4, 1.0);
    f_zpos_cathode_time_profile->SetParameter(5, 0.5);
    
    f_zpos_cathode_time_profile->SetParLimits(2, 0.0, 999.9);
    f_zpos_cathode_time_profile->SetParLimits(5, 0.0, 999.9);
    
    f_zpos_cathode_time_profile->SetNpx(1000);


    ////////////////////////////////////////////////////////////////////////////
    // DATA LOOP
    ////////////////////////////////////////////////////////////////////////////
   
    #define COUT_TIMESTAMP_GOOD 0 // 0
    #define COUT_TIMESTAMP_FAIL 0 // 1
    #define COUT_TIMESTAMP_WAIT 0 // 0
    #define WAVEFORM_PRINT_GOOD 1
    #define WAVEFORM_PRINT_FAIL 1

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
        
        h_cathode_time_plasma_propagation_time->Fill(store.cathode_time, store.plasma_propagation_time);

        // Can only fill for MC, because otherwise store.truth_position contains
        // zero
        // Note: Fill anyway, but histogram contains only zero
        if(falaise_mc_ == true)
        {
            h_truth_position->Fill(store.truth_position);
        }
        else
        {
            h_truth_position->Fill(0.0);
        }

        ////////////////////////////////////////////////////////////////////////
        // NOTE TO SELF: New varaible store.truth_position introduced. This
        // variable is only valid for MC (falaise simulation) data.
        ////////////////////////////////////////////////////////////////////////

        // check anode and cathode peaks - only for real data
        bool event_check_flag{true};
        if(falaise_mc_ == false)
        {
            event_check_flag = cut_l(store.anode_peak, 250.0) && cut_l(store.cathode_peak, 350.0);
        }


        bool event_check_flag_position = false;
        if(falaise_mc_ == false)
        {
            if(
                cut(store.plasma_propagation_time, 38.5, 42.0) &&
                cut(store.position, -0.95, 0.95)
                )
            {
                event_check_flag_position = true;
            }
        }
        else
        {
            if(
                cut(store.plasma_propagation_time, 0.1, 9999.9) && // TODO: these are temp. test values
                cut(store.position, -0.95, 0.95)
                )
            {
                event_check_flag_position = true;
            }
        }
        if(event_check_flag_position == true)
        {

            ////////////////////////////////////////////////////////////////////
            // Z POSITION CATHODE TIME CORRELATION HISTOGRAMS
            ////////////////////////////////////////////////////////////////////

            h_position->Fill(store.position);
            h_half_position->Fill(store.half_position);

            h_zpos_cathode_time->Fill(store.position, store.cathode_time);
        
            // Note: Parameters MUST be entered correctly here 
            Double_t zpos_cathode_time_mean = 1.00427; //f_zpos_cathode_time->GetParameter(1);
            Double_t zpos_cathode_time_theta = -0.0506056; //f_zpos_cathode_time->GetParameter(3);
            Double_t func_x = store.position;
            Double_t func_y = store.cathode_time;
            Double_t profile_func_eval = zpos_cathode_time_profilef(func_x, zpos_cathode_time_mean, zpos_cathode_time_theta);
            
            // create shifted histogram
            h_zpos_cathode_time_shift->Fill(func_x, func_y - profile_func_eval);
             
            // create collapsed profile
            h_zpos_cathode_time_profile->Fill(func_y - profile_func_eval);

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
                cut(store.plasma_propagation_time, 0.1, 9999.9) && //57.0, 59.0) && // TODO: these are temp. test values
                cut(store.position, -0.95, 0.95) &&
                event_check_flag
            )
        )
        {
            ++ count_accept;

            // Sum anode t0 timestamp times for computation of mean time
            mean_anode_time += store.feast_t0;
            
            //std::cout << "Good event" << std::endl;
            //std::cout << "Event index is: " << ix << std::endl;

            // TODO: shift all the timestamp variables so that they make more logical sense
            // this is currently done in the falaise module, however should be done here
            // instead

            #if COUT_TIMESTAMP_GOOD
                timestamp_print(std::cout, store);
                std::cout << std::endl;
                
                #if COUT_TIMESTAMP_WAIT
                    std::cin.get();
                #endif
            #endif
            
            //std::cout << "fill " << ix << std::endl;
            t2->Fill();
            


            ////////////////////////////////////////////////////////////////////
            // FEAST TIME DIFFERENCES
            ////////////////////////////////////////////////////////////////////

            //h_cathode_time->Fill(store.cathode_time);
            
            // TODO: moved to below
            //h_feast_t0->Fill(store.feast_t0);
            //h_feast_t1->Fill(store.feast_t1);
            //h_feast_t2->Fill(store.feast_t2);
            //h_feast_t3->Fill(store.feast_t3);
            //h_feast_t4->Fill(store.feast_t4);
            
            //h_cathode->Fill(store.cathode_time);
            
            h_feast_t0_diff->Fill(store.feast_t0 - store.cathode_time);
            h_feast_t1_diff->Fill(store.feast_t1 - store.feast_t0 - store.cathode_time);
            h_feast_t2_diff->Fill(store.feast_t2 - store.feast_t0 - store.cathode_time);
            h_feast_t3_diff->Fill(store.feast_t3 - store.feast_t0 - store.cathode_time);
            h_feast_t4_diff->Fill(store.feast_t4 - store.feast_t0 - store.cathode_time);
            
            ////////////////////////////////////////////////////////////////////
            // Correlations
            ////////////////////////////////////////////////////////////////////
            
            h_feast_t0_t1_cor->Fill(store.feast_t0, store.feast_t1);
            h_feast_t0_t2_cor->Fill(store.feast_t0, store.feast_t2);
            h_feast_t0_t3_cor->Fill(store.feast_t0, store.feast_t3);
            h_feast_t0_t4_cor->Fill(store.feast_t0, store.feast_t4);
            
            h_feast_t0_cathode_time_cor->Fill(store.feast_t0, store.cathode_time);
            //h_feast_t1_cathode_time_cor->Fill(store.feast_t1, store.cathode_time);
            
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
            // TODO: could look at correlation of feast_t0 to cathode time?
            // Or feast_t0 to feast_t1
            // Note: cathode_time already has feast_t0 subtracted
            // Hence code below subtracts feast_t0 from feast_tX before subtracting
            // cathode time
            Double_t t1_diff = store.feast_t1 - store.feast_t0 - store.cathode_time;
            Double_t t2_diff = store.feast_t2 - store.feast_t0 - store.cathode_time;
            Double_t t3_diff = store.feast_t3 - store.feast_t0 - store.cathode_time;
            Double_t t4_diff = store.feast_t4 - store.feast_t0 - store.cathode_time; // TODO: maybe these are missleading
            
            // Absolute value differences
            Double_t t0_abs_diff = std::abs(store.feast_t0 - store.cathode_time); // TODO: maybe these are not missleading?
            Double_t t1_abs_diff = std::abs(store.feast_t1 /*- store.feast_t0*/ - store.cathode_time); // TODO: I added " - store.feast_t0" on 2018-01-17
            Double_t t2_abs_diff = std::abs(store.feast_t2 /*- store.feast_t0*/ - store.cathode_time);
            Double_t t3_abs_diff = std::abs(store.feast_t3 /*- store.feast_t0*/ - store.cathode_time);
            Double_t t4_abs_diff = std::abs(store.feast_t4 /*- store.feast_t0*/ - store.cathode_time);
            
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
            
            if(t1_diff <= 0.0) // note: added <= instead of <
                sort_me_neg.push_back(t1_diff);
            else if(t1_diff > 0.0)
                sort_me_pos.push_back(t1_diff);
            
            if(t2_diff <= 0.0)
                sort_me_neg.push_back(t2_diff);
            else if(t2_diff > 0.0)
                sort_me_pos.push_back(t2_diff);
            
            if(t3_diff <= 0.0)
                sort_me_neg.push_back(t3_diff);
            else if(t3_diff > 0.0)
                sort_me_pos.push_back(t3_diff);
            
            if(t4_diff <= 0.0)
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
                std::string output_file_name_differential("differential_");
                //canvas_name_counter = ix;
                Long64_t ix_copy{ix};
                #if WAVEFORM_PRINT_GOOD
                    store.anode_histo->SetStats(0);
                    waveform_print(store.anode_histo, ix_copy /*canvas_name_counter*/, output_file_name, "anode_histo_good");

                    Double_t xlow{store.anode_histo->GetBinLowEdge(1)};
                    Double_t xhigh{store.anode_histo->GetBinLowEdge(store.anode_histo->GetNbinsX()) + store.anode_histo->GetBinWidth(store.anode_histo->GetNbinsX())};
                    Int_t nbinsx{store.anode_histo->GetNbinsX()};
                    TH1F *differential_histo = new TH1F("differential_histo", "differential_histo", nbinsx, xlow, xhigh);
                    histogram_differentiate(differential_histo, store.anode_histo, differentiate_method_simple);
                    store.differential_histo = differential_histo;
                    waveform_print(differential_histo, ix_copy, output_file_name_differential, "differential_histo_good");
                #endif
                
                
                // TODO: Note to self, moved fill of timestamps to here to
                // see if eranious events dissapear
                h_cathode_time->Fill(store.cathode_time);
                
                //if((store.feast_t0 == *sort_me_neg.rbegin()) || (store.feast_t0 == *sort_me_pos.begin()))
                //{
                    h_feast_t0->Fill(store.feast_t0);
                //}
                //if((t1_diff == *sort_me_neg.rbegin()) || (t1_diff == *sort_me_pos.begin()))
                //{
                    h_feast_t1->Fill(store.feast_t1);
                //}
                //if((t2_diff == *sort_me_neg.rbegin()) || (t2_diff == *sort_me_pos.begin()))
                //{
                    h_feast_t2->Fill(store.feast_t2);
                //}
                //if((t3_diff == *sort_me_neg.rbegin()) || (t3_diff == *sort_me_pos.begin()))
                //{
                    h_feast_t3->Fill(store.feast_t3);
                //}
                //if((t4_diff == *sort_me_neg.rbegin()) || (t4_diff == *sort_me_pos.begin()))
                //{
                    h_feast_t4->Fill(store.feast_t4);
                //}
           
                h_t1_t3->Fill(store.feast_t1, store.feast_t3);
                h_t2_t4->Fill(store.feast_t2, store.feast_t4);

                //h_cathode->Fill(store.cathode_time);
            
                // Note: It does NOT remove the eranious bumps which are seen in
                // the t3 and t1 histograms
            }
            else
            {
                ++ count_reject_timestamp; 

                std::string output_file_name("anode_");
                std::string output_file_name_differential("differential_");
                //canvas_name_counter = ix;
                Long64_t ix_copy{ix};
                //waveform_print(ix_copy /*canvas_name_counter*/, cathode_histo, output_file_name);
                #if WAVEFORM_PRINT_FAIL
                    store.anode_histo->SetStats(0);
                    waveform_print(store.anode_histo, ix_copy /*canvas_name_counter*/, output_file_name, "anode_histo_fail");

                    Double_t xlow{store.anode_histo->GetBinLowEdge(1)};
                    Double_t xhigh{store.anode_histo->GetBinLowEdge(store.anode_histo->GetNbinsX()) + store.anode_histo->GetBinWidth(store.anode_histo->GetNbinsX())};
                    Int_t nbinsx{store.anode_histo->GetNbinsX()};
                    TH1F *differential_histo = new TH1F("differential_histo", "differential_histo", nbinsx, xlow, xhigh);
                    histogram_differentiate(differential_histo, store.anode_histo, differentiate_method_simple);
                    store.differential_histo = differential_histo;
                    waveform_print(differential_histo, ix_copy, output_file_name_differential, "differential_histo_fail");
                #endif
        
                // NOTE: This isn't the fail part of the first if statement
                std::cout << "Rejecting event with no timestamp either side" << std::endl;
                std::cout << "Event index is: " << ix << std::endl;
                std::cout << "Printing event to file: " << output_file_name << std::endl;
                std::cout << "sort_me_neg.size() = " << sort_me_neg.size() << std::endl;
                std::cout << "sort_me_pos.size() = " << sort_me_pos.size() << std::endl;
                
                #if COUT_TIMESTAMP_FAIL
                    timestamp_print(std::cout, store);
                    std::cout << std::endl;

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
            // TODO: can't remember why I chose 0.065 here, so removing
            // Note: This is in case there is a close 3rd timestamp which
            // occurs close to the second timestamp - or something like this
            /*
            else if(std::abs(sort_me.at(2) - small_abs_diff_2) <= 3.0 * 0.065)
            {
                std::cout << "Rejecting event with close times" << std::endl;
            }
            */
            else
            {
                // t0 may still be smaller than small_abs_diff_2
                // however assume it is not
                h_t_smallest->Fill(small_abs_diff_1 - store.feast_t0); // TODO: no idea if this is supposed to be subtracted or not anymore
                h_t_next_smallest->Fill(small_abs_diff_2 - store.feast_t0);
                
                //h_t_correlation->Fill(small_abs_diff_1 - store.feast_t0, small_abs_diff_2 - store.feast_t0);
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
    
    //Double_t integral;


    ////////////////////////////////////////////////////////////////////////////
    // CATHODE TIME PLASMA PROPAGATION TIME
    ////////////////////////////////////////////////////////////////////////////

    canvas(h_cathode_time_plasma_propagation_time, "c_cathode_time_plasma_propagation_time");
    delete h_cathode_time_plasma_propagation_time;


    ////////////////////////////////////////////////////////////////////////////
    // BASIC RAW TIMESTAMP HISTOGRAMS AND PLASMA PROPAGATION TIME HISTOGRAM
    ////////////////////////////////////////////////////////////////////////////

    #define TIMESTAMP_CANVAS_ENABLE 1 // TODO:move
    #if TIMESTAMP_CANVAS_ENABLE
        std::cout << "\n>>> f_plasma_propagation_time" << std::endl;
        canvas_scale_fit(h_plasma_propagation_time, "c_plasma_propagation_time", "f_plasma_propagation_time");
        delete h_plasma_propagation_time;
        
        fit_param_print(std::cout, f_plasma_propagation_time);
        std::cout << std::endl;

        std::cout << "\n>>> f_cathode_time" << std::endl;
        canvas_scale_fit(h_cathode_time, "c_cathode_time", "f_cathode_time");
        delete h_cathode_time;
        
        fit_param_print(std::cout, f_cathode_time);
        std::cout << std::endl;

        std::cout << "\n>>> f_feast_t0" << std::endl;
        canvas_scale_fit(h_feast_t0, "c_feast_t0", "f_feast_t0");
        delete h_feast_t0;
       
        // Print parameters for feast_t0 fit
        fit_param_print(std::cout, f_feast_t0);
        std::cout << std::endl;

        canvas(h_feast_t1, "c_feast_t1");
        delete h_feast_t1;
        
        canvas(h_feast_t2, "c_feast_t2");
        delete h_feast_t2;
        
        canvas(h_feast_t3, "c_feast_t3");
        delete h_feast_t3;
            
        canvas(h_feast_t4, "c_feast_t4");
        delete h_feast_t4;
        
        canvas(h_t1_t3, "c_t1_t3");
        delete h_t1_t3;

        canvas(h_t2_t4, "c_t2_t4");
        delete h_t2_t4;
        
        /*TCanvas *c_cathode = new TCanvas("c_cathode", "c_cathode", 800, 600);
        h_cathode->Draw();
        c_cathode->SaveAs("c_cathode.C");
        c_cathode->SaveAs("c_cathode.png");
        c_cathode->SaveAs("c_cathode.pdf");
        h_cathode->Write();
        delete h_cathode;
        delete c_cathode;*/
    
        // Print parameters for feast_t0 fit
        //fit_param_print(std::cout, f_feast_t0, 5);
        //std::cout << std::endl;

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
    // FEAST T0 T1 CORRELATION
    ////////////////////////////////////////////////////////////////////////////

    canvas(h_feast_t0_t1_cor, "c_feast_t0_t1_cor");
    delete h_feast_t0_t1_cor;
    
    canvas(h_feast_t0_t2_cor, "c_feast_t0_t2_cor");
    delete h_feast_t0_t2_cor;
    
    canvas(h_feast_t0_t3_cor, "c_feast_t0_t3_cor");
    delete h_feast_t0_t3_cor;
    
    canvas(h_feast_t0_t4_cor, "c_feast_t0_t4_cor");
    delete h_feast_t0_t4_cor;
    
    ////////////////////////////////////////////////////////////////////////////
    // FEAST T0 TO CATHODE CORRELATION AND FEAST T1 TO CATHODE CORRELATION
    ////////////////////////////////////////////////////////////////////////////
    
    canvas(h_feast_t0_cathode_time_cor, "c_feast_t0_cathode_time_cor");
    delete h_feast_t0_cathode_time_cor;
    
    /*
    TCanvas *c_feast_t1_cathode_time_cor = new TCanvas("c_feast_t1_cathode_time_cor", "c_feast_t1_cathode_time_cor", 800, 600);
    h_feast_t1_cathode_time_cor->Draw("colz");
    c_feast_t1_cathode_time_cor->SaveAs("c_feast_t1_cathode_time_cor.C");
    c_feast_t1_cathode_time_cor->SaveAs("c_feast_t1_cathode_time_cor.png");
    c_feast_t1_cathode_time_cor->SaveAs("c_feast_t1_cathode_time_cor.pdf");
    h_feast_t1_cathode_time_cor->Write();
    delete h_feast_t1_cathode_time_cor;
    delete c_feast_t1_cathode_time_cor;
    */
    
    ////////////////////////////////////////////////////////////////////////////
    // CANVAS OUTPUT (METHOD 1)
    ////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n>>> f_t0_smallest" << std::endl;
    canvas_fit(h_t0_smallest, "c_t0_smallest", "f_t0_smallest");
    delete h_t0_smallest;
    
    fit_param_print(std::cout, f_t0_smallest);
    std::cout << std::endl;
    
    std::cout << "\n>>> f_t_smallest" << std::endl;
    canvas_fit(h_t_smallest, "c_t_smallest", "f_t_smallest");
    //delete h_t_smallest;
    
    fit_param_print(std::cout, f_t_smallest);
    std::cout << std::endl;
    
    std::cout << "\n>>> f_t_next_smallest" << std::endl;
    canvas_fit(h_t_next_smallest, "c_t_next_smallest", "f_t_next_smallest");
    //delete h_t_next_smallest;
        
    fit_param_print(std::cout, f_t_next_smallest);
    std::cout << std::endl;
    
    /*
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
    */
    
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
    
    
    canvas(h_t_smallest_residual, "c_t_smallest_residual");
    delete h_t_smallest_residual;
    
    canvas(h_t_next_smallest_residual, "c_t_next_smallest_residual");
    delete h_t_next_smallest_residual;
    
    
    delete h_t_smallest;
    delete h_t_next_smallest;
    
    delete f_t0_smallest;
    delete f_t_smallest;
    delete f_t_next_smallest;
    
    std::cout << "\n>>> f_t_correlation" << std::endl;
    canvas_fit(h_t_correlation, "c_t_correlation", "f_t_correlation");
    //delete h_t_correlation;
    
    // create residuals plot
    for(Int_t j = 1; j <= h_t_correlation->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_t_correlation->GetNbinsX(); ++ i)
        {
            if(h_t_correlation->GetBinError(i, j) > 0.0) // check "empty"
            {
                Double_t func_x = h_t_correlation_residual->GetXaxis()->GetBinCenter(i);
                Double_t func_y = h_t_correlation_residual->GetYaxis()->GetBinCenter(j);
                Double_t func_eval = f_t_correlation->Eval(func_x, func_y);
                
                h_t_correlation_residual->SetBinContent(i, j, h_t_correlation->GetBinContent(i, j) - func_eval);
                h_t_correlation_residual->SetBinError(i, j, h_t_correlation->GetBinError(i, j));
            }
        }
    }
    
    canvas(h_t_correlation_residual, "c_t_correlation_residual");
    delete h_t_correlation_residual;
    
    delete h_t_correlation;
    
    
    /*
    Double_t A{f_t_correlation->GetParameter(0)};
    Double_t x0{f_t_correlation->GetParameter(1)};
    Double_t y0{f_t_correlation->GetParameter(2)};
    Double_t a{f_t_correlation->GetParameter(3)};
    Double_t b{f_t_correlation->GetParameter(4)};
    Double_t c{f_t_correlation->GetParameter(5)};
    */
    
    fit_param_print(std::cout, f_t_correlation);
    std::cout << std::endl;
    
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
    canvas_fit(h_t_correlation_mc, "c_t_correlation_mc", "f_t_correlation_mc");
    delete h_t_correlation_mc;
    
    fit_param_print(std::cout, f_t_correlation_mc);
    std::cout << std::endl;
    
    delete f_t_correlation_mc;
    
    ////////////////////////////////////////////////////////////////////////////
    // CANVAS OUTPUT (METHOD 2)
    ////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n>>> f_t_pos" << std::endl;
    canvas_fit(h_t_pos, "c_t_pos", "f_t_pos");
    //delete h_t_pos;
    
    fit_param_print(std::cout, f_t_pos);
    std::cout << std::endl;
    
    std::cout << "\n>>> f_t_neg" << std::endl;
    canvas_fit(h_t_neg, "c_t_neg", "f_t_neg");
    //delete h_t_neg;
    
    
    fit_param_print(std::cout, f_t_neg);
    std::cout << std::endl;
    
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
    
    
    canvas(h_t_pos_residual, "c_t_pos_residual");
    delete h_t_pos_residual;

    canvas(h_t_neg_residual, "c_t_neg_residual");
    delete h_t_neg_residual;
    
    
    delete h_t_pos;
    delete h_t_neg;
    
    delete f_t_pos;
    delete f_t_neg;
    
    std::cout << "\n>>> f_t_cor" << std::endl;
    canvas_fit(h_t_cor, "c_t_cor", "f_t_cor");
    //delete h_t_cor;
    
    // create residuals plot
    for(Int_t j = 1; j <= h_t_cor->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_t_cor->GetNbinsX(); ++ i)
        {
            if(h_t_cor->GetBinError(i, j) > 0.0) // check "empty"
            {
                Double_t func_x = h_t_cor_residual->GetXaxis()->GetBinCenter(i);
                Double_t func_y = h_t_cor_residual->GetYaxis()->GetBinCenter(j);
                Double_t func_eval = f_t_cor->Eval(func_x, func_y);
                
                h_t_cor_residual->SetBinContent(i, j, h_t_cor->GetBinContent(i, j) - func_eval);
                h_t_cor_residual->SetBinError(i, j, h_t_cor->GetBinError(i, j));
            }
        }
    }
   
    canvas(h_t_cor_residual, "c_t_cor_residual");
    delete h_t_cor_residual;
    
    delete h_t_cor;
    
    
    /*
    Double_t A{f_t_cor->GetParameter(0)};
    Double_t x0{f_t_cor->GetParameter(1)};
    Double_t y0{f_t_cor->GetParameter(2)};
    Double_t a{f_t_cor->GetParameter(3)};
    Double_t b{f_t_cor->GetParameter(4)};
    Double_t c{f_t_cor->GetParameter(5)};
    */
    
    fit_param_print(std::cout, f_t_cor);
    std::cout << std::endl;
    
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
    canvas_fit(h_t_cor_mc, "c_t_cor_mc", "f_t_cor_mc");
    delete h_t_cor_mc;
    
    fit_param_print(std::cout, f_t_cor_mc);
    std::cout << std::endl;
    delete f_t_cor_mc;
    
    

    ////////////////////////////////////////////////////////////////////////////
    // Z POSITION CATHODE TIME CORRELATION HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    canvas(h_position, "c_position");
    delete h_position;

    canvas(h_truth_position, "c_truth_position");
    delete h_truth_position;

    canvas(h_half_position , "c_half_position");
    delete h_half_position;

    std::cout << "\n>>> f_zpos_cathode_time" << std::endl;
    canvas_fit(h_zpos_cathode_time, "c_zpos_cathode_time", "f_zpos_cathode_time");
    //delete h_zpos_cathode_time;

    fit_param_print(std::cout, f_zpos_cathode_time);
    std::cout << std::endl;

    // create residuals plot
    for(Int_t j = 1; j <= h_zpos_cathode_time->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_zpos_cathode_time->GetNbinsX(); ++ i)
        {
            if(h_zpos_cathode_time->GetBinError(i, j) > 0.0) // check "empty"
            {
                //std::cout << "error > 0.0" << std::endl;

                Double_t func_x = h_zpos_cathode_time_residual->GetXaxis()->GetBinCenter(i);
                Double_t func_y = h_zpos_cathode_time_residual->GetYaxis()->GetBinCenter(j);
                Double_t func_eval = f_zpos_cathode_time->Eval(func_x, func_y);
               
                // TODO: this isn't working with fill for some reason
                //std::cout << h_zpos_cathode_time->GetBinContent(i, j) - func_eval << std::endl;

                h_zpos_cathode_time_residual->SetBinContent(i, j, h_zpos_cathode_time->GetBinContent(i, j) - func_eval);
                //h_zpos_cathode_time_residual->Fill(i, j, h_zpos_cathode_time->GetBinContent(i, j) - func_eval);
                //h_zpos_cathode_time_residual->SetBinError(i, j, h_zpos_cathode_time->GetBinError(i, j));
            }
        }
    }

    // create shifted histogram
    // Note: Moved to new data loop
    // Note: Not actually moved to new loop - to save processing
    // time parameters are manually entered into above data loop
    /*
    Double_t zpos_cathode_time_mean = f_zpos_cathode_time->GetParameter(1);
    Double_t zpos_cathode_time_theta = f_zpos_cathode_time->GetParameter(3);
    for(Int_t i = 1; i <= h_zpos_cathode_time->GetNbinsX(); ++ i)
    {
        Double_t func_x = h_zpos_cathode_time->GetXaxis()->GetBinCenter(i);
        Double_t profile_func_eval = zpos_cathode_time_profilef(func_x, zpos_cathode_time_mean, zpos_cathode_time_theta);

        for(Int_t j = 1; j <= h_zpos_cathode_time->GetNbinsY(); ++ j)
        {
            Double_t func_y = h_zpos_cathode_time->GetYaxis()->GetBinCenter(j);
            Double_t content = h_zpos_cathode_time->GetBinContent(i, j);

            if(content > 0.0)
            {
                h_zpos_cathode_time_shift->Fill(func_x, func_y - profile_func_eval, content);
            }

        }
    }
    */
    

    // create collapsed profile
    // Note: Moved to new data loop
    // Note: Not actually moved to new loop - to save processing
    // time parameters are manually entered into above data loop
    //Double_t zpos_cathode_time_mean = f_zpos_cathode_time->GetParameter(1);
    //Double_t zpos_cathode_time_theta = f_zpos_cathode_time->GetParameter(3);
    /*
    for(Int_t i = 1; i <= h_zpos_cathode_time->GetNbinsX(); ++ i)
    {
        Double_t func_x = h_zpos_cathode_time->GetXaxis()->GetBinCenter(i);
        Double_t profile_func_eval = zpos_cathode_time_profilef(func_x, zpos_cathode_time_mean, zpos_cathode_time_theta);

        for(Int_t j = 1; j <= h_zpos_cathode_time->GetNbinsY(); ++ j)
        {
            Double_t func_y = h_zpos_cathode_time->GetYaxis()->GetBinCenter(j);
            Double_t content = h_zpos_cathode_time->GetBinContent(i, j);

            if(content > 0.0)
            {
                h_zpos_cathode_time_profile->Fill(func_y - profile_func_eval, content);
            }

        }
    }
    */

    delete h_zpos_cathode_time;
    delete f_zpos_cathode_time;
    
    canvas(h_zpos_cathode_time_shift, "c_zpos_cathode_time_shift");
    delete h_zpos_cathode_time_shift;

    canvas(h_zpos_cathode_time_residual, "c_z_pos_cathode_time_residual"); 
    delete h_zpos_cathode_time_residual;
   
    std::cout << "\n>>> f_zpos_cathode_time_profile" << std::endl;
    canvas_fit(h_zpos_cathode_time_profile, "c_zpos_cathode_time_profile", "f_zpos_cathode_time_profile");
    delete h_zpos_cathode_time_profile;

    fit_param_print(std::cout, f_zpos_cathode_time_profile);
    std::cout << std::endl;
    delete f_zpos_cathode_time_profile;


    f2->Close();
    delete f2;
    
    f->Close();
    delete f;

    std::cout << "\n";
    std::cout << "Stats from Preselection Cuts:" << std::endl;
    std::cout << "Number of accepted events: " << count_accept << "\nNumber of rejected events: " << count_reject << std::endl;
    std::cout << "Ratio accepted: " << (Double_t)count_accept / (Double_t)(count_accept + count_reject) << std::endl;

    std::cout << "Stats from Timestamp Selection:" << std::endl;
    std::cout << "Number of accepted events: " << count_accept_timestamp << "\nNumber of rejected events: " << count_reject_timestamp << std::endl;
    std::cout << "Ratio accepted: " << (Double_t)count_accept_timestamp / (Double_t)(count_accept_timestamp + count_reject_timestamp) << std::endl;

    return 0;
}
