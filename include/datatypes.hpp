#ifndef DATATYPES_HPP
#define DATATYPES_HPP

//#include <ROOT/RTypes.h>
#include <TH1F.h>
#include <TTree.h>

////////////////////////////////////////////////////////////////////////////////
// TEST TANK DATA STRUCTURE
////////////////////////////////////////////////////////////////////////////////

// TODO: Move this to a header file somewhere else so that both TestModule
// and TrackerTestTankAnalysis can see it (only 1 copy of code)
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

    // added 2018-02-13
    TH1F *anode_smooth_histo = (TH1F*)0;
    TH1F *anode_average_histo = (TH1F*)0;
    TH1F *anode_differential_histo = (TH1F*)0;
    TH1F *anode_average_differential_histo = (TH1F*)0;

    Float_t truth_position = 0;
};


// Helper functions
void TestTankStorage_init(TestTankStorage * const testtankstorage_);
void TestTankStorage_setbranchaddress(TestTankStorage * const testtankstorage_, TTree * const tree_, bool mc);
void TestTankStorage_branch(TestTankStorage * const testtankstorage_, TTree * const tree_, bool mc);

////////////////////////////////////////////////////////////////////////////////
// PRINT TIMESTAMPS
////////////////////////////////////////////////////////////////////////////////

void timestamp_print(std::ostream& os, const TestTankStorage & testtankstorage_);

#endif
