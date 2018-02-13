#include "datatypes.hpp"

void TestTankStorage_init(TestTankStorage * const testtankstorage_)
{

    TestTankStorage &store = *testtankstorage_;
     
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

    // added 2018-02-13
    store.differential_histo = (TH1F*)0;

    // this only for MC data
    store.truth_position = 0;
}


void TestTankStorage_setbranchaddress(TestTankStorage * const testtankstorage_, TTree * const tree_, bool mc = false)
{

    TTree * const t = tree_;
    TestTankStorage &store = *testtankstorage_;

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
    t->SetBranchAddress("anode_histo", &store.anode_histo);
    t->SetBranchAddress("deriv_histo", &store.deriv_histo);
    t->SetBranchAddress("cathode_histo", &store.cathode_histo);

    // only for mc data
    if(mc == true)
    {
        t->SetBranchAddress("truth_position", &store.truth_position);
    }
}

void TestTankStorage_branch(TestTankStorage * const testtankstorage_, TTree * const tree_, bool mc = false)
{

    TTree * const t2 = tree_;
    TestTankStorage &store = *testtankstorage_;

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

    // only for mc data
    if(mc == true)
    {
        t2->Branch("truth_position", &store.truth_position);
    }
}

////////////////////////////////////////////////////////////////////////////////
// PRINT TIMESTAMPS
////////////////////////////////////////////////////////////////////////////////

void timestamp_print(std::ostream& os, const TestTankStorage & testtankstorage_)
{

    const TestTankStorage &store = testtankstorage_;
    
    os << "cathode               : " << store.cathode_time << "\n"\
       << "t0                    : " << store.feast_t0 << "\n"\
       << "t1                    : " << store.feast_t1 << "\n"\
       << "t3                    : " << store.feast_t3 << "\n"\
       << "t2                    : " << store.feast_t2 << "\n"\
       << "t4                    : " << store.feast_t4 << "\n"\
       << "cathode + t0          : " << store.cathode_time + store.feast_t0 << "\n"\
       << "(t1 - t0) - cathode   : " << store.feast_t1 - store.feast_t0 - store.cathode_time << "\n"\
       << "(t3 - t0) - cathode   : " << store.feast_t3 - store.feast_t0 - store.cathode_time << "\n"\
       << "(t2 - t0) - cathode   : " << store.feast_t2 - store.feast_t0 - store.cathode_time << "\n"\
       << "(t4 - t0) - cathode   : " << store.feast_t4 - store.feast_t0 - store.cathode_time << "\n";

    /*
    std::cout << "cathode               : " << store.cathode_time << "\n"\
      << "t0                    : " << store.feast_t0 << "\n"\
      << "t1                    : " << store.feast_t1 << "\n"\
      << "t3                    : " << store.feast_t3 << "\n"\
      << "t2                    : " << store.feast_t2 << "\n"\
      << "t4                    : " << store.feast_t4 << "\n"\
      << "cathode + t0          : " << store.cathode_time + store.feast_t0 << "\n"\
      << "(t1 - t0) - cathode   : " << store.feast_t1 - store.feast_t0 - store.cathode_time << "\n"\
      << "(t3 - t0) - cathode   : " << store.feast_t3 - store.feast_t0 - store.cathode_time << "\n"\
      << "(t2 - t0) - cathode   : " << store.feast_t2 - store.feast_t0 - store.cathode_time << "\n"\
      << "(t4 - t0) - cathode   : " << store.feast_t4 - store.feast_t0 - store.cathode_time << "\n";
  */
  
}


