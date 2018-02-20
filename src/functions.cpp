
#include "functions.hpp"


////////////////////////////////////////////////////////////////////////////////
// HISTOGRAM DIFFERENTIATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////


TH1F* histogram_create_copy_limits(TH1F* histo, const std::string& name)
{
    // get limits and number of bins
    Double_t xlow{histo->GetBinLowEdge(1)};
    Double_t xhigh{histo->GetBinLowEdge(histo->GetNbinsX()) + histo->GetBinWidth(histo->GetNbinsX())};
    Int_t nbinsx{histo->GetNbinsX()};
    // create new histogram
    return new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xhigh);
}


TH1F *histogram_create_copy_limits_average(TH1F* histo, const std::string& name, const Int_t level)
{
    // get limits and number of bins
    Double_t xlow{histo->GetBinLowEdge(1)};
    Double_t xhigh{histo->GetBinLowEdge(histo->GetNbinsX()) + histo->GetBinWidth(histo->GetNbinsX())};
    Int_t nbinsx{histo->GetNbinsX() / level};
    // create new histogram
    return new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xhigh);
}


void histogram_destroy(TH1F* histo)
{
    delete histo;
    histo = nullptr;
}


void histogram_average(TH1F* const output, TH1F* const input, const Int_t level)
{
    // TODO: bin range check
    Int_t ix_max{output->GetNbinsX()};
    if((ix_max - 1) * level + level > input->GetNbinsX())
    {
        std::cout << "correcting for overflow" << std::endl;
        -- ix_max; // TODO: can this still overflow?
        if(ix_max & level - 1 > input->GetNbinsX())
        {
            std::cout << "error in histogram_average" << std::endl;
        }
    }
    // TODO: fix this (bin input/output ranges) [re-enable couts]
    //std::cout << "input->GetNbinsX()=" << input->GetNbinsX() << std::endl;
    //std::cout << "output->GetNbinsX()=" << output->GetNbinsX() << std::endl;
    for(Int_t ix{1}; ix <= ix_max; ++ ix)
    {
        Double_t content{0.0};
        Double_t content_square{0.0};
        for(Int_t jx{1}; jx <= level; ++ jx)
        {
            Int_t input_bin{jx + (ix - 1) * level};
            //std::cout << "input_bin=" << input_bin << " ";
            content += input->GetBinContent(input_bin);
            content_square += input->GetBinContent(input_bin) * input->GetBinContent(input_bin);
        }
        //std::cout << "\noutput_bin=" << ix << std::endl;
        content /= (Double_t)(level);
        output->SetBinContent(ix, content);
        //output->SetBinError(ix, std::sqrt(content_square - content * content));
    }
}


void histogram_smooth(TH1F* const output, TH1F* const input, const Int_t level)
{
    if(input->GetNbinsX() != output->GetNbinsX()) throw "error histogram nbinsx differ";
    if(output->GetNbinsX() < level) throw "error histogram nbinsx < level";
    
    const Int_t half_level{level / 2};
    const Int_t level_minus_half_level{level - half_level};
    for(Int_t ix{1 + level_minus_half_level}; ix <= input->GetNbinsX() - half_level; ++ ix)
    {
        Double_t sum{0.0};
        for(Int_t jx{-level_minus_half_level}; jx <= half_level; ++ jx)
        {
            sum += input->GetBinContent(ix + jx);
        }
        sum /= (Double_t)level;
        output->SetBinContent(ix, sum);
    }
    for(Int_t ix{1}; ix < 1 - level_minus_half_level; ++ ix)
    {
        output->SetBinContent(ix, 0.0);
    }
    for(Int_t ix{input->GetNbinsX() - half_level + 1}; ix <= input->GetNbinsX(); ++ ix)
    {
        output->SetBinContent(ix, 0.0);
    }
}


void histogram_differentiate(TH1F* const differential_histo, const TH1F* const histo, Int_t method)
{
    if(method == differentiate_method_simple) histogram_differentiate_simple(differential_histo, histo);
    else if(method == differentiate_method_simple_with_smooth) histogram_differentiate_simple_with_smooth(differential_histo, histo);
    else if(method == differentiate_method_filter_simulation) histogram_differentiate_filter_simulation(differential_histo, histo);
    else throw "invalid differentiation method flag";
}


// trivial / simple differentation method
// dy/dt = Delta_y / Delta_t
// first order
void histogram_differentiate_simple(TH1F* const differential_histo, const TH1F* const histo)
{
    if(differential_histo->GetNbinsX() != histo->GetNbinsX()) throw "error histogram nbinsx differ";
    if(histo->GetNbinsX() < 2) throw "error histogram nbinsx < 2";

    Double_t prev_y{histo->GetBinContent(1)};
    const Double_t delta_t{histo->GetBinLowEdge(2) - histo->GetBinLowEdge(1)};
    for(Int_t ix{1 + 1}; ix <= histo->GetNbinsX(); ++ ix)
    {
        Double_t this_y{histo->GetBinContent(ix)};
        //Double_t content{(this_y - prev_y) / delta_t}; // TODO: delta_t here is small
        Double_t content{(this_y - prev_y)};
        differential_histo->SetBinContent(ix - 1, content);
        prev_y = this_y;
    }
    // set final bin to zero
    differential_histo->SetBinContent(histo->GetNbinsX(), 0.0);
}


// compute derivative using trivial differentiation method with smoothing
void histogram_differentiate_simple_with_smooth(TH1F* const differential_histo, const TH1F* const histo)
{
    // TODO: CHECK INDICES
    std::cout << "a" << std::endl;

    if(differential_histo->GetNbinsX() != histo->GetNbinsX()) throw "error histogram nbinsx differ";
    if(histo->GetNbinsX() < 2) throw "error histogram nbinsx < 2";

    // average of 3 numbers either side
    const Int_t average_count{3};
    Double_t sum{0};
    const Double_t delta_t{histo->GetBinLowEdge(2) - histo->GetBinLowEdge(1)};
    for(Int_t ix{1 + average_count}; ix <= histo->GetNbinsX() - average_count; ++ ix)
    {
        for(Int_t jx{1}; jx <= average_count; ++ jx)
        {
            sum += histo->GetBinContent(ix + jx) - histo->GetBinContent(ix - jx);
        }
        differential_histo->SetBinContent(ix, sum / (Double_t)(2 * average_count));
    }
    // set edge bins to zero
    for(Int_t ix{1}; ix < 1 + average_count; ++ ix)
    {
        differential_histo->SetBinContent(ix, 0.0);
    }
    for(Int_t ix{histo->GetNbinsX() - average_count + 1}; ix < histo->GetNbinsX(); ++ ix)
    {
        differential_histo->SetBinContent(ix, 0.0);
    }

    std::cout << "b" << std::endl;
}


// compute differential using high pass filter method
void histogram_differentiate_filter_simulation(TH1F* const differential_histo, const TH1F* const histo)
{

    // differentiator circuit
    // (high pass filter)
    double R; // RC filter resistance
    double C; // RC filter impedance
}


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

// Fit function for cathode time
Double_t cathode_time_fitf(Double_t *x_, Double_t *par)
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

// Plasma propagation time
Double_t ppt_fitf(Double_t *x_, Double_t *par)
{

    Double_t x = *x_;
    
    // parameters
    Double_t A{par[0]}; // amplitude
    Double_t t{par[1]}; // top start
    Double_t a{par[2]}; // turn on point
    Double_t b{par[3]}; // decay constant
    Double_t d{par[4]}; // constant
    //Double_t c{par[3]}; // x shift
    //Double_t d{par[4]}; // temp
    //Double_t l{par[6]}; // top length
    
    if(t <= x && x < a)
    {
        return A + d;
    }
    else if(a <= x)
    {
        //return A * std::exp(-b * (x - c)) + d;
        return A * std::exp(-b * (x - a)) + d;
    }
    else
    {
        return 0.0 + d;
    }
    
}

// Z Position Cathode Time
Double_t zpos_cathode_time_fitf(Double_t *x_, Double_t *par)
{
    Double_t A{par[0]}; // amplitude
    Double_t mean{par[1]};
    Double_t sigma{par[2]};
    Double_t theta{par[3]}; // rotation

    Double_t x{x_[0]};
    Double_t y{x_[1]};

    // Reject central region, and edges
    if(x <= -0.95)
    {
        TF1::RejectPoint();
        return 0;
    }
    else if(-0.05 <= x && x <= 0.05)
    {
        TF1::RejectPoint();
        return 0;
    }
    else if(0.95 <= x)
    {
        TF1::RejectPoint();
        return 0;
    }

    Double_t xx{x * std::cos(theta) - y * std::sin(theta)};
    return A * std::exp(-std::pow((xx - mean) / sigma, 2.0));

}

// Z Position Cathode Time Profile

// Z Position Cathode Time Profile (non-fit)
Double_t zpos_cathode_time_profilef(Double_t x, Double_t mean, Double_t theta)
{
    return (x * std::cos(theta) - mean) / std::sin(theta);
}

// Gaussian distribution
Double_t gaussian_fitf(Double_t *x_, Double_t *par)
{
    Double_t x{x_[0]};
    Double_t A{par[0]};
    Double_t mean{par[1]};
    Double_t sigma{par[2]};
    return A * std::exp(-std::pow((x - mean) / sigma, 2.0));
}

// Double Gaussian distribution
Double_t double_gaussian_fitf(Double_t *x_, Double_t *par)
{
    Double_t x{x_[0]};
    Double_t A{par[0]};
    Double_t mean{par[1]};
    Double_t sigma{par[2]};
    Double_t A2{par[3]};
    Double_t mean2{par[4]};
    Double_t sigma2{par[5]};
    return A * std::exp(-std::pow((x - mean) / sigma, 2.0)) + A2 * std::exp(-std::pow((x - mean2) / sigma2, 2.0));
}


////////////////////////////////////////////////////////////////////////////////
// TCANVAS FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

// TODO: get histogram name auto

void canvas(TH1F* const histogram_, const std::string& canvas_name_, const std::string& draw_opt_)
{
    TCanvas *c_c = new TCanvas(canvas_name_.c_str(), canvas_name_.c_str(), 800, 600);
    histogram_->Draw(draw_opt_.c_str());
    c_c->SaveAs((canvas_name_ + std::string(".C")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".eps")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".png")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".pdf")).c_str());
    histogram_->Write();
    delete c_c;

    return;
}

void canvas(TH2F* const histogram_, const std::string& canvas_name_)
{
    TCanvas *c_c = new TCanvas(canvas_name_.c_str(), canvas_name_.c_str(), 800, 600);
    histogram_->Draw("colz");
    c_c->SaveAs((canvas_name_ + std::string(".C")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".eps")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".png")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".pdf")).c_str());
    histogram_->Write();  
    delete c_c;

    return;
}

void canvas_fit(TH1F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_)
{
    TCanvas *c_c = new TCanvas(canvas_name_.c_str(), canvas_name_.c_str(), 800, 600);
    histogram_->Fit(fit_name_.c_str());
    histogram_->Draw();
    c_c->SaveAs((canvas_name_ + std::string(".C")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".eps")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".png")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".pdf")).c_str());
    histogram_->Write();
    delete c_c;

    return;
}

void canvas_fit(TH2F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_)
{
    TCanvas *c_c = new TCanvas(canvas_name_.c_str(), canvas_name_.c_str(), 800, 600);
    histogram_->Fit(fit_name_.c_str());
    histogram_->Draw("colz");
    c_c->SaveAs((canvas_name_ + std::string(".C")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".eps")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".png")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".pdf")).c_str());
    histogram_->Write();
    delete c_c;

    return;
}

void canvas_scale_fit(TH1F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_)
{
    TCanvas *c_c = new TCanvas(canvas_name_.c_str(), canvas_name_.c_str(), 800, 600);
    Double_t integral = histogram_->Integral();
    histogram_->Scale(1.0 / integral);
    histogram_->Fit(fit_name_.c_str());
    histogram_->Draw();
    c_c->SaveAs((canvas_name_ + std::string(".C")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".eps")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".png")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".pdf")).c_str());
    histogram_->Write();
    delete c_c;

    return;
}

void canvas_scale_fit(TH2F* const histogram_, const std::string& canvas_name_, const std::string& fit_name_)
{
    TCanvas *c_c = new TCanvas(canvas_name_.c_str(), canvas_name_.c_str(), 800, 600);
    Double_t integral = histogram_->Integral();
    histogram_->Scale(1.0 / integral);
    histogram_->Fit(fit_name_.c_str());
    histogram_->Draw();
    c_c->SaveAs((canvas_name_ + std::string(".C")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".eps")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".png")).c_str());
    c_c->SaveAs((canvas_name_ + std::string(".pdf")).c_str());
    histogram_->Write();
    delete c_c;

    return;
}


////////////////////////////////////////////////////////////////////////////////
// DIFFERENTIAL HISTOGRAM FIT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

Double_t get_timestamps_and_peaks(const TH1F* const histo, const Double_t VLN, const Double_t VHN, const Double_t VHP, Double_t &R0, Double_t &R1, Double_t &R2, Double_t &R3, Double_t &R4, Double_t &A0, Double_t& A1, Double_t& A2, Double_t& A3, Double_t& A4)
{
    // set to invalid values
    R0 = -1.0;
    R1 = -1.0;
    R2 = -1.0;
    R3 = -1.0;
    R4 = -1.0;

    bool have_R0{false};
    bool have_R1{false};
    bool have_R2{false};
    bool have_R3{false};
    bool have_R4{false};

    Int_t ix{1};
    // get R0
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(!have_R0)
        {
            if(content <= VLN)
            {
                R0 = histo->GetBinCenter(ix);
                A0 = content; // set A0 to some initial value
                //std::cout << "R0: VLN=" << VLN << " VALUE=" << content << " TIME=" << R0 << std::endl; 
                //++ ix;
                //break;
                have_R0 = true;
            }
        }
        else
        {
            // search for minimum A0, stop searching
            // when content reaches threshold for below loop
            if(content < A0)
            {
                A0 = content;
            }

            if(content > VLN)
            {
                //++ ix;
                break;
            }
        }
        ++ ix;
    }
    // seek to when signal passes VLN again
    //for(; ix <= histo->GetNbinsX(); )
    //{
    //    const Double_t content{histo->GetBinContent(ix)};
    //    if(content > VLN)
    //    {
    //        ++ ix;
    //        break;
    //    }
    //    ++ ix;
    //}
    // get R1
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(!have_R1)
        {
            if(content <= VHN)
            {
                R1 = histo->GetBinCenter(ix);
                A1 = content;
                have_R1 = true;
                //std::cout << "R1: VHN=" << VHN << " VALUE=" << content << " TIME=" << R1 << std::endl; 
                //++ ix;
                //break;
            }
        }
        else
        {
            if(content < A1)
            {
                A1 = content;
            }

            if(content > 0.0)
            {
                //++ ix; // TODO: should this be here?
                break;
            }
        }
        ++ ix;
    }
    // get R3
    for(; ix <= histo->GetNbinsX(); )
    {
        //std::cout << "enter to continue..." << std::endl;
        //std::cin.get();
        const Double_t content{histo->GetBinContent(ix)};
        if(!have_R3)
        {
            //std::cout << "don't have R3" << std::endl;
            if(content >= VHP)
            {
                //std::cout << ">= VHP" << std::endl;
                R3 = histo->GetBinCenter(ix);
                A3 = content;
                //std::cout << "setting A3=" << A3 << std::endl;
                have_R3 = true;
                //std::cout << "R3: VHP=" << VHP << " VALUE=" << content << " TIME=" << R3 << std::endl; 
                //++ ix;
                //break;
            }
        }
        else
        {
            //std::cout << "have R3" << std::endl;
            if(content > A3)
            {
                //std::cout << "A3=" << content << std::endl;
                A3 = content;
            }

            if(content < 0.0)
            {
                //std::cout << "break" << std::endl;
                //++ ix;
                break; // TODO
            }
        }
        ++ ix;
    }
    // get R2
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(!have_R2)
        {
            if(content <= VHN)
            {
                R2 = histo->GetBinCenter(ix);
                A2 = content;
                have_R2 = true;
                //std::cout << "R2: VHN=" << VHN << " VALUE=" << content << " TIME=" << R2 << std::endl; 
                //++ ix;
                //break;
            }
        }
        else
        {
            if(content < A2)
            {
                A2 = content;
            }

            if(content > 0.0)
            {
                //++ ix;
                break; // TODO
            }
        }
        ++ ix;
    }
    // get R4
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(!have_R4)
        {
            if(content >= VHP)
            {
                R4 = histo->GetBinCenter(ix);
                A4 = content;
                have_R4 = true;
                //std::cout << "R4: VLN=" << VHP << " VALUE=" << content << " TIME=" << R4 << std::endl; 
                //++ ix;
                //break;
            }
        }
        else
        {
            if(content > A4)
            {
                A4 = content;
            }

            // no break statement here - search to end
        }
        ++ ix;
    }
    // R0 to R4 are either invalid or set to valid values

    //std::cout << "R0: VLN=" << VLN << " VALUE=" << A0 << " TIME=" << R0 << std::endl; 
    //std::cout << "R1: VHN=" << VHN << " VALUE=" << A1 << " TIME=" << R1 << std::endl; 
    //std::cout << "R3: VHP=" << VHP << " VALUE=" << A3 << " TIME=" << R3 << std::endl; 
    //std::cout << "R2: VHN=" << VHN << " VALUE=" << A2 << " TIME=" << R2 << std::endl; 
    //std::cout << "R4: VLN=" << VHP << " VALUE=" << A4 << " TIME=" << R4 << std::endl; 
}


// this is not robust to noise! TODO
Double_t get_timestamps(const TH1F* const histo, const Double_t VLN, const Double_t VHN, const Double_t VHP, Double_t &R0, Double_t &R1, Double_t &R2, Double_t &R3, Double_t &R4)
{
    // set to invalid values
    R0 = -1.0;
    R1 = -1.0;
    R2 = -1.0;
    R3 = -1.0;
    R4 = -1.0;

    Int_t ix{1};
    // get R0
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(content <= VLN)
        {
            R0 = histo->GetBinCenter(ix);
            std::cout << "R0: VLN=" << VLN << " VALUE=" << content << " TIME=" << R0 << std::endl; 
            ++ ix;
            break;
        }
        ++ ix;
    }
    // seek to when signal passes VLN again
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(content > VLN)
        {
            ++ ix;
            break;
        }
        ++ ix;
    }
    // get R1
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(content <= VHN)
        {
            R1 = histo->GetBinCenter(ix);
            std::cout << "R1: VHN=" << VHN << " VALUE=" << content << " TIME=" << R1 << std::endl; 
            ++ ix;
            break;
        }
        ++ ix;
    }
    // get R3
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(content >= VHP)
        {
            R3 = histo->GetBinCenter(ix);
            std::cout << "R3: VHP=" << VHP << " VALUE=" << content << " TIME=" << R3 << std::endl; 
            ++ ix;
            break;
        }
        ++ ix;
    }
    // get R2
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(content <= VHN)
        {
            R2 = histo->GetBinCenter(ix);
            std::cout << "R2: VHN=" << VHN << " VALUE=" << content << " TIME=" << R2 << std::endl; 
            ++ ix;
            break;
        }
        ++ ix;
    }
    // get R4
    for(; ix <= histo->GetNbinsX(); )
    {
        const Double_t content{histo->GetBinContent(ix)};
        if(content >= VHP)
        {
            R4 = histo->GetBinCenter(ix);
            std::cout << "R4: VLN=" << VHP << " VALUE=" << content << " TIME=" << R4 << std::endl; 
            ++ ix;
            break;
        }
        ++ ix;
    }
    // R0 to R4 are either invalid or set to valid values

}


////////////////////////////////////////////////////////////////////////////////
// FIT FUNCTIONS FOR ANODE HISTOGRAM
////////////////////////////////////////////////////////////////////////////////

Double_t anode_fitf(Double_t *x_, Double_t *p_)
{

    Double_t x{*x_};
    Double_t *p{p_};
    Int_t p_ix{0};

    Double_t A{p[p_ix ++]}; // start
    Double_t B{p[p_ix ++]}; // first discontinuous change (cathode 1)
    Double_t C{p[p_ix ++]}; // second discontinuous change (cathode 2)
    Double_t D{p[p_ix ++]}; // end

    // first function
    Double_t A0{p[p_ix ++]};
    Double_t k0{p[p_ix ++]};
    Double_t x0{p[p_ix ++]};
    
    // second function
    Double_t A1{p[p_ix ++]};
    Double_t k1{p[p_ix ++]};
    Double_t x1{p[p_ix ++]};
    Double_t C1{p[p_ix ++]};
    
    // third function
    Double_t A2{p[p_ix ++]};
    Double_t k2{p[p_ix ++]};
    Double_t x2{p[p_ix ++]};

    if(A <= x && x < B)
    {
        //std::cout << "A0=" << A0 << std::endl;
        return anode_fitf_0(x, A0, k0, x0/*, -1.0*/); // C0 = 1.0
    }
    else if(B <= x && x < C)
    {
        return anode_fitf_1(x, A1, k1, x1, C1);
    }
    else if(C <= x && x < D)
    {
        //std::cout << "evaluating between C and D" << std::endl;
        //std::cout << "x=" << x << std::endl;
        Double_t sum{anode_fitf_2(x, A2, k2, x2, -C1)};
        //std::cout << "sum=" << sum << std::endl;
        sum += anode_fitf_1(x, A1, k1, x1, C1);
        //std::cout << "sum=" << sum << std::endl;
        //std::cin.get();
        return sum;
    }
    else
    {
        //std::cout << "ERROR" << std::endl;
    }
    return 0.0;

}

Double_t anode_fitf_0(Double_t x, Double_t A0, Double_t k0, Double_t x0/*, Double_t C0*/)
{
    //std::cout << A0 << " " << k0 << " " << x0 << std::endl;
    Double_t arg{k0 * (x - x0)};
    //std::cout << "arg=" << arg << std::endl;
    //std::cout << "exp=" << std::exp(-arg) << std::endl;
    //return A0 * (std::exp(-arg) + C0);
    //return A0 * (std::exp(-arg) - 1.0);
    return A0 * (1.0 - std::exp(-arg));
}

Double_t anode_fitf_1(Double_t x, Double_t A1, Double_t k1, Double_t x1, Double_t C1)
{
    Double_t arg{k1 * (x - x1)};
    //return A1 * (std::exp(-arg) + C1);
    return A1 * (std::exp(-arg)) + C1;
}

// TODO: oscillatory version
Double_t anode_fitf_2(Double_t x, Double_t A2, Double_t k2, Double_t x2, Double_t C2)
{
    Double_t arg{k2 * (x - x2)};
    //return A2 * (std::exp(-arg) + C2); 
    return A2 * (std::exp(-arg) /*- 1.0*/) + C2; 
}


////////////////////////////////////////////////////////////////////////////////
// FIT FUNCTIONS FOR DIFFERENTIAL HISTOGRAM
////////////////////////////////////////////////////////////////////////////////

Double_t differential_fitf(Double_t *x_, Double_t *par)
{
    Int_t par_ix{0};
    Double_t A{par[par_ix]}; ++ par_ix; // first peak point
    Double_t B{par[par_ix]}; ++ par_ix; // before second peak point
    Double_t C{par[par_ix]}; ++ par_ix; // second peak point
    Double_t D{par[par_ix]}; ++ par_ix; // before third peak point
    Double_t E{par[par_ix]}; ++ par_ix; // third peak point

    // TODO: change from piecewise to continuous
    Double_t x{x_[0]};
    
    Double_t x0{A};
    Double_t A0{par[par_ix]}; ++ par_ix;
    Double_t k0{par[par_ix]}; ++ par_ix;
    //std::cout << "x0=" << x0 << " A0=" << A0 << " k0=" << k0 << std::endl;
    
    Double_t x1{C};
    Double_t A1{par[par_ix]}; ++ par_ix;
    Double_t k1{par[par_ix]}; ++ par_ix;
    //std::cout << "x1=" << x1 << " A1=" << A1 << " k1=" << k1 << std::endl;

    Double_t x2{C};
    Double_t A2{par[par_ix]}; ++ par_ix;
    Double_t k2{par[par_ix]}; ++ par_ix;
    //std::cout << "x2=" << x2 << " A2=" << A2 << " k2=" << k2 << std::endl;

    Double_t x3{E};
    Double_t A3{par[par_ix]}; ++ par_ix;
    Double_t k3{par[par_ix]}; ++ par_ix;
    //std::cout << "x3=" << x3 << " A3=" << A3 << " k3=" << k3 << std::endl;

    Double_t x4{E};
    Double_t A4{par[par_ix]}; ++ par_ix;
    Double_t k4{par[par_ix]}; ++ par_ix;
    //std::cout << "x4=" << x4 << " A4=" << A4 << " k4=" << k4 << std::endl;

    /*
    if((A <= x) && (x < B))
    {
        Double_t exp_decay_par[3] = {x0, A0, k0};
        return exp_decay(x_, exp_decay_par);
    }
    else if((B <= x) && (x < C))
    {
        Double_t exp_decay_par[3] = {x1, A1, k1};
        return exp_decay(x_, exp_decay_par);
    }
    else if((C <= x) && (x < D))
    {
        Double_t exp_decay_par[3] = {x2, A2, k2};
        return exp_decay(x_, exp_decay_par);
    }
    else if((D <= x) && (x < E))
    {
        Double_t exp_decay_par[3] = {x3, A3, k3};
        return exp_decay(x_, exp_decay_par);
    }
    else if((E <= x) && (x < 120.0)) // TODO: should not be 120.0
    {
        Double_t exp_decay_par[3] = {x4, A4, k4};
        return exp_decay(x_, exp_decay_par);
    }
    return 0.0;
    */
    
    // change to continuous functions between A - C and C - E
    if((A <= x) && (x < C))
    {
        Double_t exp_decay_par_0[3] = {x0, A0, k0};
        Double_t exp_decay_par_1[3] = {x1, A1, k1};
        Double_t exp_decay_0 = exp_decay(x_, exp_decay_par_0);
        Double_t exp_decay_1 = exp_decay(x_, exp_decay_par_1);
        return exp_decay_0 + exp_decay_1;
    }
    else if((C <= x) && (x < E))
    {
        Double_t exp_decay_par_2[3] = {x2, A2, k2};
        Double_t exp_decay_par_3[3] = {x3, A3, k3};
        Double_t exp_decay_2 = exp_decay(x_, exp_decay_par_2);
        Double_t exp_decay_3 = exp_decay(x_, exp_decay_par_3);
        return exp_decay_2 + exp_decay_3;
    }
    else if((E <= x) && (x < 120.0)) // TODO: should not be 120.0
    {
        Double_t exp_decay_par[3] = {x4, A4, k4};
        return exp_decay(x_, exp_decay_par);
    }
    return 0.0;
}


Double_t exp_decay(Double_t *x_, Double_t *par)
{
    Double_t x{x_[0]};
    Double_t x0{par[0]};
    Double_t A0{par[1]};
    Double_t k0{par[2]};
    return A0 * std::exp(k0 * (x - x0));
}



////////////////////////////////////////////////////////////////////////////////
// PRINT FIT FUNCTION OUTPUT PARAMETERS
////////////////////////////////////////////////////////////////////////////////

void fit_param_print(std::ostream& os, TF1* func)
{
    int n_param = func->GetNpar();
    
    os << "chisq = " << func->GetChisquare() / func->GetNDF() << "\n";
    for(int px = 0; px < n_param; ) //++ px)
    {
        os << "p" << px << " = " << func->GetParameter(px) << " +- " << func->GetParError(px);
    
        if(++ px < n_param)
            os << "\n";
    }

}

void fit_param_print(std::ostream& os, TF2* func)
{
    int n_param = func->GetNpar();

    os << "chisq = " << func->GetChisquare() / func->GetNDF() << "\n";
    for(int px = 0; px < n_param; ) //++ px)
    {
        os << "p" << px << " = " << func->GetParameter(px) << " +- " << func->GetParError(px);
    
        if(++ px < n_param)
            os << "\n";
    }

}

////////////////////////////////////////////////////////////////////////////////
// WAVEFORM OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

std::string waveform_print(TH1F* histo, const Long64_t canvas_name_counter, const std::string& canvas_name, const std::string& canvas_directory, const std::string& draw_opt_, const Int_t width, const Int_t height)
{

    if(histo != nullptr)
    {
        // print out the waveform of the t0 events for inspection
        //const std::string canvas_name_base("cathode_event_");
        //const std::string canvas_name_base(output_name);
        //std::string canvas_name(canvas_name_base + int_to_string(canvas_name_counter, 6));
        std::string canvas_name_base(canvas_name + int_to_string(canvas_name_counter, 6));
        TCanvas *c = new TCanvas(canvas_name_base.c_str(), canvas_name_base.c_str(), width, height);
        histo->Draw(draw_opt_.c_str());
        std::string output_file_name(std::string("./") + canvas_directory + std::string("/") + canvas_name_base);
        c->SaveAs((output_file_name + std::string(".png")).c_str());
        delete c;
        //++ canvas_name_counter;
    
        return output_file_name;
    }

}

