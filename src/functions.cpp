
#include "functions.hpp"


////////////////////////////////////////////////////////////////////////////////
// HISTOGRAM DIFFERENTIATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

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
        Double_t content{(this_y - prev_y) / delta_t};
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

