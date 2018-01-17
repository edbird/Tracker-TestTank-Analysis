
#include "functions.hpp"

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
// PRINT FIT FUNCTION OUTPUT PARAMETERS
////////////////////////////////////////////////////////////////////////////////

void fit_param_print(std::ostream& os, TF1* func, int n_param)
{

    os << "chisq = " << f_feast_t0->GetChisquare() / func->GetNDF() << "\n";
    for(int px = 0; px < n_param, ) //++ px)
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

