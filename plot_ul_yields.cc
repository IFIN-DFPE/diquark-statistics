#include <chrono>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>

#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TObject.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"

struct ul_point {
    std::vector<double> masses;
    std::vector<double> obs, med;
    std::vector<double> sig1_lo, sig1_hi;
    std::vector<double> sig2_lo, sig2_hi;
};

struct DataPoint {
    std::vector<double> m_s;
    std::vector<double> sig, sigma_sig;
    std::vector<double> bkg, sigma_bkg;
};

ul_point yuu02_limit, yuu04_limit;
DataPoint yuu02_yield, yuu04_yield;

DataPoint read_yields(const char* inputFile) {
    // Check if user has entered the path to the data file
    if(!inputFile) {
        throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");
    }

    // Check if the file can be opened or not
    std::ifstream csvFile(inputFile);
    if (!csvFile.is_open()) {
        throw std::runtime_error("Error: Please check the path of the input file.\n");
    }

    std::cout << "Reading data from: " << inputFile << '\n';

    DataPoint result;
    std::string line;

    // Skip header
    getline(csvFile, line);

    // Parse csv file by line
    while(getline(csvFile, line)) {
        std::stringstream str(line);
        std::string cell;

        // Mass values
        if(getline(str, cell, ',')) result.m_s.push_back(stod(cell));
        // Signal values
        if(getline(str, cell, ',')) result.sig.push_back(stod(cell));
        // Signal uncertainties
        if(getline(str, cell, ',')) result.sigma_sig.push_back(stod(cell));
        // Background values
        if(getline(str, cell, ',')) result.bkg.push_back(stod(cell));
        // Background uncertainties
        if(getline(str, cell, ',')) result.sigma_bkg.push_back(stod(cell));

    }

    csvFile.close();

    return result;
}


/*
    Function taking as input the path of a .csv datafile and reads it
*/
ul_point read_CSV(const char* inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (!inputFile) throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");

    // Check if the file can be opened or not
    std::ifstream csvFile(inputFile);
    if (!csvFile.is_open()) throw std::runtime_error("Error: Please check the path of the input file.\n");

    std::cout << "Reading data from: " << inputFile << '\n';

    std::string line;

    ul_point result;

    // Skip header
    std::getline(csvFile, line);
    
    // Parse csv file by line
    while(std::getline(csvFile, line)) {
        std::stringstream str(line);
        std::string cell;

        // Mass values
        if(std::getline(str, cell, ',')) result.masses.push_back(stod(cell));
        // Observed median
        if(std::getline(str, cell, ',')) result.obs.push_back(stod(cell));
        // -2 sigma
        if(std::getline(str, cell, ',')) result.sig2_lo.push_back(stod(cell));
        // -1 sigma
        if(std::getline(str, cell, ',')) result.sig1_lo.push_back(stod(cell));
        // Expected median
        if(std::getline(str, cell, ',')) result.med.push_back(stod(cell));
        // +1 sigma
        if(std::getline(str, cell, ',')) result.sig1_hi.push_back(stod(cell));
        // +2 sigma
        if(std::getline(str, cell, ',')) result.sig2_hi.push_back(stod(cell));    

    }

    csvFile.close();

    return result;

}

void generate_plot() {

    gROOT->SetBatch(1);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetLabelSize(0.04, "XYZ");

    TCanvas* c = new TCanvas("c_banana", "Upper Limit", 800, 600);
    c->SetLogy();


    TGraph* yuu02_lim = new TGraph(yuu02_limit.masses.size(), &yuu02_limit.masses[0], &yuu02_limit.med[0]);
    TGraph* yuu02_yld = new TGraph(yuu02_yield.m_s.size(), &yuu02_yield.m_s[0], &yuu02_yield.sig[0]);
    TGraph* yuu04_lim = new TGraph(yuu04_limit.masses.size(), &yuu04_limit.masses[0], &yuu04_limit.med[0]);
    TGraph* yuu04_yld = new TGraph(yuu04_yield.m_s.size(), &yuu04_yield.m_s[0], &yuu04_yield.sig[0]);


    yuu04_yld->SetTitle(";M_{S} [TeV];Event counts");
    yuu04_yld->SetMarkerStyle(22);
    yuu04_yld->SetMarkerSize(1.5);
    yuu04_yld->SetLineWidth(2);
    yuu04_yld->SetLineColor(kMagenta+1);
    yuu04_yld->SetMarkerColor(kMagenta+1);

    yuu04_lim->SetMarkerStyle(22);
    yuu04_lim->SetMarkerSize(1.5);
    yuu04_lim->SetLineWidth(2);
    yuu04_lim->SetLineStyle(2);
    yuu04_lim->SetLineColor(kMagenta+1);
    yuu04_lim->SetMarkerColor(kMagenta+1);
    
    yuu02_yld->SetMarkerStyle(21);
    yuu02_yld->SetMarkerSize(1.5);
    yuu02_yld->SetLineWidth(2);
    yuu02_yld->SetLineColor(kBlue+1);
    yuu02_yld->SetMarkerColor(kBlue+1);

    yuu02_lim->SetMarkerStyle(21);
    yuu02_lim->SetMarkerSize(1.5);
    yuu02_lim->SetLineWidth(2);
    yuu02_lim->SetLineStyle(2);
    yuu02_lim->SetLineColor(kBlue+1);
    yuu02_lim->SetMarkerColor(kBlue+1);


    yuu04_yld->GetYaxis()->SetRangeUser(5e-1, 3e2);
    yuu04_yld->GetXaxis()->SetRangeUser(6.75, 10.);
    yuu04_yld->GetXaxis()->SetNdivisions(515);

    TLine *line = new TLine(8.17, 5e-1, 8.17, 3e2);
    line->SetLineStyle(9);
    line->SetLineColor(kRed-4);
    line->SetLineWidth(2);

    TLine *line1 = new TLine(9.16, 5e-1, 9.16, 3e2);
    line1->SetLineStyle(10);
    line1->SetLineColor(kRed-4);
    line1->SetLineWidth(2);



    TLegend* legend = new TLegend(0.5, 0.7, 0.9, 0.9);
    legend->AddEntry(yuu02_lim, "#mu^{95}#times S_{ev}, M_{#chi} = 2 TeV, y_{uu} = 0.2", "pl");
    legend->AddEntry(yuu02_yld, "S_{ev}, M_{#chi} = 2 TeV, y_{uu} = 0.2", "pl");
    legend->AddEntry(yuu04_lim, "#mu^{95}#times S_{ev}, M_{#chi} = 2 TeV, y_{uu} = 0.4", "pl");
    legend->AddEntry(yuu04_yld, "S_{ev}, M_{#chi} = 2 TeV, y_{uu} = 0.4", "pl");
    legend->SetTextSize(0.03);
    // legend->SetFillStyle(0);
    // legend->SetFillColor(0);
    // legend->SetBorderSize(0);

    yuu04_yld->Draw("APL");
    yuu04_lim->Draw("PL SAME");
    yuu02_yld->Draw("PL SAME");
    yuu02_lim->Draw("PL SAME");
    line->Draw("SAME");
    line1->Draw("SAME");
    legend->Draw("SAME");
    

    c->Update();
    c->Draw();

    c->SaveAs("ChiChi/yuu_02/mChi2/graphs/new_chan_upper_limit_yields.pdf");
    c->SaveAs("ChiChi/yuu_02/mChi2/graphs/new_chan_upper_limit_yields.png");

}

void plot_ul_yields() {

    gROOT->SetBatch(1);

    try{
        const char* inputLimYuu02 = "ChiChi/yuu_02/mChi2/roostats_results/out_D900/upper_limits.csv";
        const char* inputLimYuu04 = "ChiChi/yuu_04/mChi2/roostats_results/out_D900/upper_limits.csv";
        const char* inputYldYuu02 = "ChiChi/yuu_02/mChi2/signal_yields/sig_bkg_D900.csv";
        const char* inputYldYuu04 = "ChiChi/yuu_04/mChi2/signal_yields/sig_bkg_D900.csv";
        yuu02_limit = read_CSV(inputLimYuu02);
        yuu04_limit = read_CSV(inputLimYuu04);
        yuu02_yield = read_yields(inputYldYuu02);
        yuu04_yield = read_yields(inputYldYuu04);

        generate_plot();
    
    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        exit(EXIT_FAILURE);
    }


}