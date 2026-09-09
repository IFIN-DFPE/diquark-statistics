#include <chrono>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>

#include "nlohmann/json.hpp"
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

ul_point pointLimit1;
DataPoint pointYield1;
std::string path1, process1;
int d1;

DataPoint read_yields(std::string inputFile) {
    // Check if user has entered the path to the data file
    if(inputFile.empty()) {
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
ul_point read_CSV(std::string inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (inputFile.empty()) throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");

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

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetLabelSize(0.04, "XYZ");

    TCanvas* c = new TCanvas("c_banana", "Upper Limit", 800, 600);
    c->SetLogy();

    // Normal
    TGraph* gLim1 = new TGraph(pointLimit1.masses.size(), &pointLimit1.masses[0], &pointLimit1.med[0]);
    TGraph* gYld1 = new TGraph(pointLimit1.masses.size(), &pointYield1.m_s[0], &pointYield1.sig[0]);
    TGraphAsymmErrors* gSig1 = new TGraphAsymmErrors(pointLimit1.masses.size());
    TGraphAsymmErrors* gSig2 = new TGraphAsymmErrors(pointLimit1.masses.size());    
    for(int i = 0; i < pointLimit1.masses.size(); i++) {        
        gSig1->SetPoint(i, pointLimit1.masses[i], pointLimit1.med[i]);
        gSig1->SetPointError(i, 0., 0., pointLimit1.med[i]-pointLimit1.sig1_lo[i], pointLimit1.sig1_hi[i]-pointLimit1.med[i]);

        gSig2->SetPoint(i, pointLimit1.masses[i], pointLimit1.med[i]);
        gSig2->SetPointError(i, 0., 0., pointLimit1.med[i]-pointLimit1.sig2_lo[i], pointLimit1.sig2_hi[i]-pointLimit1.med[i]);
    
    }


    gYld1->SetTitle(";M_{S} [TeV];Event counts");
    gYld1->GetXaxis()->SetTitleOffset(1.2);
    // gYld1->GetXaxis()->SetNdivisions(310);
    gYld1->SetMarkerStyle(22);
    gYld1->SetMarkerSize(1.5);
    gYld1->SetLineWidth(2);
    gYld1->SetLineStyle(9);
    gYld1->SetLineColor(kRed-4);
    gYld1->SetMarkerColor(kRed-4);

    gLim1->SetMarkerStyle(22);
    gLim1->SetMarkerSize(1.5);
    gLim1->SetLineWidth(2);
    gLim1->SetLineStyle(7);
    gLim1->SetLineColor(kBlue+1);
    gLim1->SetMarkerColor(kBlue+1);

    gSig1->SetFillColor(38);
    gSig2->SetTitle(";M_{S} [TeV];CL_{s}");
    gSig2->GetXaxis()->SetTitleOffset(1.2);
    gSig2->SetFillColor(kOrange-4);

    gYld1->GetYaxis()->SetRangeUser(3e-1, 1e2); // Normal
    // gYld1->GetYaxis()->SetRangeUser(1, 20); // Discriminator limits
    // gYld1->GetXaxis()->SetRangeUser(7.35, 8.65); // Discriminator limits


    TLegend* legend = new TLegend(0.7, 0.6, 0.9, 0.9);
    legend->AddEntry(gLim1, "#mu^{95}#times S_{ev}, u(h^{0}t)", "pl");
    legend->AddEntry(gYld1, "S_{ev}, u(h^{0}t)", "pl");       
    legend->AddEntry((TObject*)0, "m_{#chi} = 2.0 TeV", "");
    // legend->AddEntry((TObject*)0, "y_{uu} = 0.2", "");

    legend->SetTextSize(0.03);
    // legend->SetFillStyle(0);
    // legend->SetFillColor(0);
    // legend->SetBorderSize(0);

    gSig2->GetYaxis()->SetRangeUser(1e0, 5e1);
    gSig2->Draw("A3");
    gSig1->Draw("3 SAME");
    gYld1->Draw("C SAME");
    gLim1->Draw("L SAME");
    legend->Draw("SAME");
    

    c->Update();
    c->Draw();

    std::string outPdf = path1 + "/graphs/" + process1 + "_upper_limit_errors.pdf";
    std::string outPng = path1 + "/graphs/" + process1 + "_upper_limit_errors.png";
    c->SaveAs(outPdf.c_str());
    c->SaveAs(outPng.c_str());

}

void plot_limits_errors() {

    gROOT->SetBatch(1);

    // Load the file containing the paths
    std::ifstream pathFile("analysis_paths.json");
    nlohmann::json paths = nlohmann::json::parse(pathFile);

    // Load the configuration file
    std::ifstream configFile("config_plot.json");
    nlohmann::json config = nlohmann::json::parse(configFile);

    auto d1 = int(config["discriminator_1"].get<float>()*1000);
    process1 = config["process_1"].get<std::string>();
    path1 = paths[process1].get<std::string>();
    std::string inputLim1 = path1 + Form("/roostats_results/out_D%d/upper_limits.csv", d1);
    std::string inputYld1 = path1 + Form("/signal_yields/sig_bkg_D%d.csv", d1);

    try{
        pointLimit1 = read_CSV(inputLim1);
        pointYield1 = read_yields(inputYld1);

        generate_plot();
    
    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        exit(EXIT_FAILURE);
    }


}