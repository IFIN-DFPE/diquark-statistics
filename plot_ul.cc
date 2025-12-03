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

std::vector<double> masses;
std::vector<double> obs, med;
std::vector<double> sig1_lo, sig1_hi;
std::vector<double> sig2_lo, sig2_hi;

std::string path, process;
int discriminator;

/*
    Function taking as input the path of a .csv datafile and reads it
*/
void read_CSV(std::string inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (inputFile.empty()) throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");

    // Check if the file can be opened or not
    std::ifstream csvFile(inputFile);
    if (!csvFile.is_open()) throw std::runtime_error("Error: Please check the path of the input file.\n");

    std::cout << "Reading data from: " << inputFile << '\n';

    std::string line;

    // Skip header
    std::getline(csvFile, line);
    
    // Parse csv file by line
    while(std::getline(csvFile, line)) {
        std::stringstream str(line);
        std::string cell;

        // Mass values
        if(std::getline(str, cell, ',')) masses.push_back(stod(cell));
        // Observed median
        if(std::getline(str, cell, ',')) obs.push_back(stod(cell));
        // -2 sigma
        if(std::getline(str, cell, ',')) sig2_lo.push_back(stod(cell));
        // -1 sigma
        if(std::getline(str, cell, ',')) sig1_lo.push_back(stod(cell));
        // Expected median
        if(std::getline(str, cell, ',')) med.push_back(stod(cell));
        // +1 sigma
        if(std::getline(str, cell, ',')) sig1_hi.push_back(stod(cell));
        // +2 sigma
        if(std::getline(str, cell, ',')) sig2_hi.push_back(stod(cell));    

    }

    csvFile.close();

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


    TGraph* gObs = new TGraph(masses.size());
    TGraph* gMed = new TGraph(masses.size());

    TGraphAsymmErrors* gSig1 = new TGraphAsymmErrors(masses.size());
    TGraphAsymmErrors* gSig2 = new TGraphAsymmErrors(masses.size());

    for(int i = 0; i < masses.size(); i++) {
        gObs->SetPoint(i, masses[i], obs[i]);
        gMed->SetPoint(i, masses[i], med[i]);

        gSig1->SetPoint(i, masses[i], med[i]);
        gSig1->SetPointError(i, 0., 0., med[i]-sig1_lo[i], sig1_hi[i]-med[i]);

        gSig2->SetPoint(i, masses[i], med[i]);
        gSig2->SetPointError(i, 0., 0., med[i]-sig2_lo[i], sig2_hi[i]-med[i]);
    }

    gObs->SetLineColor(46);
    gObs->SetLineStyle(1);
    gObs->SetLineWidth(2);
    gObs->SetMarkerStyle(24);
    gObs->SetMarkerColor(46);
    gObs->SetMarkerSize(1.5);
    gMed->SetLineColor(kBlack);
    gMed->SetLineStyle(2);
    gMed->SetLineWidth(2);
    gSig1->SetFillColor(38);
    gSig2->SetTitle(";M_{S} [TeV];Upper Limit on #mu^{95}#times S_{ev}");
    gSig2->GetYaxis()->SetRangeUser(5e-1, 6e1);
    gSig2->SetFillColor(kOrange-4);

    TLegend* legend = new TLegend(0.6, 0.6, 0.8, 0.89);
    // legend->AddEntry(gObs, "Observed U.L.", "l");
    legend->AddEntry(gMed, "Expected U.L.", "l");
    legend->AddEntry(gSig1, "#pm1#sigma", "f");
    legend->AddEntry(gSig2, "#pm2#sigma", "f");
    legend->AddEntry((TObject*)0, "M_{#chi} = 1.5 TeV", "");
    legend->AddEntry((TObject*)0, "D = 0.9", "");
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);

    gSig2->Draw("A3 SAME");
    gSig1->Draw("3 SAME");
    gMed->Draw("L SAME");
    // gObs->Draw("PL SAME");
    legend->Draw("SAME");

    c->Update();
    c->Draw();

    std::string outPdf = path + "/graphs/" + process + "_upper_limit.pdf";
    std::string outPng = path + "/graphs/" + process + "_upper_limit.png";
    c->SaveAs(outPdf.c_str());
    c->SaveAs(outPng.c_str());

}

void plot_ul() {

    gROOT->SetBatch(1);

    // Load the file containing the paths
    std::ifstream pathFile("analysis_paths.json");
    nlohmann::json paths = nlohmann::json::parse(pathFile);

    // Load the configuration file
    std::ifstream configFile("config_plot.json");
    nlohmann::json config = nlohmann::json::parse(configFile);

    auto discriminator = int(config["discriminator_1"].get<float>()*1000);
    process = config["process_1"].get<std::string>();
    path = paths[process].get<std::string>();
    std::string inputFile = path + Form("/roostats_results/out_D%d/upper_limits.csv", discriminator);


    try{
        read_CSV(inputFile);
        generate_plot();
    
    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        exit(EXIT_FAILURE);
    }


}