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

/*
    Function taking as input the path of a .csv datafile and reads it
*/
void read_CSV(const char* inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (!inputFile) throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");

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
    gSig2->SetTitle(";M_{S} [TeV];Model independent 95\% C.L. limits on #mu");
    gSig2->GetYaxis()->SetRangeUser(0.05,6);
    gSig2->SetFillColor(kOrange-4);

    TLegend* legend = new TLegend(0.15, 0.7, 0.4, 0.88);
    legend->AddEntry(gObs, "Observed U.L.", "l");
    legend->AddEntry(gMed, "Expected U.L.", "l");
    legend->AddEntry(gSig1, "#pm1#sigma", "f");
    legend->AddEntry(gSig2, "#pm2#sigma", "f");
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);


    TLatex* chiMass = new TLatex(7., 1.4, "M_{#chi} = 2.0 TeV");
    chiMass->SetTextSize(0.04);

    gSig2->Draw("A3");
    gSig1->Draw("3 SAME");
    gMed->Draw("L SAME");
    gObs->Draw("PL SAME");
    // c->SetGrid();
    // c->RedrawAxis("g");
    legend->Draw("SAME");
    chiMass->Draw("SAME");

    c->Update();
    c->Draw();

    if(std::filesystem::create_directories("results/mChi2/graphs/out_D900"))
    ;
    c->SaveAs("results/mChi2/graphs/out_D900/chi2_upper_limit.pdf");

}

void plot_ul() {

    gROOT->SetBatch(1);

    try{
        const char* inputFile = "results/mChi2/roostats_results/out_D900/upper_limits.csv";

        read_CSV(inputFile);
        generate_plot();
    
    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        exit(EXIT_FAILURE);
    }


}