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

ul_point best, worst;

/*
    Function taking as input the path of a .csv datafile and reads it
*/
struct ul_point read_CSV(const char* inputFile) {
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


    TGraph* gMed_best = new TGraph(best.masses.size(), &best.masses[0], &best.med[0]);
    TGraph* gMed_worst = new TGraph(worst.masses.size(), &worst.masses[0], &worst.med[0]);

    gMed_best->SetTitle(";M_{S} [TeV];Expected Upper Limit on #mu^{95}#times S_{ev}");
    gMed_best->SetMarkerStyle(25);
    gMed_best->SetMarkerSize(1.5);
    gMed_best->SetLineWidth(2);
    gMed_best->SetLineColor(38);
    gMed_best->SetMarkerColor(38);
    
    gMed_worst->SetMarkerStyle(24);
    gMed_worst->SetMarkerSize(1.5);
    gMed_worst->SetLineWidth(2);

    gMed_best->GetYaxis()->SetRangeUser(5e-1, 6e1);

    TLegend* legend = new TLegend(0.4, 0.65, 0.7, 0.85);
    legend->AddEntry(gMed_best, "M_{#chi} = 1.5 TeV, y_{uu} = 0.4, D = 0.925", "pl");
    legend->AddEntry(gMed_worst, "M_{#chi} = 2 TeV, y_{uu} = 0.2, D = 0.9", "pl");
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);

    gMed_best->Draw("APL");
    gMed_worst->Draw("PL SAME");
    legend->Draw("SAME");

    c->Update();
    c->Draw();

    c->SaveAs("ChiChi/yuu_02/mChi2/graphs/upper_limit_bw.pdf");
    c->SaveAs("ChiChi/yuu_02/mChi2/graphs/upper_limit_bw.png");

}

void plot_ul_bw() {

    gROOT->SetBatch(1);

    try{
        const char* inputBest = "ChiChi/yuu_04/mChi1_5/roostats_results/out_D925/upper_limits.csv";
        const char* inputWorst = "ChiChi/yuu_02/mChi2/roostats_results/out_D900/upper_limits.csv";
        best = read_CSV(inputBest);
        worst = read_CSV(inputWorst);
        generate_plot();
    
    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        exit(EXIT_FAILURE);
    }


}