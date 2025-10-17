#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>

#include "TGraph.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"

struct p_val_points{
    std::vector<double> masses;
    std::vector<double> p_vals;
};

p_val_points D900, D925;

/*
    Function taking as input the path of a .csv datafile and reads it
*/
p_val_points read_CSV(const char* inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (!inputFile) throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");

    // Check if the file can be opened or not
    std::ifstream csvFile(inputFile);
    if (!csvFile.is_open()) throw std::runtime_error("Error: Please check the path of the input file.\n");

    std::cout << "Reading data from: " << inputFile << '\n';

    std::string line;

    // Skip header
    std::getline(csvFile, line);
    
    p_val_points result;

    // Parse csv file by line
    while(std::getline(csvFile, line)) {
        std::stringstream str(line);
        std::string cell;

        // Mass values
        if(std::getline(str, cell, ',')) result.masses.push_back(stod(cell));
        // Local p-value
        if(std::getline(str, cell, ',')) result.p_vals.push_back(stod(cell));   

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

    TCanvas* c = new TCanvas("c_pval", "Local p-values", 800, 600);
    c->SetLogy();

    TGraph* gP_val_900 = new TGraph(D900.masses.size(), &D900.masses[0], &D900.p_vals[0]);
    TGraph* gP_val_925 = new TGraph(D925.masses.size(), &D925.masses[0], &D925.p_vals[0]);

    gP_val_900->SetTitle(";M_{S} [TeV];Local p-value");
    gP_val_900->SetMarkerStyle(24);
    gP_val_900->SetMarkerSize(1.5);
    gP_val_900->SetLineWidth(2);

    gP_val_900->Draw("APL");
    gP_val_900->GetXaxis()->SetRangeUser(6.8, 8.95);
    gP_val_900->GetYaxis()->SetRangeUser(5e-8, 1.);

    double x_min = 6.8;
    double x_max = 8.95;
    for (int sigma = 1; sigma <= 5; ++sigma) {
        double z = sigma;
        double pval = 0.5 * TMath::Erfc(z / TMath::Sqrt2());

        TLine *line = new TLine(x_min, pval, x_max, pval);
        line->SetLineStyle(4);
        line->SetLineColor(46);
        line->SetLineWidth(3);
        line->Draw("SAME");

        // Add label
        TLatex *label = new TLatex(x_max + 0.005, pval, Form("%d#sigma", sigma));
        label->SetTextSize(0.04);
        label->SetTextAlign(12);
        label->SetTextColor(46);
        label->Draw("SAME");
    }

    TLatex* chiMass = new TLatex(8.4, 1e-6, "M_{#chi} = 1.5 TeV");
    chiMass->SetTextSize(0.04);
    chiMass->Draw("SAME");

    c->Update();
    c->Draw();


    if(std::filesystem::create_directories("ChiChi/yuu_02/mChi1_5/graphs/out_D925"))
    ;
    c->SaveAs("ChiChi/yuu_02/mChi1_5/graphs/out_D925/local_pvals.pdf");
    c->SaveAs("ChiChi/yuu_02/mChi1_5/graphs/out_D925/local_pvals.png");


}


void plot_pval() {

    gROOT->SetBatch(1);

    try{
        const char* inputD900 = "ChiChi/yuu_02/mChi2/roofit_results/out_D900/p_values.csv";
        const char* inputD925 = "ChiChi/yuu_02/mChi2/roofit_results/out_D925/p_values.csv";

        D900 = read_CSV(inputD900);
        D925 = read_CSV(inputD925);
        generate_plot();
    }
    catch(const std::exception& exc) {
        std::cerr <<  exc.what() << '\n';
        exit(EXIT_FAILURE);
    }

}