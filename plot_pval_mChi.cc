#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>

#include "nlohmann/json.hpp"
#include "TGraph.h"
#include "TObject.h"
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

p_val_points chi1_5, chi2;
std::string path1, process1, path2, process2;
int d1, d2;

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

    TGraph* gP_val_chi15 = new TGraph(chi1_5.masses.size(), &chi1_5.masses[0], &chi1_5.p_vals[0]);
    TGraph* gP_val_chi2 = new TGraph(chi2.masses.size(), &chi2.masses[0], &chi2.p_vals[0]);

    gP_val_chi15->SetTitle(";M_{S} [TeV];Local p-value");
    gP_val_chi15->GetXaxis()->SetTitleOffset(1.2);
    gP_val_chi2->SetMarkerStyle(22);
    gP_val_chi2->SetMarkerSize(1.5);
    gP_val_chi2->SetLineWidth(2);
    gP_val_chi2->SetLineColor(kBlue+1);
    gP_val_chi2->SetMarkerColor(kBlue+1);
    
    gP_val_chi15->SetMarkerStyle(20);
    gP_val_chi15->SetMarkerSize(1.5);
    gP_val_chi15->SetLineWidth(2);
    gP_val_chi15->SetLineColor(kBlack);
    gP_val_chi15->SetMarkerColor(kBlack);

    gP_val_chi15->GetXaxis()->SetRangeUser(6.8, 8.95);
    gP_val_chi15->GetYaxis()->SetRangeUser(1e-8, 1.);

    gP_val_chi15->Draw("APL");
    gP_val_chi2->Draw("PL SAME");


    double x_min = gP_val_chi15->GetXaxis()->GetXmin();
    double x_max = gP_val_chi15->GetXaxis()->GetXmax();
    for (int sigma = 1; sigma <= 5; ++sigma) {
        double z = sigma;
        double pval = 0.5 * TMath::Erfc(z / TMath::Sqrt2());

        TLine *line = new TLine(x_min, pval, x_max, pval);
        line->SetLineStyle(4);
        line->SetLineColor(46);
        line->SetLineWidth(3);
        line->Draw("SAME");

        // Add label
        TLatex *label = new TLatex(x_max + 0.025, pval, Form("%d#sigma", sigma));
        label->SetTextSize(0.04);
        label->SetTextAlign(12);
        label->SetTextColor(46);
        label->Draw("SAME");
    }

    TLegend* legend = new TLegend(0.65, 0.25, 0.85, 0.43);
    legend->AddEntry(gP_val_chi15, "m_{#chi} = 1.5 TeV", "pl");
    legend->AddEntry(gP_val_chi2, "m_{#chi} = 2 TeV", "pl");
    legend->AddEntry((TObject*)0, "D = 0.9", "");
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);
    legend->Draw("SAME");

    c->Update();
    c->Draw();

    std::string outPdf = path1 + "/graphs/" + process1 + "_local_p_vals_mChi.pdf";
    std::string outPng = path1 + "/graphs/" + process1 + "_local_p_vals_mChi.png";
    c->SaveAs(outPdf.c_str());
    c->SaveAs(outPng.c_str());
}


void plot_pval_mChi() {

    gROOT->SetBatch(1);

    // Load the file containing the paths
    std::ifstream pathFile("analysis_paths.json");
    nlohmann::json paths = nlohmann::json::parse(pathFile);

    // Load the configuration file
    std::ifstream configFile("config_plot.json");
    nlohmann::json config = nlohmann::json::parse(configFile);

    d1 = int(config["discriminator_1"].get<float>()*1000);
    d2 = int(config["discriminator_2"].get<float>()*1000);
    process1 = config["process_1"].get<std::string>();
    path1 = paths[process1].get<std::string>();
    std::string input1 = path1 + Form("/roofit_results/out_D%d/p_values.csv", d1);
    process2 = config["process_2"].get<std::string>();
    path2 = paths[process2].get<std::string>();
    std::string input2 = path2 + Form("/roofit_results/out_D%d/p_values.csv", d2);

    try{
        chi1_5 = read_CSV(input2.c_str());
        chi2 = read_CSV(input1.c_str());
        generate_plot();
    }
    catch(const std::exception& exc) {
        std::cerr <<  exc.what() << '\n';
        exit(EXIT_FAILURE);
    }

}