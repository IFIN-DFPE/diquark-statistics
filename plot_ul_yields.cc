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

ul_point pointLimit1, pointLimit2, pointLimit3;
DataPoint pointYield1, pointYield2, pointYield3;
std::string path1, process1, path2, process2, path3, process3;
int d1, d2, d3;

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


    TGraph* gLim1 = new TGraph(pointLimit1.masses.size(), &pointLimit1.masses[0], &pointLimit1.med[0]);
    TGraph* gLim2 = new TGraph(pointLimit2.masses.size(), &pointLimit2.masses[0], &pointLimit2.med[0]);
    TGraph* gLim3 = new TGraph(pointLimit3.masses.size(), &pointLimit3.masses[0], &pointLimit3.med[0]);
    TGraph* gYld1 = new TGraph(pointYield1.m_s.size(), &pointYield1.m_s[0], &pointYield1.sig[0]);
    TGraph* gYld2 = new TGraph(pointYield2.m_s.size(), &pointYield2.m_s[0], &pointYield2.sig[0]);
    TGraph* gYld3 = new TGraph(pointYield3.m_s.size(), &pointYield3.m_s[0], &pointYield3.sig[0]);


    gYld1->SetTitle(";M_{S} [TeV];Event counts");
    gYld1->GetXaxis()->SetTitleOffset(1.2);
    // gYld1->GetXaxis()->SetNdivisions(310);
    gYld1->SetMarkerStyle(22);
    gYld1->SetMarkerSize(1.5);
    gYld1->SetLineWidth(2);
    gYld1->SetLineColor(kBlue+1);
    gYld1->SetMarkerColor(kBlue+1);

    gLim1->SetMarkerStyle(22);
    gLim1->SetMarkerSize(1.5);
    gLim1->SetLineWidth(2);
    gLim1->SetLineStyle(7);
    gLim1->SetLineColor(kBlue+1);
    gLim1->SetMarkerColor(kBlue+1);

    gYld2->SetMarkerStyle(20);
    gYld2->SetMarkerSize(1.5);
    gYld2->SetLineWidth(2);
    gYld2->SetLineColor(kBlack);
    gYld2->SetMarkerColor(kBlack);

    gLim2->SetMarkerStyle(20);
    gLim2->SetMarkerSize(1.5);
    gLim2->SetLineWidth(2);
    gLim2->SetLineStyle(7);
    gLim2->SetLineColor(kBlack);
    gLim2->SetMarkerColor(kBlack);

    gYld3->SetMarkerStyle(25);
    gYld3->SetMarkerSize(1.5);
    gYld3->SetLineWidth(2);
    gYld3->SetLineColor(kRed-4);
    gYld3->SetMarkerColor(kRed-4);

    gLim3->SetMarkerStyle(25);
    gLim3->SetMarkerSize(1.5);
    gLim3->SetLineWidth(2);
    gLim3->SetLineStyle(8);
    gLim3->SetLineColor(kRed-4);
    gLim3->SetMarkerColor(kRed-4);
    


    gYld1->GetYaxis()->SetRangeUser(5e-1, 3e2);
    // gYld1->GetXaxis()->SetRangeUser(8.75, 9.75);

    TLine *line1 = new TLine(8.15, 5e-1, 8.15, 8e1);
    // TLine *line1 = new TLine(9.26, 5e-1, 9.26, 3e2);
    // TLine *line1 = new TLine(9.17, 5e-1, 9.17, 3e2);
    line1->SetLineStyle(9);
    line1->SetLineColor(kRed-4);
    line1->SetLineWidth(2);

    // TLine *line2 = new TLine(9.26, 1e0, 9.26, 4.5e0);
    TLine *line2 = new TLine(8.18, 5e-1, 8.18, 3e2);
    // TLine *line2 = new TLine(8.25, 5e-1, 8.25, 8e1);
    line2->SetLineStyle(10);
    line2->SetLineColor(kRed-4);
    line2->SetLineWidth(2);

    TLine *line3 = new TLine(9.282, 1e0, 9.282, 4.5e0);
    // TLine *line2 = new TLine(9.26, 5e-1, 9.26, 3e2);
    line3->SetLineStyle(8);
    line3->SetLineColor(kRed-4);
    line3->SetLineWidth(2);


    TLegend* legend = new TLegend(0.65, 0.65, 0.9, 0.9);
    legend->AddEntry(gLim1, "#mu^{95}#times S_{ev}, y_{u#chi} = 0.5", "pl");
    legend->AddEntry(gYld1, "S_{ev}, y_{u#chi} = 0.5", "pl");    
    // legend->AddEntry(gLim2, "#mu^{95}#times S_{ev}, y_{uu} = 0.2", "pl");
    // legend->AddEntry(gYld2, "S_{ev}, y_{uu} = 0.2", "pl");
    // legend->AddEntry(gLim3, "#mu^{95}#times S_{ev}, D = 0.925", "pl");
    // legend->AddEntry(gYld3, "S_{ev}, D = 0.925", "pl");    
    legend->AddEntry((TObject*)0, "m_{#chi} = 2.0 TeV", "");
    // legend->AddEntry((TObject*)0, "y_{uu} = 0.4", "");

    legend->SetTextSize(0.03);
    // legend->SetFillStyle(0);
    // legend->SetFillColor(0);
    // legend->SetBorderSize(0);

    gYld1->Draw("APL");
    gLim1->Draw("PL SAME");
    // gYld2->Draw("PL SAME");
    // gLim2->Draw("PL SAME");
    // gYld3->Draw("PL SAME");
    // gLim3->Draw("PL SAME");
    line1->Draw("SAME");
    // line2->Draw("SAME");
    // line3->Draw("SAME");
    legend->Draw("SAME");
    

    c->Update();
    c->Draw();

    std::string outPdf = path1 + "/graphs/" + process1 + "_upper_limit_yields.pdf";
    std::string outPng = path1 + "/graphs/" + process1 + "_upper_limit_yields.png";
    c->SaveAs(outPdf.c_str());
    c->SaveAs(outPng.c_str());

}

void plot_ul_yields() {

    gROOT->SetBatch(1);

    // Load the file containing the paths
    std::ifstream pathFile("analysis_paths.json");
    nlohmann::json paths = nlohmann::json::parse(pathFile);

    // Load the configuration file
    std::ifstream configFile("config_plot.json");
    nlohmann::json config = nlohmann::json::parse(configFile);

    auto d1 = int(config["discriminator_1"].get<float>()*1000);
    auto d2 = int(config["discriminator_2"].get<float>()*1000);
    auto d3 = int(config["discriminator_3"].get<float>()*1000);
    process1 = config["process_1"].get<std::string>();
    process2 = config["process_2"].get<std::string>();
    process3 = config["process_3"].get<std::string>();
    path1 = paths[process1].get<std::string>();
    path2 = paths[process2].get<std::string>();
    path3 = paths[process3].get<std::string>();
    std::string inputLim1 = path1 + Form("/roostats_results/out_D%d/upper_limits.csv", d1);
    std::string inputLim2 = path2 + Form("/roostats_results/out_D%d/upper_limits.csv", d2);
    std::string inputLim3 = path3 + Form("/roostats_results/out_D%d/upper_limits.csv", d3);
    std::string inputYld1 = path1 + Form("/signal_yields/sig_bkg_D%d.csv", d1);
    std::string inputYld2 = path2 + Form("/signal_yields/sig_bkg_D%d.csv", d2);
    std::string inputYld3 = path3 + Form("/signal_yields/sig_bkg_D%d.csv", d3);

    try{
        pointLimit1 = read_CSV(inputLim1);
        pointLimit2 = read_CSV(inputLim2);
        pointLimit3 = read_CSV(inputLim3);
        pointYield1 = read_yields(inputYld1);
        pointYield2 = read_yields(inputYld2);
        pointYield3 = read_yields(inputYld3);

        generate_plot();
    
    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        exit(EXIT_FAILURE);
    }


}