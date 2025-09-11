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


std::vector<double> masses;
std::vector<double> p_vals;


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
        // Local p-value
        if(std::getline(str, cell, ',')) p_vals.push_back(stod(cell));   

    }

    csvFile.close();

    // Print the read data for debugging
    std::cout << "\n=== Summary by Mass Points ===\n";
    std::cout << "M_s [TeV]\t"
         << "p_val\n";
    for (int i = 0 ; i < masses.size(); i++)
        std::cout << std::fixed << std::setprecision(2) << masses[i] << "\t\t" << std::scientific << p_vals[i] << "\n";
    std::cout << "\n";

}


void generate_plot() {

    TCanvas* c = new TCanvas("c_pval", "Local p-values", 800, 600);
    c->SetLogy();

    TGraph* gP_val = new TGraph(masses.size(), &masses[0], &p_vals[0]);

    gP_val->SetTitle(";M_{S} [TeV];Local p-value");
    gP_val->SetMarkerStyle(24);
    gP_val->SetMarkerSize(1.5);
    gP_val->SetLineWidth(2);

    gP_val->Draw("APL");
    gP_val->GetXaxis()->SetRangeUser(6.8, 8.95);
    gP_val->GetYaxis()->SetRangeUser(1e-6, 1.);

    double x_min = 6.8;
    double x_max = 8.95;
    for (int sigma = 1; sigma <= 4; ++sigma) {
        double z = sigma;
        double pval = 0.5 * TMath::Erfc(z / TMath::Sqrt2());

        TLine *line = new TLine(x_min, pval, x_max, pval);
        line->SetLineStyle(4);
        line->SetLineColor(46);
        line->SetLineWidth(3);
        line->Draw();

        // Add label
        TLatex *label = new TLatex(x_max + 0.005, pval, Form("%d#sigma", sigma));
        label->SetTextSize(0.03);
        label->SetTextAlign(12);
        label->SetTextColor(46);
        label->Draw();
    }

    c->Update();
    c->Draw();


    if(std::filesystem::create_directories("results/mChi2/graphs"))
    ;
    c->SaveAs("results/mChi2/graphs/local_pvals.pdf");


}


void plot_pval() {

    gROOT->SetBatch(1);

    try{
        const char* inputFile = "results/mChi2/roofit_results/out_D900/p_values.csv";

        read_CSV(inputFile);
        generate_plot();
    }
    catch(const std::exception& exc) {
        std::cerr <<  exc.what() << '\n';
        exit(EXIT_FAILURE);
    }

}