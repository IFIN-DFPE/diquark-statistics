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
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"

std::vector<double> masses;
std::vector<double> medians;
std::vector<double> sig1_lo, sig1_hi;
std::vector<double> sig2_lo, sig2_hi;

std::string path, process;

void plot_cls() {

    // Load the file containing the paths
    std::ifstream pathFile("analysis_paths.json");
    nlohmann::json paths = nlohmann::json::parse(pathFile);

    // Load the configuration file
    std::ifstream configFile("config_plot.json");
    nlohmann::json config = nlohmann::json::parse(configFile);
    auto discriminator = int(config["discriminator_1"].get<float>()*1000);
    process = config["process_1"].get<std::string>();
    path = paths[process].get<std::string>();



    gROOT->SetBatch(1);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetLabelSize(0.04, "XYZ");

    for(int i = 700; i <= 875; i += 25) {
        std::string inputPath = path + Form("/roofit_results/out_D%d/output_S%d.root", discriminator, i);
        TFile* in_file = TFile::Open(inputPath.c_str(), "READ");
        masses.push_back(double(i)/100);

        TH1F* hCLS = (TH1F*)in_file->Get("hCLS");
        double quantiles[5] = {0.0228, 0.1587, 0.5, 0.8413, 0.9772};
        double expected_cls[5];
        hCLS->GetQuantiles(5, expected_cls, quantiles);

        sig2_lo.push_back(expected_cls[0]);
        sig1_lo.push_back(expected_cls[1]);
        medians.push_back(expected_cls[2]);
        sig1_hi.push_back(expected_cls[3]);
        sig2_hi.push_back(expected_cls[4]);
    }

    TCanvas* c = new TCanvas("c_banana", "CLS", 800, 600);

    TGraph* gMedian = new TGraph(masses.size());
    
    TGraphAsymmErrors* gSig1 = new TGraphAsymmErrors(masses.size());
    TGraphAsymmErrors* gSig2 = new TGraphAsymmErrors(masses.size());
    

    for(int i = 0; i < masses.size(); i++) {
        gMedian->SetPoint(i, masses[i], medians[i]);
        
        gSig1->SetPoint(i, masses[i], medians[i]);
        gSig1->SetPointError(i, 0., 0., medians[i]-sig1_lo[i], sig1_hi[i]-medians[i]);

        gSig2->SetPoint(i, masses[i], medians[i]);
        gSig2->SetPointError(i, 0., 0., medians[i]-sig2_lo[i], sig2_hi[i]-medians[i]);
    
    }

    TLine* line_lim = new TLine(gMedian->GetXaxis()->GetXmin(), 0.05, gMedian->GetXaxis()->GetXmax(), 0.05);
    TLatex* label = new TLatex(masses.front()-0.05, 0.055, "95\% C.L.");

    gMedian->SetLineColor(kBlack);
    gMedian->SetLineStyle(2);
    gMedian->SetLineWidth(2);
    gSig1->SetFillColor(38);
    gSig2->SetTitle(";M_{S} [TeV];CL_{s}");
    gSig2->GetXaxis()->SetTitleOffset(1.2);
    gSig2->SetFillColor(kOrange-4);
    line_lim->SetLineColor(46);
    line_lim->SetLineWidth(2);
    line_lim->SetLineStyle(1);
    label->SetTextColor(46);
    label->SetTextSize(0.04);

    gSig2->GetYaxis()->SetRangeUser(0., 1.);
    gSig2->Draw("A3");
    gSig1->Draw("3 SAME");
    gMedian->Draw("L SAME");
    line_lim->Draw("SAME");
    label->Draw("SAME");


    TLegend* legend = new TLegend(0.15, 0.6, 0.4, 0.88);
    legend->AddEntry(gMedian, "Median expected CL_{s}", "l");
    legend->AddEntry(gSig1, "#pm1#sigma", "f");
    legend->AddEntry(gSig2, "#pm2#sigma", "f");
    legend->AddEntry(line_lim, "#alpha = 0.05", "l");
    legend->AddEntry((TObject*)0, "m_{#chi} = 1.5 TeV", "");
    legend->AddEntry((TObject*)0, "y_{uu} = 0.2", "");
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);

    legend->Draw();
    c->Update();
    c->Draw();
    
    std::string outputPdf = path + "/graphs/" + process + "_cls_scan.pdf";
    std::string outputPng = path + "/graphs/" + process + "_cls_scan.png";
    c->SaveAs(outputPdf.c_str());
    c->SaveAs(outputPng.c_str());

}