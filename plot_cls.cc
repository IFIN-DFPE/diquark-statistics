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

std::vector<double> masses;
std::vector<double> medians;
std::vector<double> sig1_lo, sig1_hi;
std::vector<double> sig2_lo, sig2_hi;


void plot_cls() {

    gROOT->SetBatch(1);

    for(int i = 700; i <= 875; i += 25) {
        TFile* in_file = TFile::Open(Form("results/mChi2/roofit_results/out_D900/output_S%d.root", i), "READ");
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
    TLatex* label = new TLatex(masses.front()-0.1, 0.055, "95\% C.L.");

    gMedian->SetLineColor(kBlack);
    gMedian->SetLineStyle(9);
    gMedian->SetLineWidth(2);
    gSig1->SetFillColor(kGreen);
    gSig2->SetTitle("CL_{S} scan;M_{S} [TeV/c^{2}];CL_{S}");
    gSig2->SetFillColor(kYellow);
    line_lim->SetLineColor(kRed);
    line_lim->SetLineWidth(2);
    line_lim->SetLineStyle(1);
    label->SetTextSize(0.03);
    label->SetTextColor(kRed);

    gSig2->Draw("A3");
    gSig1->Draw("3 SAME");
    gMedian->Draw("L SAME");
    line_lim->Draw("SAME");
    label->Draw("SAME");


    TLegend* legend = new TLegend(0.15, 0.7, 0.4, 0.88);
    legend->AddEntry(gMedian, "Median expected CL_{S}", "l");
    legend->AddEntry(gSig1, "#pm1#sigma", "f");
    legend->AddEntry(gSig2, "#pm2#sigma", "f");
    legend->AddEntry(line_lim, "#alpha = 0.05", "l");
    legend->Draw();

    c->Update();
    c->Draw();
    if(std::filesystem::create_directories("results/mChi2/graphs"))
    ;
    c->SaveAs("results/mChi2/graphs/cls_scan.pdf");

}