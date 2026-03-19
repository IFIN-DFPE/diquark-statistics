#include "nlohmann/json.hpp"

#include "TCanvas.h"
#include "TF1.h"    
#include "TLine.h"
#include "TROOT.h"

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooRandom.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooLognormal.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"  

using namespace RooStats;
using namespace RooFit;

void combined_example() {

    gROOT->SetBatch(1);
    
    RooRandom::randomGenerator()->SetSeed(42);

    RooWorkspace* wspace = new RooWorkspace("wspace");
    
    RooRealVar mass("mass", "mass", 7.5, 7., 10.);
    RooUniform sig_shape("sig_shape", "sig_shape", mass);
    RooUniform bkg_shape("bkg_shape", "bkg_shape", mass);

    RooRealVar mu("mu", "signal multiplier", 1.0, 0.0, 10.);
    // mu.setConstant(true);


    // uChi channel 1
    RooRealVar S0_obs_1("S0_obs_1", "Signal yield from simulation first channel", 15.80408, .5, 100.);
    S0_obs_1.setConstant();
    RooRealVar S0_true_1("S0_true_1", "True signal yield first channel", 15.80408, .5, 100.);
    RooRealVar sigmaS_1("sigmaS_1", "Std dev of true signal yield first channel", 1. + 0.02425918 / 15.80408, 1.0001, 100.);
    sigmaS_1.setConstant();
    RooLognormal constraint_yield_sig_1("constraint_yield_sig_1", "constraint_yield_sig_1", S0_obs_1, S0_true_1, sigmaS_1);
    RooFormulaVar total_signal_1("total_signal_1", "mu*S0_true_1", {mu, S0_true_1});

    RooRealVar B0_obs_1("B0_obs_1", "Background yield from simulation first channel", 8.5966, .5, 1000.);
    B0_obs_1.setConstant();
    RooRealVar B0_true_1("B0_true_1", "True Background yield first channel", 8.5966, .5, 1000.);
    RooRealVar sigmaB_1("sigmaB_1", "std dev of background yield first channel", 1.0 + 0.96473 / 8.5966, 1.01, 100.);
    sigmaB_1.setConstant();
    RooLognormal constraint_yield_bkg_1("constraint_yield_bkg_1", "Yield constraint shape for bkg first channel", B0_obs_1, B0_true_1, sigmaB_1);

    RooExtendPdf ext_sig_1("ext_sig_1", "Extended Signal PDF first channel", sig_shape, total_signal_1);
    RooExtendPdf ext_bkg_1("ext_bkg_1", "Extended Bkg PDF first channel", bkg_shape, B0_true_1);

    RooAddPdf sb_tmp_1("sb_tmp_1", "sb_tmp_1", {ext_sig_1, ext_bkg_1});
    RooProdPdf sb_full_1("sb_full_1", "sb_full_1", {sb_tmp_1, constraint_yield_sig_1, constraint_yield_bkg_1});


    // ChiChi channel 2
    RooRealVar S0_obs_2("S0_obs_2", "Signal yield from simulation second channel", 18.30226, 1E-6, 100.);
    RooRealVar S0_true_2("S0_true_2", "True signal yield second channel", 18.30226, 1E-6, 100.);
    RooRealVar sigmaS_2("sigmaS_2", "Std dev of true signal yield second channel", 1. + 0.01462763 / 18.30226, 1.0001, 100.);
    sigmaS_2.setConstant();
    RooLognormal constraint_yield_sig_2("constraint_yield_sig_2", "constraint_yield_sig_2", S0_obs_2, S0_true_2, sigmaS_2);
    RooFormulaVar total_signal_2("total_signal_2", "mu*S0_true_2", {mu, S0_true_2});

    RooRealVar B0_obs_2("B0_obs_2", "Background yield from simulation second channel", 8.004, 1E-6, 1000.);
    RooRealVar B0_true_2("B0_true_2", "True Background yield second channel", 8.004, 1E-6, 1000.);
    RooRealVar sigmaB_2("sigmaB_2", "std dev of background yield second channel", 1.0 + 4.64 / 8.004, 1.01, 100.);
    sigmaB_2.setConstant();
    RooLognormal constraint_yield_bkg_2("constraint_yield_bkg_2", "Yield constraint shape for bkg second channel", B0_obs_2, B0_true_2, sigmaB_2);

    RooExtendPdf ext_sig_2("ext_sig_2", "Extended Signal PDF second channel", sig_shape, total_signal_2);
    RooExtendPdf ext_bkg_2("ext_bkg_2", "Extended Bkg PDF second channel", bkg_shape, B0_true_2);

    RooAddPdf sb_tmp_2("sb_tmp_2", "sb_tmp_2", {ext_sig_2, ext_bkg_2});
    RooProdPdf sb_full_2("sb_full_2", "sb_full_2", {sb_tmp_2, constraint_yield_sig_2, constraint_yield_bkg_2});


    RooCategory channels("channels", "channels");
    channels.defineType("channel1");
    channels.defineType("channel2");

    RooSimultaneous combined("combined", "combined", channels);
    combined.addPdf(sb_full_1, "channel1");
    combined.addPdf(sb_full_2, "channel2");


    ModelConfig* sbModel = new ModelConfig("sbModel", wspace);
    sbModel->SetPdf(combined);
    sbModel->SetParametersOfInterest(mu);
    sbModel->SetObservables({mass, channels});
    sbModel->SetNuisanceParameters({S0_true_1, B0_true_1, S0_true_2, B0_true_2});
    sbModel->SetGlobalObservables({S0_obs_1, B0_obs_1, S0_obs_2, B0_obs_2});
    sbModel->SetSnapshot(mu);
    
    ModelConfig* bModel = new ModelConfig("bModel", wspace);
    bModel->SetPdf(combined);
    RooRealVar bPoi = RooRealVar("bPoi", "signal multiplier in b-only model", 0.);
    bPoi.setConstant();
    bModel->SetParametersOfInterest(bPoi);
    bModel->SetNuisanceParameters({S0_true_1, B0_true_1, S0_true_2, B0_true_2});
    bModel->SetObservables({mass, channels});
    bModel->SetGlobalObservables({S0_obs_1, B0_obs_1, S0_obs_2, B0_obs_2});    
    bModel->SetSnapshot(bPoi);


    mu.setVal(0.0);
    mu.setConstant(kTRUE);
    RooDataSet *data = combined.generate(RooArgSet(mass, channels), 5000);
    mu.setConstant(kFALSE);


    AsymptoticCalculator* asympCalc = new AsymptoticCalculator(*data, *bModel, *sbModel);
    asympCalc->SetOneSided(true);
    asympCalc->SetPrintLevel(0);

    HypoTestInverter* inverter = new HypoTestInverter(*asympCalc);
    inverter->SetConfidenceLevel(0.95);
    inverter->UseCLs(true);
    inverter->SetVerbose(false);
    inverter->SetFixedScan(50, 0, 6);

    HypoTestInverterResult* result = inverter->GetInterval();

    TCanvas* c2 = new TCanvas("CLs_mu95", "CLs_mu95", 800, 600);
    HypoTestInverterPlot* plot = new HypoTestInverterPlot("HTI_Result_Plot", "95\% Upper Limit scan on #mu;#mu;CL_{S}", result);
    plot->Draw("CLB 2CL");
    c2->Draw();
    c2->SaveAs("combined/CLs_mu95.pdf");


}