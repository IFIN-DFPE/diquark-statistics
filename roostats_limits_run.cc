/*
    Run the program with the path to the .csv as an argument
    E.g.: root 'roostats_analysis.cc("path/to/file/.csv")'
*/

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooPoisson.h"
#include "RooDataSet.h"
#include "RooEffProd.h"
#include "RooAbsReal.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooUniform.h"
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooRandom.h"
#include "RooMinimizer.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "TRandom.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/ToyMCSampler.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>

using namespace std;
using namespace RooFit;
using namespace RooStats;


/*
    Data type holding:
    - the mass value
    - the theoretical signal yield and its uncertainty
    - the theoretical background yield and its uncertainty 
*/
struct DataPoint {
    double m_s; 
    double sig, sigma_sig;
    double bkg, sigma_bkg;
};


/*
    Function taking as input the path of a .csv datafile and 
    returns a vector of DataPoint type values.
*/
vector<DataPoint> read_CSV(const char* inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (!inputFile) {
        cerr << "Error: Please enter the name of the data file to be read.\n";
        cerr << "Usage: root \'roostats_analysis(\"filename\")\'\n";
        exit(EXIT_FAILURE);
    }

    // Check if the file can be opened or not
    ifstream csvFile(inputFile);
    if (!csvFile.is_open()) {
        cerr << "Error: Please check the path of the input file.\n";
        exit(EXIT_FAILURE);
    }

    cout << "Reading data from: " << inputFile << '\n';

    vector<DataPoint> data;
    string line;

    // Skip header
    getline(csvFile, line);
    
    // Parse csv file by line
    while(getline(csvFile, line)) {
        stringstream str(line);
        string cell;
        DataPoint point;

        // Mass values
        if(getline(str, cell, ',')) point.m_s = stod(cell);
        // Signal values
        if(getline(str, cell, ',')) point.sig = stod(cell);
        // Signal uncertainties
        if(getline(str, cell, ',')) point.sigma_sig = stod(cell);
        // Background values
        if(getline(str, cell, ',')) point.bkg = stod(cell);
        // Background uncertainties
        if(getline(str, cell, ',')) point.sigma_bkg = stod(cell);
            
        data.push_back(point);
    }

    csvFile.close();

    // Check if the program actually read something
    if(data.empty()) {
        cerr << "Error: The file contains no data.\n";
        exit(EXIT_FAILURE);
    }

    // Print the read data for debugging
    cout << "\n=== Summary by Mass Points ===\n";
    cout << "M_s [TeV]\t"
         << "SIG\t\t\t" 
         << "BKG\n";
    for (const DataPoint point : data)
        cout << fixed << setprecision(2) << point.m_s << "\t\t" << scientific 
             << point.sig << " ± " << point.sigma_sig << "\t"
             << point.bkg << " ± " << point.sigma_bkg << "\n";
    cout << "\n";

    return data;
}



/*
    Function taking as input a DataPoint variable and computes the 
    exclusion limit for the specified point.

    Returns a vector of double type, containing:
        - exclusion_limits[0] = 95% upper limit
        - exclusion_limits[1] = expected -2 sigma limit
        - exclusion_limits[2] = expected -1 sigma limit
        - exclusion_limits[3] = expected median limit
        - exclusion_limits[4] = expected +1 sigma limit
        - exclusion_limits[5] = expected +2 sigma limit
*/
vector<double> point_exclusion(DataPoint point, TFile* output_file) {
    vector<double> exclusion_limits;

    RooRandom::randomGenerator()->SetSeed(69);

    // Create workspace
    RooWorkspace* wspace = new RooWorkspace("wspace");
    
    // Define a dummy variable for mass and create the shapes for sig and bkg
    RooRealVar mass("mass", "mass", point.m_s, 7., 10.);
    RooUniform sig_shape("sig_shape", "sig_shape", mass);
    RooUniform bkg_shape("bkg_shape", "bkg_shape", mass);

// "Observed" yield from the simulation
    RooRealVar S0_obs("S0_obs", "Signal yield from simulation", point.sig, 1E-6, 100.);
    // "True" yield from fitting
    RooRealVar S0_true("S0_true", "True signal yield", point.sig, 1E-6, 100.);
    // Uncertainty of the yield, using lognormal so we don't run into numerical issues
    RooRealVar sigmaS("sigmaS", "Std dev of true signal yield", 1. + point.sigma_sig / point.sig, 1.0001, 100.);
    sigmaS.setConstant();
    RooLognormal constraint_yield_sig("constraint_yield_sig", "constraint_yield_sig", S0_obs, S0_true, sigmaS);

    
    // Signal strength multiplier
    RooRealVar mu("mu", "signal multiplier", 1.0, 0.0, 10.);
    mu.setConstant(true);
    // Complete signal expression
    RooFormulaVar total_signal("total_signal", "mu*S0_true", RooArgList(mu, S0_true));


    // "Observed" yield from the simulation
    RooRealVar B0_obs("B0_obs", "Background yield from simulation", point.bkg, 1E-6, 1000.);
    // "True" yield from fitting
    RooRealVar B0_true("B0_true", "True Background yield", point.bkg, 1E-6, 1000.);
    // Uncertainty of the yield, using lognormal so we don't run into numerical issues
    RooRealVar sigmaB("sigmaB", "std dev of background yield", 1.0 + point.sigma_bkg / point.bkg, 1.01, 100.);
    sigmaB.setConstant();
    RooLognormal constraint_yield_bkg("constraint_yield_bkg", "Yield constraint shape for bkg", B0_obs, B0_true, sigmaB);

    // Extend the pdfs over the entire mass region
    RooExtendPdf ext_sig("ext_sig", "Extended Signal PDF", sig_shape, total_signal);
    RooExtendPdf ext_bkg("ext_bkg", "Extedned Bkg PDF", bkg_shape, B0_true);

    // Add up pdfs into the full model
    RooAddPdf sb_tmp("sb_tmp", "sb_tmp", RooArgSet(ext_sig, ext_bkg));
    RooProdPdf sb_full("sb_full", "sb_full", RooArgSet(sb_tmp, constraint_yield_sig, constraint_yield_bkg));
    

    // Create s+b model configuration
    ModelConfig* sbModel = new ModelConfig("sbModel", wspace);
    // Add parameters to the model configuration
    sbModel->SetPdf(sb_full);
    sbModel->SetParametersOfInterest(mu);
    sbModel->SetNuisanceParameters({S0_true, B0_true});
    sbModel->SetObservables(mass);
    sbModel->SetGlobalObservables({S0_obs, B0_obs});
    sbModel->SetSnapshot(mu);


    // Create b-only model configuration 
    ModelConfig* bModel = new ModelConfig("bModel", wspace);
    bModel->SetPdf(sb_full);
    RooRealVar bPoi = mu;
    bPoi.setVal(0.);
    bPoi.setConstant(kTRUE);
    bModel->SetParametersOfInterest(bPoi);
    bModel->SetNuisanceParameters({S0_true, B0_true});
    bModel->SetObservables(mass);
    bModel->SetGlobalObservables({S0_obs, B0_obs});
    bModel->SetSnapshot(bPoi);

    // Create dataset       
    RooDataSet* toyData = sb_full.generate(RooArgSet(S0_obs, B0_obs), 1);

    // Create asymptotic calculator 
    AsymptoticCalculator* asympCalc = new AsymptoticCalculator(*toyData, *bModel, *sbModel);
    asympCalc->SetOneSided(true);
    asympCalc->SetPrintLevel(0);
    // asympCalc->GenerateAsimovData(sb_full, RooArgSet(S0_obs, B0_obs));

    HypoTestInverter* inverter = new HypoTestInverter(*asympCalc);
    inverter->SetConfidenceLevel(0.95);
    inverter->UseCLs(true);
    inverter->SetVerbose(false);
    inverter->SetFixedScan(100, 0, 6); 
    // inverter->SetAutoScan();

    HypoTestInverterResult* result = inverter->GetInterval();
    
    exclusion_limits.push_back(result->UpperLimit());
    exclusion_limits.push_back(result->GetExpectedUpperLimit(-2));
    exclusion_limits.push_back(result->GetExpectedUpperLimit(-1));
    exclusion_limits.push_back(result->GetExpectedUpperLimit(0));
    exclusion_limits.push_back(result->GetExpectedUpperLimit(1));
    exclusion_limits.push_back(result->GetExpectedUpperLimit(2));

    // Plot results
    TCanvas* c = new TCanvas(Form("CLs_mu95_S%d", int(point.m_s*100)), Form("CLs_mu95_S%d", int(point.m_s*100)), 800, 600);
    HypoTestInverterPlot* plot = new HypoTestInverterPlot("HTI_Result_Plot", Form("95\% Upper Limit scan on #mu for M_{S} = %.2f TeV;#mu;CL_{S}", point.m_s), result);
    plot->Draw("CLB 2CL");
    c->Draw();    
    output_file->cd();
    c->Write();      
    
    
    
    // Do some clean-up
    delete wspace;
    wspace = nullptr;
    delete sbModel;
    sbModel = nullptr;
    delete bModel;
    bModel = nullptr;
    delete toyData;
    toyData = nullptr;
    delete asympCalc;
    asympCalc = nullptr;
    delete inverter;
    inverter = nullptr;
    delete result;
    result = nullptr;
    delete plot;
    plot = nullptr;
    delete c;
    c = nullptr;

    return exclusion_limits;

}


/*
    Main analysis function
*/
void roostats_limits_run(const char* inputFile = nullptr) {
    
    gROOT->SetBatch(1);

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    vector<DataPoint> data = read_CSV(inputFile);
    vector<vector <double>> limits;

    TFile* output_file = TFile::Open("results/mChi1_5/roostats_results/out_D900/mu95_limits.root", "RECREATE");
    ofstream upper_file("results/mChi1_5/roostats_results/out_D900/upper_limits.csv");

    for(auto point : data) 
        limits.push_back(point_exclusion(point, output_file));

    upper_file << "M_S,obs_med,sig2_lo,sig1_lo,exp_med,sig1_hi,sig2_hi\n";
    
    // Print results
    for(int i = 0; i < limits.size(); i++) {
        upper_file << std::fixed << std::setprecision(2) << data[i].m_s << std::scientific;
        upper_file << Form(",%f,%f,%f,%f,%f,%f\n", limits[i][0], limits[i][1], limits[i][2], limits[i][3], limits[i][4], limits[i][5]);
    }

    upper_file.close();
    output_file->Write();
    output_file->Close();

}