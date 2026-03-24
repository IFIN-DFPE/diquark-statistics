/*
    Run the program with the path to the .csv as an argument
    E.g.: root 'roostats_analysis.cc("path/to/file/.csv")'
*/
#include "nlohmann/json.hpp"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooPoisson.h"
#include "RooBifurGauss.h"
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
#include "RooStats/HypoTestResult.h"
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


std::string path;


/*
    Data type holding:
    - the mass value
    - the theoretical signal yield and its uncertainty
    - the theoretical background yield and its uncertainty 
*/
struct DataPoint {
    double m_s; 
    double sig, sig_ml_uncrt;
    double hjj, hjj_ml_uncrt;
    double wj, wj_ml_uncrt;
    double qq2gg, qq2gg_ml_uncrt;
};

struct SignalUncertainties {
    double m_s;
    double lumi_uncrt;
    double JER_uncrt, JES_uncrt;
    double sig_PDF_uncrt, sig_scale_uncrt_hi, sig_scale_uncrt_lo;
    double hjj_PDF_uncrt, hjj_scale_uncrt_hi, hjj_scale_uncrt_lo;
    double wj_PDF_uncrt, wj_scale_uncrt_hi, wj_scale_uncrt_lo;
    double qq2gg_PDF_uncrt, qq2gg_scale_uncrt_hi, qq2gg_scale_uncrt_lo;
};

/*
    Function taking as input the path of a .csv datafile and 
    returns a vector of DataPoint type values.
*/
vector<DataPoint> read_CSV(std::string inputFile) {
    // Check if user has entered the path to the data file when running the macro
    if (inputFile.empty()) {
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
        if(getline(str, cell, ',')) point.sig_ml_uncrt= stod(cell);
        // Background values
        if(getline(str, cell, ',')) point.hjj = stod(cell);
        // Background uncertainties
        if(getline(str, cell, ',')) point.hjj_ml_uncrt = stod(cell);
        // Background values
        if(getline(str, cell, ',')) point.wj = stod(cell);
        // Background uncertainties
        if(getline(str, cell, ',')) point.wj_ml_uncrt = stod(cell);
        // Background values
        if(getline(str, cell, ',')) point.qq2gg = stod(cell);
        // Background uncertainties
        if(getline(str, cell, ',')) point.qq2gg_ml_uncrt = stod(cell);
            
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
         << "SIG\t\t\t";
    for (const DataPoint point : data)
        cout << fixed << setprecision(2) << point.m_s << "\t\t" << scientific 
             << point.sig << " ± " << point.sig_ml_uncrt << '\n';

    return data;
}


vector<SignalUncertainties> read_uncrt(std::string inputFile) {
// Check if user has entered the path to the data file when running the macro
    if (inputFile.empty()) {
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

    vector<SignalUncertainties> data;
    string line;

    // Skip header
    getline(csvFile, line);
    
    // Parse csv file by line
    while(getline(csvFile, line)) {
        stringstream str(line);
        string cell;
        SignalUncertainties point;

        // Mass values
        if(getline(str, cell, ',')) point.m_s = stod(cell);
        // PDF uncertainty
        if(getline(str, cell, ',')) point.sig_PDF_uncrt = stod(cell);
        // Scale uncertainty (high)
        if(getline(str, cell, ',')) point.sig_scale_uncrt_hi = stod(cell);
        // Scale uncertainty (low)
        if(getline(str, cell, ',')) point.sig_scale_uncrt_lo = stod(cell);
        // PDF uncertainty
        if(getline(str, cell, ',')) point.hjj_PDF_uncrt = stod(cell);
        // Scale uncertainty (high)
        if(getline(str, cell, ',')) point.hjj_scale_uncrt_hi = stod(cell);
        // Scale uncertainty (low)
        if(getline(str, cell, ',')) point.hjj_scale_uncrt_lo = stod(cell);
        // PDF uncertainty
        if(getline(str, cell, ',')) point.wj_PDF_uncrt = stod(cell);
        // Scale uncertainty (high)
        if(getline(str, cell, ',')) point.wj_scale_uncrt_hi = stod(cell);
        // Scale uncertainty (low)
        if(getline(str, cell, ',')) point.wj_scale_uncrt_lo = stod(cell);
        // PDF uncertainty
        if(getline(str, cell, ',')) point.qq2gg_PDF_uncrt = stod(cell);
        // Scale uncertainty (high)
        if(getline(str, cell, ',')) point.qq2gg_scale_uncrt_hi = stod(cell);
        // Scale uncertainty (low)
        if(getline(str, cell, ',')) point.qq2gg_scale_uncrt_lo = stod(cell);
        // JER uncertainty
        if(getline(str, cell, ',')) point.JER_uncrt = stod(cell);
        // JES uncertainty
        if(getline(str, cell, ',')) point.JES_uncrt = stod(cell);
        // Luminosity uncertainty
        if(getline(str, cell, ',')) point.lumi_uncrt = stod(cell);
            
        data.push_back(point);
    }

    csvFile.close();

    // Check if the program actually read something
    if(data.empty()) {
        cerr << "Error: The file contains no data.\n";
        exit(EXIT_FAILURE);
    }

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
vector<double> point_exclusion(DataPoint point, SignalUncertainties uncrt, TFile* output_file) {
    vector<double> exclusion_limits;

    RooRandom::randomGenerator()->SetSeed(69);

    // Create workspace
    RooWorkspace* wspace = new RooWorkspace("wspace");
    
    // Define a dummy variable for mass and create the shapes for sig and bkg
    RooRealVar mass("mass", "mass", point.m_s, 7., 10.);
    RooUniform sig_shape("sig_shape", "sig_shape", mass);
    RooUniform bkg_shape("bkg_shape", "bkg_shape", mass);


    /*
        !!!!!!!!!!!!!!!!!!!!!!!!
        !!!                  !!!
        !!!   Signal Model   !!!
        !!!                  !!!
        !!!!!!!!!!!!!!!!!!!!!!!!
    */


    // "Observed" yield from the simulation
    RooRealVar S0_obs("S0_obs", "Signal yield from simulation", point.sig, 1E-6, point.sig*5);
    // "True" yield from fitting
    RooRealVar S0_true("S0_true", "True signal yield", point.sig, 1E-6, point.sig*5);

    // ML uncertainty of the signal yield, using lognormal so we don't run into numerical issues
    RooRealVar sigmaML_sig("sigmaML", "Std dev of true signal yield", 1. + point.sig_ml_uncrt / point.sig, 1.0001, 100.);
    sigmaML_sig.setConstant();
    RooLognormal constraint_ML_sig("constraint_ML_sig", "constraint_ML_sig", S0_obs, S0_true, sigmaML_sig);

    // Luminosity uncertainty, both for signal and background 
    RooRealVar sigmaLumi("sigmaLumi", "std dev of lumi uncertainty", 1.0 + uncrt.lumi_uncrt/100, 1.0001, 100.);
    sigmaLumi.setConstant();
    RooRealVar theta_lumi("theta_lumi", "lumi uncertainty", 1., 1E-6, 5.);
    RooRealVar glob_lumi("glob_lumi", "global observable for lumi uncertainty", 1., 1E-6, 5.);
    glob_lumi.setConstant();
    RooLognormal constraint_lumi("constraint_lumi", "constraint_lumi", glob_lumi, theta_lumi, sigmaLumi);

    // JER uncertainty, both for signal and background
    RooRealVar sigmaJER("sigmaJER", "std dev of JER uncertainty", 1.0 + uncrt.JER_uncrt/100, 1.0001, 100.);
    sigmaJER.setConstant();
    RooRealVar theta_JER("theta_JER", "JER uncertainty", 1., 1E-6, 5.);
    RooRealVar glob_JER("glob_JER", "global observable for JER uncertainty", 1., 1E-6, 5.);
    glob_JER.setConstant();
    RooLognormal constraint_JER("constraint_JER", "constraint_JER", glob_JER, theta_JER, sigmaJER);

    // JES uncertainty, both for signal and background
    RooRealVar sigmaJES("sigmaJES", "std dev of JES uncertainty", 1.0 + uncrt.JES_uncrt/100, 1.0001, 100.);
    sigmaJES.setConstant();
    RooRealVar theta_JES("theta_JES", "JES uncertainty", 1., 1E-6, 5.);
    RooRealVar glob_JES("glob_JES", "global observable for JES uncertainty", 1., 1E-6, 5.);
    glob_JES.setConstant();
    RooLognormal constraint_JES("constraint_JES", "constraint_JES", glob_JES, theta_JES, sigmaJES);

    // PDF uncertainty for signal 
    RooRealVar sigmaPDF_sig("sigmaPDF_sig", "std dev of PDF uncertainty for signal", 1.0 + uncrt.sig_PDF_uncrt/100, 1.0001, 100.);
    sigmaPDF_sig.setConstant();
    RooRealVar theta_PDF_sig("theta_PDF_sig", "PDF uncertainty for signal", 1., 1E-6, 5.);
    RooRealVar glob_PDF_sig("glob_PDF_sig", "global observable for PDF uncertainty for signal", 1., 1E-6, 5.);
    glob_PDF_sig.setConstant();
    RooLognormal constraint_PDF_sig("constraint_PDF_sig", "constraint_PDF_sig", glob_PDF_sig, theta_PDF_sig, sigmaPDF_sig);

    /*
        Scale uncertainty for signal
        This is asymmetric, so we have to do some mumbo-jumbo to implement it.
        Easier to implement as a bifurcated Gaussian distribution, pulling on the parameters differently to the left and to the right
    */
    RooRealVar theta_scale_sig("theta_scale_sig", "scale uncertainty for signal", 1., 1E-6, 5.);
    RooRealVar sigma_scale_sig_hi("sigma_scale_sig_hi", "std dev of scale uncertainty for signal (upward)", uncrt.sig_scale_uncrt_hi/100, 1E-6, 5.);
    sigma_scale_sig_hi.setConstant();
    RooRealVar sigma_scale_sig_lo("sigma_scale_sig_lo", "std dev of scale uncertainty for signal (downward)", uncrt.sig_scale_uncrt_lo/100, 1E-6, 5.);
    sigma_scale_sig_lo.setConstant();
    RooRealVar mean_scale_sig("mean_scale_sig", "mean of scale uncertainty for signal", 1., 1E-6, 5.);
    mean_scale_sig.setConstant();
    RooBifurGauss constraint_scale_sig("constraint_scale_sig", "constraint_scale_sig", theta_scale_sig, mean_scale_sig, sigma_scale_sig_lo, sigma_scale_sig_hi);

    // Signal strength multiplier
    RooRealVar mu("mu", "signal multiplier", 1.0, 0.0, 20.);
    // Complete signal expression
    RooFormulaVar total_signal("total_signal", "@0*@1*@2*@3*@4*@5*@6", RooArgList(mu, S0_true, theta_lumi, theta_JER, theta_JES, theta_PDF_sig, theta_scale_sig));

    /*
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!                           !!!
        !!!   H+jj Background Model   !!!
        !!!                           !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    
    RooRealVar bkg_jjh_obs("bkg_jjh_obs", "Background yield from simulation for jjh", point.hjj, 1E-6, point.hjj*5);
    RooRealVar bkg_jjh_true("bkg_jjh_true", "True Background yield for jjh", point.hjj, 1E-6, point.hjj*5);
    
    RooRealVar sigmaML_bkg_jjh("sigmaML_bkg_jjh", "std dev of background yield for jjh", 1.0 + point.hjj_ml_uncrt / point.hjj, 1.0001, 100.);
    sigmaML_bkg_jjh.setConstant();
    RooLognormal constraint_ML_bkg_jjh("constraint_ML_bkg_jjh", "Yield constraint shape for jjh bkg", bkg_jjh_obs, bkg_jjh_true, sigmaML_bkg_jjh);

    RooRealVar sigmaPDF_bkg_jjh("sigmaPDF_bkg_jjh", "std dev of PDF uncertainty for jjh", 1.0 + uncrt.hjj_PDF_uncrt/100, 1.0001, 100.);
    sigmaPDF_bkg_jjh.setConstant();
    RooRealVar theta_PDF_bkg_jjh("theta_PDF_bkg_jjh", "PDF uncertainty for jjh", 1., 1E-6, 5.);
    RooRealVar glob_PDF_bkg_jjh("glob_PDF_bkg_jjh", "global observable for PDF uncertainty for jjh", 1., 1E-6, 5.);
    glob_PDF_bkg_jjh.setConstant();
    RooLognormal constraint_PDF_bkg_jjh("constraint_PDF_bkg_jjh", "constraint_PDF_bkg_jjh", glob_PDF_bkg_jjh, theta_PDF_bkg_jjh, sigmaPDF_bkg_jjh);

    RooRealVar theta_scale_bkg_jjh("theta_scale_bkg_jjh", "scale uncertainty for jjh", 1., 1E-6, 5.);
    RooRealVar sigma_scale_bkg_jjh_hi("sigma_scale_bkg_jjh_hi", "std dev of scale uncertainty for jjh (upward)", uncrt.hjj_scale_uncrt_hi/100, 1E-6, 5.);
    sigma_scale_bkg_jjh_hi.setConstant();
    RooRealVar sigma_scale_bkg_jjh_lo("sigma_scale_bkg_jjh_lo", "std dev of scale uncertainty for jjh (downward)", uncrt.hjj_scale_uncrt_lo/100, 1E-6, 5.);
    sigma_scale_bkg_jjh_lo.setConstant();
    RooRealVar mean_scale_bkg_jjh("mean_scale_bkg_jjh", "mean of scale uncertainty for jjh", 1., 1E-6, 5.);
    mean_scale_bkg_jjh.setConstant();
    RooBifurGauss constraint_scale_bkg_jjh("constraint_scale_bkg_jjh", "constraint_scale_bkg_jjh", theta_scale_bkg_jjh, mean_scale_bkg_jjh, sigma_scale_bkg_jjh_lo, sigma_scale_bkg_jjh_hi);

    RooFormulaVar total_bkg_jjh("total_bkg_jjh", "@0*@1*@2*@3*@4*@5", RooArgList(bkg_jjh_true, theta_lumi, theta_JER, theta_JES, theta_PDF_bkg_jjh, theta_scale_bkg_jjh));


    /*
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!                          !!!
        !!!   W+j Background Model   !!!
        !!!                          !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */

    RooRealVar bkg_wj_obs("bkg_wj_obs", "Background yield from simulation for wj", point.wj, 1E-6, point.wj*5);
    RooRealVar bkg_wj_true("bkg_wj_true", "True Background yield for wj", point.wj, 1E-6, point.wj*5);

    RooRealVar sigmaML_bkg_wj("sigmaML_bkg_wj", "std dev of background yield for wj", 1.0 + point.wj_ml_uncrt / point.wj, 1.0001, 100.);
    sigmaML_bkg_wj.setConstant();
    RooLognormal constraint_ML_bkg_wj("constraint_ML_bkg_wj", "Yield constraint shape for wj bkg", bkg_wj_obs, bkg_wj_true, sigmaML_bkg_wj);

    RooRealVar sigmaPDF_bkg_wj("sigmaPDF_bkg_wj", "std dev of PDF uncertainty for wj", 1.0 + uncrt.wj_PDF_uncrt/100, 1.0001, 100.);
    sigmaPDF_bkg_wj.setConstant();
    RooRealVar theta_PDF_bkg_wj("theta_PDF_bkg_wj", "PDF uncertainty for wj", 1., 1E-6, 5.);
    RooRealVar glob_PDF_bkg_wj("glob_PDF_bkg_wj", "global observable for PDF uncertainty for wj", 1., 1E-6, 5.);
    glob_PDF_bkg_wj.setConstant();
    RooLognormal constraint_PDF_bkg_wj("constraint_PDF_bkg_wj", "constraint_PDF_bkg_wj", glob_PDF_bkg_wj, theta_PDF_bkg_wj, sigmaPDF_bkg_wj);

    RooRealVar theta_scale_bkg_wj("theta_scale_bkg_wj", "scale uncertainty for wj", 1., 1E-6, 5.);
    RooRealVar sigma_scale_bkg_wj_hi("sigma_scale_bkg_wj_hi", "std dev of scale uncertainty for wj (upward)", uncrt.wj_scale_uncrt_hi/100, 1E-6, 5.);
    sigma_scale_bkg_wj_hi.setConstant();
    RooRealVar sigma_scale_bkg_wj_lo("sigma_scale_bkg_wj_lo", "std dev of scale uncertainty for wj (downward)", uncrt.wj_scale_uncrt_lo/100, 1E-6, 5.);
    sigma_scale_bkg_wj_lo.setConstant();
    RooRealVar mean_scale_bkg_wj("mean_scale_bkg_wj", "mean of scale uncertainty for wj", 1., 1E-6, 5.);
    mean_scale_bkg_wj.setConstant();
    RooBifurGauss constraint_scale_bkg_wj("constraint_scale_bkg_wj", "constraint_scale_bkg_wj", theta_scale_bkg_wj, mean_scale_bkg_wj, sigma_scale_bkg_wj_lo, sigma_scale_bkg_wj_hi);

    RooFormulaVar total_bkg_wj("total_bkg_wj", "@0*@1*@2*@3*@4*@5", RooArgList(bkg_wj_true, theta_lumi, theta_JER, theta_JES, theta_PDF_bkg_wj, theta_scale_bkg_wj));


    /*
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!                             !!!
        !!!   qq->gg Background Model   !!!
        !!!                             !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */

    RooRealVar bkg_qq2gg_obs("bkg_qq2gg_obs", "Background yield from simulation for qq->gg", point.qq2gg, 1E-6, point.qq2gg*5);
    RooRealVar bkg_qq2gg_true("bkg_qq2gg_true", "True Background yield for qq->gg", point.qq2gg, 1E-6, point.qq2gg*5);

    RooRealVar sigmaML_bkg_qq2gg("sigmaML_bkg_qq2gg", "std dev of background yield for qq->gg", 1.0 + point.qq2gg_ml_uncrt / point.qq2gg, 1.0001, 100.);
    sigmaML_bkg_qq2gg.setConstant();
    RooLognormal constraint_ML_bkg_qq2gg("constraint_ML_bkg_qq2gg", "Yield constraint shape for qq->gg bkg", bkg_qq2gg_obs, bkg_qq2gg_true, sigmaML_bkg_qq2gg);

    RooRealVar sigmaPDF_bkg_qq2gg("sigmaPDF_bkg_qq2gg", "std dev of PDF uncertainty for qq->gg", 1.0 + uncrt.qq2gg_PDF_uncrt/100, 1.0001, 100.);
    sigmaPDF_bkg_qq2gg.setConstant();
    RooRealVar theta_PDF_bkg_qq2gg("theta_PDF_bkg_qq2gg", "PDF uncertainty for qq->gg", 1., 1E-6, 5.);
    RooRealVar glob_PDF_bkg_qq2gg("glob_PDF_bkg_qq2gg", "global observable for PDF uncertainty for qq->gg", 1., 1E-6, 5.);
    glob_PDF_bkg_qq2gg.setConstant();
    RooLognormal constraint_PDF_bkg_qq2gg("constraint_PDF_bkg_qq2gg", "constraint_PDF_bkg_qq2gg", glob_PDF_bkg_qq2gg, theta_PDF_bkg_qq2gg, sigmaPDF_bkg_qq2gg);

    RooRealVar theta_scale_bkg_qq2gg("theta_scale_bkg_qq2gg", "scale uncertainty for qq->gg", 1., 1E-6, 5.);
    RooRealVar sigma_scale_bkg_qq2gg_hi("sigma_scale_bkg_qq2gg_hi", "std dev of scale uncertainty for qq->gg (upward)", uncrt.qq2gg_scale_uncrt_hi/100, 1E-6, 5.);
    sigma_scale_bkg_qq2gg_hi.setConstant();
    RooRealVar sigma_scale_bkg_qq2gg_lo("sigma_scale_bkg_qq2gg_lo", "std dev of scale uncertainty for qq->gg (downward)", uncrt.qq2gg_scale_uncrt_lo/100, 1E-6, 5.);
    sigma_scale_bkg_qq2gg_lo.setConstant();
    RooRealVar mean_scale_bkg_qq2gg("mean_scale_bkg_qq2gg", "mean of scale uncertainty for qq->gg", 1., 1E-6, 5.);
    mean_scale_bkg_qq2gg.setConstant();
    RooBifurGauss constraint_scale_bkg_qq2gg("constraint_scale_bkg_qq2gg", "constraint_scale_bkg_qq2gg", theta_scale_bkg_qq2gg, mean_scale_bkg_qq2gg, sigma_scale_bkg_qq2gg_lo, sigma_scale_bkg_qq2gg_hi);

    RooFormulaVar total_bkg_qq2gg("total_bkg_qq2gg", "@0*@1*@2*@3*@4*@5", RooArgList(bkg_qq2gg_true, theta_lumi, theta_JER, theta_JES, theta_PDF_bkg_qq2gg, theta_scale_bkg_qq2gg));



    // Complete background expression
    RooRealVar bkg_obs("bkg_obs", "Total background yield from simulation", point.hjj + point.wj + point.qq2gg, 1E-6, (point.hjj + point.wj + point.qq2gg)*5);
    RooFormulaVar total_bkg("total_bkg", "@0+@1+@2", RooArgList(total_bkg_jjh, total_bkg_wj, total_bkg_qq2gg));


    // Extend the pdfs over the entire mass region
    RooExtendPdf ext_sig("ext_sig", "Extended Signal PDF", sig_shape, total_signal);
    RooExtendPdf ext_bkg("ext_bkg", "Extedned Bkg PDF", bkg_shape, total_bkg);

    // Add up pdfs into the full model
    RooAddPdf sb_tmp("sb_tmp", "sb_tmp", RooArgSet(ext_sig, ext_bkg));
    RooProdPdf sb_full("sb_full", "sb_full", RooArgSet(sb_tmp, constraint_ML_sig, constraint_ML_bkg_jjh, constraint_ML_bkg_wj, constraint_ML_bkg_qq2gg, 
                                                        constraint_lumi, constraint_JER, constraint_JES, 
                                                        constraint_PDF_sig, constraint_scale_sig, 
                                                        constraint_PDF_bkg_jjh, constraint_scale_bkg_jjh,
                                                        constraint_PDF_bkg_wj, constraint_scale_bkg_wj,
                                                        constraint_PDF_bkg_qq2gg, constraint_scale_bkg_qq2gg));
    

    // Create s+b model configuration
    ModelConfig* sbModel = new ModelConfig("sbModel", wspace);
    // Add parameters to the model configuration
    sbModel->SetPdf(sb_full);
    sbModel->SetParametersOfInterest(mu);
    sbModel->SetNuisanceParameters({S0_true, bkg_jjh_true, bkg_wj_true, bkg_qq2gg_true, 
                                    theta_lumi, theta_JER, theta_JES, 
                                    theta_PDF_sig, theta_scale_sig, 
                                    theta_PDF_bkg_jjh, theta_scale_bkg_jjh,
                                    theta_PDF_bkg_wj, theta_scale_bkg_wj,
                                    theta_PDF_bkg_qq2gg, theta_scale_bkg_qq2gg});
    sbModel->SetObservables(mass);
    sbModel->SetGlobalObservables({S0_obs, bkg_jjh_obs, bkg_wj_obs, bkg_qq2gg_obs,
                                    glob_lumi, glob_JER, glob_JES, 
                                    glob_PDF_sig, mean_scale_sig, 
                                    glob_PDF_bkg_jjh, mean_scale_bkg_jjh,
                                    glob_PDF_bkg_wj, mean_scale_bkg_wj,
                                    glob_PDF_bkg_qq2gg, mean_scale_bkg_qq2gg});
    sbModel->SetSnapshot(mu);


    // Create b-only model configuration 
    ModelConfig* bModel = new ModelConfig("bModel", wspace);
    bModel->SetPdf(sb_full);
    RooRealVar bPoi = RooRealVar("bPoi", "signal multiplier in b-only model", 0.);
    bPoi.setConstant();
    bModel->SetParametersOfInterest(bPoi);
    bModel->SetNuisanceParameters({S0_true, bkg_jjh_true, bkg_wj_true, bkg_qq2gg_true, 
                                    theta_lumi, theta_JER, theta_JES, 
                                    theta_PDF_sig, theta_scale_sig, 
                                    theta_PDF_bkg_jjh, theta_scale_bkg_jjh,
                                    theta_PDF_bkg_wj, theta_scale_bkg_wj,
                                    theta_PDF_bkg_qq2gg, theta_scale_bkg_qq2gg});
    bModel->SetObservables(mass);
    bModel->SetGlobalObservables({S0_obs, bkg_jjh_obs, bkg_wj_obs, bkg_qq2gg_obs,
                                    glob_lumi, glob_JER, glob_JES, 
                                    glob_PDF_sig, mean_scale_sig, 
                                    glob_PDF_bkg_jjh, mean_scale_bkg_jjh,
                                    glob_PDF_bkg_wj, mean_scale_bkg_wj,
                                    glob_PDF_bkg_qq2gg, mean_scale_bkg_qq2gg});
    bModel->SetSnapshot(bPoi);

    // Create dataset       
    RooDataSet* toyData = sb_full.generate(RooArgSet(S0_obs, bkg_jjh_obs, bkg_wj_obs, bkg_qq2gg_obs), 1);
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
    
    exclusion_limits.push_back(result->UpperLimit()*point.sig);
    exclusion_limits.push_back(result->GetExpectedUpperLimit(-2)*point.sig);
    exclusion_limits.push_back(result->GetExpectedUpperLimit(-1)*point.sig);
    exclusion_limits.push_back(result->GetExpectedUpperLimit(0)*point.sig);
    exclusion_limits.push_back(result->GetExpectedUpperLimit(1)*point.sig);
    exclusion_limits.push_back(result->GetExpectedUpperLimit(2)*point.sig);

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
void roostats_limits_run(const char* configFilePath = nullptr) {
    
    gROOT->SetBatch(1);

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    
    // Load the file containing the paths
    std::ifstream pathFile("analysis_paths.json");
    nlohmann::json paths = nlohmann::json::parse(pathFile);

    // Load configuration file
    std::ifstream configFile(configFilePath);
    nlohmann::json config = nlohmann::json::parse(configFile);

    auto discriminator = int(config["discriminator"].get<float>()*1000);
    auto process = config["process"].get<std::string>();
    path = paths[process].get<std::string>();
    std::string inputFilePath = path + Form("/signal_yields/sig_bkg_D%d.csv", discriminator);
    std::string uncrtFilePath = path + Form("/signal_yields/uncrt.csv");

    // Read signal yields and uncertainties 
    vector<DataPoint> data = read_CSV(inputFilePath);
    vector<SignalUncertainties> uncrt_data = read_uncrt(uncrtFilePath);
    vector<vector <double>> limits;

    std::string out_path = path + Form("/roostats_results/out_D%d/mu95_limits.root", discriminator);
    TFile* output_file = TFile::Open(out_path.c_str(), "RECREATE");
    ofstream upper_file(path + Form("/roostats_results/out_D%d/upper_limits.csv", discriminator));

    for(int i = 0; i < data.size(); i++) {
        limits.push_back(point_exclusion(data[i], uncrt_data[i], output_file));
    }

    upper_file << "M_S,obs_med,sig2_lo,sig1_lo,exp_med,sig1_hi,sig2_hi\n";
    
    // Print results
    for(int i = 0; i < limits.size(); i++) {
        upper_file << std::fixed << std::setprecision(2) << data[i].m_s << std::scientific;
        upper_file << Form(",%f,%f,%f,%f,%f,%f\n", limits[i][0], limits[i][1], limits[i][2], limits[i][3], limits[i][4], limits[i][5]);
    }

    pathFile.close();
    configFile.close();
    upper_file.close();
    output_file->Write();
    output_file->Close();

}