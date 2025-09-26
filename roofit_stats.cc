#include <chrono>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>
#include <omp.h>


#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooPoisson.h"
#include "RooDataSet.h"
#include "RooUniform.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooNumIntConfig.h"
#include "RooProdPdf.h"
#include "RooRealConstant.h"
#include "RooExtendPdf.h"
#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "TApplication.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"
#include "TROOT.h"



/*
======================== RUNNING FROM SCRIPT ================================
    It is enough to run the command:
    ./run_roofit.sh
    and the program will be compiled and run with the parameters given in stats.ini
*/



/*
======================== COMPILATION INSTRUCTIONS ================================
    g++ -O3 `root-config --cflags --libs` -l RooStats -l RooFitCore -l RooFit -o roofit_stats roofit_stats.cc
======================= or if it doesn't work, try: ==============================
    g++ -O3 roofit_stats.cc `root-config --cflags` -o roofit_stats -L$ROOTSYS/lib -lRooFitCore -lRooFit -lRooStats `root-config --glibs`
*/



/*
========================== RUNNING INSTRUCTIONS ==================================
    ./roofit_stats.cc stats.ini
    stats.ini contains the following initialization parameters
        - discriminator: the discriminator D used in machine learning
        - runType: can take two values depending on the working mode desired
            - full: runs the analysis over all masses
            - point: runs the analysis on a single mass point
        - mass: if runType = point, specify the mass of Suu for the analysis
        - nPseudoExps: the number of pseudo-experiments to be run
        - nToys: the number of toys generated for each pseudo-experiment
        
*/ 



/*
    In what follows the terminology is:
        - PSEUDOEXPERIMENT = Monte Carlo (MC) sampled under the bkg-only hypothesis WITH DEFAULT YIELD.
                            We treat this "first layer" as real data, in the absence of signal,
                            i.e. what one would do with experimental data once one has it.
                            These must be performed many times for each parameter-space point in order
                            to count how often such a point is excluded at a certain confidence level.
        - TOY = MC sampled under either S+B or B-only hypothesis with YIELD FIXED TO THE ONE OBTAINED FROM
                FITTING the pseudoexperiment (there can be variations here). The point is that "toys" are
                always sampled, even when one has real data, to construct the distribution of the test statistic
                under the two competing hypotheses. Each pseudo-experiment is then excluded or not via
                CLS using these distributions.

    The test statistic is q_mu as defined here:
        - https://indico.cern.ch/event/126652/contributions/1343592/attachments/80222/115004/Frequentist_Limit_Recommendation.pdf
        - https://cds.cern.ch/record/1379837
    but with mu = 1. N.B. mu_hat is still there, just normal "mu" is set to 1.
*/



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


struct HistResult {
    TH1D hCLS;
    TH1D hCLSB;
    TH1D hCLB;
    double nTimesExcluded;
    double nTotalSB;
};


// Value of the discriminator used in ML
int discriminator;
// Number of PSEUDOEXPERIMENTS. Recommend 100 at least, maybe 1000
int nPseudoExps = 100;
// Number of TOYS per PSEUDOEXPERIMENTS. Recommend 1000 at least, maybe 10000
int nToys = 1000;



/*
    Function taking as argument the path to the config file, parsing the 
    config file and returning a map containing:
        - key: the defined property
        - value: the value of said property
*/
std::map<std::string, std::string> load_Config(const char* inputFile) {
    // Check if user has entered the path to the config file
    if(!inputFile) {
        throw std::runtime_error("Error: Please enter the name of the configuration file.\n");
    }

    std::ifstream configFile(inputFile);
    std::map<std::string, std::string> conf_Map;

    // Parse config file
    std::string line;
    while(getline(configFile, line)) {
        // Trim white spaces
        auto trim = [](std::string& s) {
            size_t start = s.find_first_not_of(" \t");
            size_t end = s.find_last_not_of(" \t");
            if(start == std::string::npos) {s.clear(); return;}
            s = s.substr(start, end-start+1);
        };
        
        // Check if line is empty
        trim(line);
        if(line.empty()) continue;

        // Find the position of the '=' sign
        size_t limPos = line.find('=');
        if(limPos == std::string::npos) continue;

        // Read the key and value from the line and trim white spaces
        std::string key = line.substr(0, limPos);
        trim(key);
        std::string value = line.substr(limPos+1);
        trim(value);

        conf_Map[key] = value;

    }
    configFile.close();

    if(conf_Map.empty()) {
        throw std::runtime_error("Error: Configuration file is empty, please check it.\n");
    }

    return conf_Map;

}



/*
    Function taking as argument the path of a .csv datafile and 
    returns a vector of DataPoint type values.
*/
std::vector<DataPoint> read_CSV(const char* inputFile) {
    // Check if user has entered the path to the data file 
    if(!inputFile) {
        throw std::runtime_error("Error: Please enter the name of the data file to be read.\n");
    }

    // Check if the file can be opened or not
    std::ifstream csvFile(inputFile);
    if (!csvFile.is_open()) {
        throw std::runtime_error("Error: Please check the path of the input file.\n");
    }

    std::cout << "Reading data from: " << inputFile << '\n';  

    std::vector<DataPoint> data;
    std::string line;

    // Skip header
    getline(csvFile, line);
    
    // Parse csv file by line
    while(getline(csvFile, line)) {
        std::stringstream str(line);
        std::string cell;
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
        throw std::runtime_error("Error: The data file is empty.\n");
    }

    // Print the read data for debugging
    std::cout << "\n=== Summary by Mass Points ===\n";
    std::cout << "M_s [TeV]\t"
         << "SIG\t\t\t" 
         << "BKG\n";
    for (const DataPoint point : data)
        std::cout << std::fixed << std::setprecision(2) << point.m_s << "\t\t" << std::scientific 
             << point.sig << " ± " << point.sigma_sig << "\t"
             << point.bkg << " ± " << point.sigma_bkg << "\n";
    std::cout << "\n" << std::fixed;

    return data;
}



void generate(DataPoint point, int num_point) {
    RooWorkspace* wspace;

    #pragma omp critical
    {
        wspace = new RooWorkspace(Form("wspace_%d", num_point));
    }

    #pragma omp critical
    {
        RooRealVar mass(Form("mass_%d", num_point), Form("mass_%d", num_point), point.m_s, 7., 10.);
        RooUniform sig_shape(Form("sig_shape_%d", num_point), Form("sig_shape_%d", num_point), mass);
        RooUniform bkg_shape(Form("bkg_shape_%d", num_point), Form("bkg_shape_%d", num_point), mass);

        wspace->import(sig_shape);
        wspace->import(bkg_shape);
    }

    #pragma omp critical
    {
        wspace->Print("v");
    }


    delete wspace;
    wspace = nullptr;

}



HistResult runPseudoExp(const RooArgSet& sb_params_default_vals,
                const RooArgSet& in_globals, const RooArgSet& sb_global_default_vals,
                RooRealVar mu, RooRealVar mass, const RooProdPdf& in_sb_full,
                const std::chrono::_V2::system_clock::time_point startFull,
                int num_thread) {

    const auto startPseudoExp = std::chrono::high_resolution_clock::now();
    
    RooProdPdf* sb_full;
    RooArgSet* globals;
    RooArgSet sb_params;
                    
    #pragma omp critical 
    {
    // OPTION 1. Set the XX_true parameters to the default values in order to generate PSEUDOEXPERIMENTS the same way every-time.
    
    sb_params.assign(sb_params_default_vals);

    sb_full = (RooProdPdf*)in_sb_full.cloneTree(Form("sb_full_%d", num_thread));
    // sb_full->Print("v");

    // 1. Global observables
    // We regenerate those because if we were to "repeat" the experiment, you would redo the simulations as well, obtaining different
    // S0_obs and B0_obs.
    globals = (RooArgSet*)in_globals.snapshot();
    RooDataSet* ds_global = sb_full->generate(*globals, 1);

    globals->assign(*ds_global->get(0));

    }

    // 2. Main data
    mu.setVal(0.0); // we're generating bkg-only PSEUDOEXPERIMENT, so turn off the signal
    mu.setConstant(true);
    RooDataSet *ds = sb_full->generate({mass}, RooFit::Extended());
    ds->setGlobalObservables(*globals); // set global observables so that RooFit does not fit S0_obs and B0_obs


    // fit to s+b and record parameters and nll
    // set and fix mu to 1.0, meaning full signal strength
    sb_params.assign(sb_params_default_vals);
    mu.setVal(1.0);
    mu.setConstant(true);
    
    RooFitResult *result_mu;
    #pragma omp critical 
    {
    result_mu = sb_full->fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                            RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                            RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
    }
                                            // std::cout << "================ mu=1 fit ==============================" << std::endl;
    // result_mu->Print("V");
    Double_t nll_mu = result_mu->minNll();
    RooArgSet *params_fit_sb = sb_params.snapshot();
    // params_fit_sb->Print("V");
    

    // fit to b-only and record parameters and nll
    // allow mu to float as described in the documents.
    mu.setVal(1E-5);
    mu.setConstant(false);
    RooFitResult *result_mu_hat = sb_full->fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
    // std::cout << "================ mu_hat fit ==============================" << std::endl;
    // result_mu_hat->Print("V");
    Double_t nll_mu_hat = result_mu_hat->minNll();
    // construct the test statistic of the "observed" data
    // Q = -2 ln (L(mu=1)/L(mu_hat)) = -2 * (ln(L(mu=1)) - ln(L(mu_hat))) = 2*(nll(mu=1) - nll(mu_hat))
    Double_t q_obs = 2. * (nll_mu - nll_mu_hat);
    

    // fit under pure background to find nuissance parameters
    mu.setVal(0.0);
    mu.setConstant(true);
    RooFitResult *result_b = sb_full->fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                            RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                            RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
    // std::cout << "================ bkg fit ==============================" << std::endl;
    // result_b->Print("V");
    RooArgSet *params_fit_b = sb_params.snapshot();
    // params_fit_b->Print("V");


    // Build distribution of q under b-only and count how many times the TOY has a higher q than the PSEUDOEXPERIMENT.
    // This number is 1-CLB (as defined in https://cds.cern.ch/record/1379837/files/NOTE2011_005.pdf)
    Double_t n_higher_bkg = 0.0;
    // histogram for the distribution of q under B-only
    // This will show that when using the "unconstrained" ensemble, the q_mu dsitributions do not depend on the pseudoexperiment
    // So one can derive them only once. In this code we do not use this property, but it's very nice to have in mind.
    for (Int_t i_toy = 0; i_toy < nToys; ++i_toy)
    {
        const auto startToy = std::chrono::high_resolution_clock::now();
        // OPTION 1. Set the bkg yield to the one obtained from the PSEUDOEXPERIMENT fit.
        // sb_params.assign(*params_fit_b);
        sb_params.assign(*params_fit_b);
        

        // generate the "XX_obs" variables
        RooDataSet *ds_global_toy = sb_full->generate(*globals, 1);
        globals->assign(*ds_global_toy->get(0));
        
        mu.setVal(0.0);
        RooDataSet *ds_toy = sb_full->generate({mass}, RooFit::Extended());
        ds_toy->setGlobalObservables(*globals);
        mu.setVal(1.0);
        mu.setConstant(true);
        RooFitResult *result_mu_toy = sb_full->fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                    RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                    RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
        Double_t nll_mu_toy = result_mu_toy->minNll();
        
        // mu.setVal(1E-5);
        mu.setConstant(false);
        RooFitResult *result_mu_hat_toy = sb_full->fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                        RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                        RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
        Double_t nll_mu_hat_toy = result_mu_hat_toy->minNll();

        Double_t q_toy = 2. * (nll_mu_toy - nll_mu_hat_toy);
        delete result_mu_toy;
        result_mu_toy = nullptr;
        delete result_mu_hat_toy;
        result_mu_hat_toy=nullptr;
        delete ds_toy;
        ds_toy = nullptr;
        delete ds_global_toy;
        ds_global_toy = nullptr;
        
        
        if (q_toy >= q_obs)
        {
            n_higher_bkg += 1.0;
        }


        // remove the line below if you want to do debug printing.
        const auto stopToy = std::chrono::high_resolution_clock::now();
        auto durationToy = std::chrono::duration_cast<std::chrono::microseconds>(stopToy - startToy);
        auto durationPE = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startPseudoExp);
        auto durationNow = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startFull);

        Double_t durationToyNice = durationToy.count() / 1000.; // ms
        Double_t durationPENice = durationNow.count() / 1000.;  // s
        Double_t durationNowNice = durationNow.count() / 1000.; // s
        if ((i_toy+1) % (nToys / 2) == 0)
        {
            std::cout << " Duration: " << durationNowNice << "[s]\n" << std::flush;
        }
    }

    // Build distribution of q under S+B and count how many times the TOY has a higher q than the PSEUDOEXPERIMENT.
    // This number is CLSB (as defined in https://cds.cern.ch/record/1379837/files/NOTE2011_005.pdf)
    Double_t n_higher_sb = 0.0;
    for (Int_t i_toy = 0; i_toy < nToys; ++i_toy)
    {
        const auto startToy = std::chrono::high_resolution_clock::now();
        // OPTION 1. set the signal and bkg yield to the one obtained from the PSEUDOEXPERIMENT fit.
        sb_params.assign(sb_params_default_vals);

        // sb_full->generate(b_params);
        // OPTION 2. Sample the signal and bkg yield before generating the PSEUDOEXPERIMENTS.
        // NOT a frequentist construction, but we should not fear this.
        RooDataSet *ds_global_toy = sb_full->generate(*globals, 1);
        globals->assign(*ds_global_toy->get(0));
        mu.setVal(1.0);
        RooDataSet *ds_toy = sb_full->generate({mass}, RooFit::Extended());
        ds_toy->setGlobalObservables(*globals);
        mu.setVal(1.0);
        mu.setConstant(true);
        RooFitResult *result_mu_toy = sb_full->fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                    RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                    RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
        Double_t nll_mu_toy = result_mu_toy->minNll();
        mu.setVal(1E-5);
        mu.setConstant(false);
        RooFitResult *result_mu_hat_toy = sb_full->fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                        RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                        RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(*globals));
        Double_t nll_mu_hat_toy = result_mu_hat_toy->minNll();

        Double_t q_toy = 2. * (nll_mu_toy - nll_mu_hat_toy);

        delete result_mu_toy;
        result_mu_toy = nullptr;
        delete result_mu_hat_toy;
        result_mu_hat_toy=nullptr;
        delete ds_toy;
        ds_toy = nullptr;
        delete ds_global_toy;
        ds_global_toy = nullptr;
        

        if (q_toy >= q_obs)
        {
            n_higher_sb += 1.0;
        }

        // remove the line below if you want to do debug printing.
        const auto stopToy = std::chrono::high_resolution_clock::now();
        auto durationToy = std::chrono::duration_cast<std::chrono::microseconds>(stopToy - startToy);
        auto durationPE = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startPseudoExp);
        auto durationNow = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startFull);

        Double_t durationToyNice = durationToy.count() / 1000.; // ms
        Double_t durationPENice = durationNow.count() / 1000.;  // s
        Double_t durationNowNice = durationNow.count() / 1000.; // s
        if ((i_toy+1) % (nToys / 2) == 0)
        {
            std::cout << " Duration: " << durationNowNice << "[s]\n" << std::flush;
        }
    }

    // Computation of CLs. Exactly as described in the notes
    Double_t clSB = (1.0 * n_higher_sb) / (1.0 * nToys);
    Double_t clB = (1.0 * n_higher_bkg) / (1.0 * nToys);
    Double_t CLS = clB > 0. ? clSB / clB : 0.0;
    HistResult result;
    result.hCLS.Fill(CLS);
    result.hCLSB.Fill(clSB);
    result.hCLB.Fill(clB);
    // Assume 95% confindence level
    if (CLS < 0.05)
    {
        result.nTimesExcluded = 1.0;
    }
    result.nTotalSB = n_higher_sb;
    
    delete globals;


    return result;

}



/*
    Function that performs the statistical analysis and checks if the points are excluded or not
    As an argument, pass the DataPoint associated to a given mass of Suu
*/
void analysisRun(DataPoint point, std::ofstream &prob_file) {

    std::cout << Form("=== Running analysis for M_S = %.2f TeV ===\n", point.m_s);

    /*
        Define a dummy "discriminator" variable, for example reconstructed mass.
        We will assume it is uniformly distributed for both signal and background.
        In this case it will simply not matter that we have it, but the construction
        is easier to udnerstand when all ingredients are there.
    */

    RooRealVar mass("mass", "mass", point.m_s, 7., 10.);
    RooUniform sig_shape("sig_shape", "sig_shape", mass);
    RooUniform bkg_shape("bkg_shape", "bkg_shape", mass);

    // define constraints on yields as lognormal distribution (to not allow negative yields).
    // Also plot the p.d.f. of the yields.
    
    // ================ signal yield ===============================
    // "Observed" yield, i.e. the value you obtained from simulation
    RooRealVar S0_obs("S0_obs", "Signal yield from simulation", point.sig, 1E-6, 100.);
    // "True" yield, i.e. the nuisance parameter we will obtain from the fit. This one enters the poisson rate
    RooRealVar S0_true("S0_true", "True signal yield", point.sig, 1E-6, 100.);
    // The sigma of the yield. This you obtained from simulation. We will treat is as a constant.
    // N.B.: I'm using the lognormal distribution so the sigma is a scale variable, larger than one.
    //       A value of 1.1 means "S0_obs will be within 10% of S0_true most of the time".
    RooRealVar sigmaS("sigmaS", "Std dev of true signal yield", 1. + point.sigma_sig / point.sig, 1.0001, 100.);
    sigmaS.setConstant();
    // Define shape
        RooLognormal constraint_yield_sig("constraint_yield_sig", "constraint_yield_sig", S0_obs, S0_true, sigmaS);

    // The famous scale "mu". mu=1 means full effect of signal. mu=0 means no signal.
    RooRealVar mu("mu", "signal multiplier", 1.0, 0.0, 1.0);
    mu.setConstant(true); // for the moment
    // define "full" signal yield
    RooFormulaVar total_signal("total_signal", "mu*S0_true", RooArgList(mu, S0_true));

    // Plot the shape of the signal PDF
    RooPlot *frameSigConstraint = S0_obs.frame(0., 30., 50);
    frameSigConstraint->SetTitle("Signal yield constraint shape");
    TCanvas *cSigConstraint = new TCanvas("cSigConstraint", "cSigConstraint");
    constraint_yield_sig.plotOn(frameSigConstraint);
    frameSigConstraint->Draw();
    cSigConstraint->Update();

    // ================ background yield ===============================
    // "Observed" background yield. This you obtained from simulation.
    RooRealVar B0_obs("B0_obs", "Background yield from simulation", point.bkg, 1E-6, 1000.);
    // "True" yield, i.e. the nuisance parameter we will obtain from the fit. This one enters the poisson rate
    RooRealVar B0_true("B0_true", "True Background yield", point.bkg, 1E-6, 1000.);
    // The sigma of the yield. This you obtained from simulation. We will treat is as a constant.
    // See note above for details.
    RooRealVar sigmaB("sigmaB", "std dev of background yield", 1.0 + point.sigma_bkg / point.bkg, 1.01, 100.);
    sigmaB.setConstant();
    // shape the constraint for the background
    RooLognormal constraint_yield_bkg("constraint_yield_bkg", "Yield constraint shape for bkg", B0_obs, B0_true, sigmaB);

    // Plot the shape of the background PDF
    RooPlot *frameBkgConstraint = B0_obs.frame(0., 20., 50);
    frameBkgConstraint->SetTitle("Bkg yield constraint shape");
    TCanvas *cBkgConstraint = new TCanvas("cBkgConstraint", "cBkgConstraint");
    constraint_yield_bkg.plotOn(frameBkgConstraint);
    frameBkgConstraint->Draw();
    cBkgConstraint->Update();

    // define the extended pdfs (i.e. shapes x poisson term)
    // Notice they are extended on the XX_true variable
    RooExtendPdf ext_sig("ext_sig", "Extended Signal PDF", sig_shape, total_signal);
    RooExtendPdf ext_bkg("ext_bkg", "Extedned Bkg PDF", bkg_shape, B0_true);

    // define full pdfs (i.e. extended x constraint terms)
    RooAddPdf sb_tmp("sb_tmp", "sb_tmp", RooArgList(ext_sig, ext_bkg));
    RooProdPdf sb_full("sb_full", "sb_full", RooArgSet(sb_tmp, constraint_yield_sig, constraint_yield_bkg));

    // Default values of fit parameters.
    // these are used to reset the fit parameters to default values before
    // generating the PSEUDOEXPERIMENTS and TOYS.
    RooArgSet sb_params;
    sb_params.add(mu);
    sb_params.add(S0_true);
    sb_params.add(B0_true);
    RooArgSet *sb_params_default_vals = sb_params.snapshot();

    // Global observables
    // We will sample S0_obs and B0_obs in each PSEUDOEXPERIMENT and in each TOY, before generating
    // the respective PSEUDOEXPERIMENT/TOY.
    RooArgSet globals;
    globals.add(S0_obs);
    globals.add(B0_obs);
    RooArgSet *sb_global_default_vals = globals.snapshot();

    if(std::filesystem::create_directories(Form("results/mChi2/roofit_results/out_D%d", discriminator)))
    ;
    TFile *output_file = TFile::Open(Form("results/mChi2/roofit_results/out_D%d/output_S%d.root", discriminator, int(point.m_s*100)), "RECREATE");    
    // dsitributions of CLS, CLSB and CLB
    TH1D *hCLS = new TH1D("hCLS", "CLS", 100, 0., 1.);
    TH1D *hCLSB = new TH1D("hCLSB", "CLSB", 100, 0., 1.);
    TH1D *hCLB = new TH1D("hCLB", "1-CLB", 100, 0., 1.);

    // Number of PSEUDOEXPERIMENTS. Recommend 100 at least, maybe 1000
    // Number of TOYS per PSEUDOEXPERIMENTS. Recommend 1000 at least, maybe 10000
    TH1D *hTimeBkgToy = new TH1D("hTimeBkgToy", "Bkg toys; Toy time [ms]", 100, 0., 20.);
    TH1D *hTimeSBToy = new TH1D("hTimeSBToy", "SB toys; Toy time [ms]", 100, 0., 20.);
    TH1D *hTimePseudoExp = new TH1D("hTimePseudoExp", "Pseudo-experiments; PE time [s]", 100, 0., 20. * nToys / 1000.);

    // Counter to check how often we exclude the parameter point (in which the signal yield is what it is above).
    Double_t nTimesExcluded = 0.0;
    Double_t nTotalSB = 0.;
    const auto startFull = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for 
    for(int i = 0; i < nPseudoExps; i++) {
        HistResult result = runPseudoExp(*sb_params_default_vals, globals, *sb_global_default_vals,
                   mu, mass, sb_full, startFull, i);
        
        // generate(point, i);
        
        #pragma omp critical 
        {
            std::cout << i << '\n';
        }
    }

    /*

    for (Int_t i = 0; i < nPseudoExps; ++i) {
        const auto startPseudoExp = std::chrono::high_resolution_clock::now();
        // OPTION 1. Set the XX_true parameters to the default values in order to generate PSEUDOEXPERIMENTS the same way every-time.
        sb_params.assign(*sb_params_default_vals);

        // ############### GENERATING REAL DATA ####################

        // 1. Global observables
        // We regenerate those because if we were to "repeat" the experiment, you would redo the simulations as well, obtaining different
        // S0_obs and B0_obs.
        RooDataSet *ds_global = sb_full.generate(globals, 1);
        globals.assign(*ds_global->get(0));
        // std::cout << "================ globals ==============================" << std::endl;
        // globals.Print("V");
        // An alternative to the above is to the S0_obs and B0_obs as you obtained them from the initial simulation.
        // https://indico.cern.ch/event/126652/contributions/1343592/attachments/80222/115004/Frequentist_Limit_Recommendation.pdf
        // The doc above says such an approach does not have good asymptotic properties.
        // globals.assign(*sb_global_default_vals);

        // 2. Main data
        mu.setVal(0.0); // we're generating bkg-only PSEUDOEXPERIMENT, so turn off the signal
        mu.setConstant(true);
        RooDataSet *ds = sb_full.generate({mass}, RooFit::Extended());
        ds->setGlobalObservables(globals); // set global observables so that RooFit does not fit S0_obs and B0_obs
        // std::cout << "================ ds ==============================" << std::endl;
        // ds->Print("V");


        // fit to s+b and record parameters and nll
        // set and fix mu to 1.0, meaning full signal strength
        sb_params.assign(*sb_params_default_vals);
        mu.setVal(1.0);
        mu.setConstant(true);
        RooFitResult *result_mu = sb_full.fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
        // std::cout << "================ mu=1 fit ==============================" << std::endl;
        // result_mu->Print("V");
        Double_t nll_mu = result_mu->minNll();
        RooArgSet *params_fit_sb = sb_params.snapshot();
        // params_fit_sb->Print("V");


        // fit to b-only and record parameters and nll
        // allow mu to float as described in the documents.
        mu.setVal(1E-5);
        mu.setConstant(false);
        RooFitResult *result_mu_hat = sb_full.fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                    RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                    RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
        // std::cout << "================ mu_hat fit ==============================" << std::endl;
        // result_mu_hat->Print("V");
        Double_t nll_mu_hat = result_mu_hat->minNll();
        // construct the test statistic of the "observed" data
        // Q = -2 ln (L(mu=1)/L(mu_hat)) = -2 * (ln(L(mu=1)) - ln(L(mu_hat))) = 2*(nll(mu=1) - nll(mu_hat))
        Double_t q_obs = 2. * (nll_mu - nll_mu_hat);
        

        // fit under pure background to find nuissance parameters
        mu.setVal(0.0);
        mu.setConstant(true);
        RooFitResult *result_b = sb_full.fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                               RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                               RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
        // std::cout << "================ bkg fit ==============================" << std::endl;
        // result_b->Print("V");
        RooArgSet *params_fit_b = sb_params.snapshot();
        // params_fit_b->Print("V");


        // Build distribution of q under b-only and count how many times the TOY has a higher q than the PSEUDOEXPERIMENT.
        // This number is 1-CLB (as defined in https://cds.cern.ch/record/1379837/files/NOTE2011_005.pdf)
        Double_t n_higher_bkg = 0.0;
        // histogram for the distribution of q under B-only
        // This will show that when using the "unconstrained" ensemble, the q_mu dsitributions do not depend on the pseudoexperiment
        // So one can derive them only once. In this code we do not use this property, but it's very nice to have in mind.
        TH1D *hqB_check = new TH1D(Form("hqB_check_%d", i), Form("q distr under bkg only, pseudoexp %d", i), 50, 0., 20.);
        for (Int_t i_toy = 0; i_toy < nToys; ++i_toy)
        {
            const auto startToy = std::chrono::high_resolution_clock::now();
            // OPTION 1. Set the bkg yield to the one obtained from the PSEUDOEXPERIMENT fit.
            // sb_params.assign(*params_fit_b);
            sb_params.assign(*params_fit_b);
            

            // generate the "XX_obs" variables
            RooDataSet *ds_global_toy = sb_full.generate(globals, 1);
            globals.assign(*ds_global_toy->get(0));
            
            mu.setVal(0.0);
            RooDataSet *ds_toy = sb_full.generate({mass}, RooFit::Extended());
            ds_toy->setGlobalObservables(globals);
            mu.setVal(1.0);
            mu.setConstant(true);
            RooFitResult *result_mu_toy = sb_full.fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                        RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                        RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
            Double_t nll_mu_toy = result_mu_toy->minNll();
            
            // mu.setVal(1E-5);
            mu.setConstant(false);
            RooFitResult *result_mu_hat_toy = sb_full.fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                            RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                            RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
            Double_t nll_mu_hat_toy = result_mu_hat_toy->minNll();

            Double_t q_toy = 2. * (nll_mu_toy - nll_mu_hat_toy);
            delete result_mu_toy;
            result_mu_toy = nullptr;
            delete result_mu_hat_toy;
            result_mu_hat_toy=nullptr;
            delete ds_toy;
            ds_toy = nullptr;
            delete ds_global_toy;
            ds_global_toy = nullptr;
            
            
            if (q_toy >= q_obs)
            {
                n_higher_bkg += 1.0;
            }

            hqB_check->Fill(q_toy);

            // remove the line below if you want to do debug printing.
            const auto stopToy = std::chrono::high_resolution_clock::now();
            auto durationToy = std::chrono::duration_cast<std::chrono::microseconds>(stopToy - startToy);
            auto durationPE = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startPseudoExp);
            auto durationNow = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startFull);

            Double_t durationToyNice = durationToy.count() / 1000.; // ms
            Double_t durationPENice = durationNow.count() / 1000.;  // s
            Double_t durationNowNice = durationNow.count() / 1000.; // s
            if ((i_toy+1) % (nToys / 2) == 0)
            {
                std::cout << "\r Ran pseudo-experiment " << i+1 << "/" << nPseudoExps
                          << ": bkg toy " << i_toy+1 << "/" << nToys
                          << " Duration: " << durationNowNice << "[s]" << std::flush;
            }
            hTimeBkgToy->Fill(durationToyNice);
        }

        TCanvas *cqB = new TCanvas(Form("cqB_%i", i),Form("cqB_%i", i), 800, 600);
        cqB->cd();
        hqB_check->Draw();
        TLine *line_obs = new TLine(q_obs, 0., q_obs, hqB_check->GetMaximum());
        line_obs->SetLineColor(kRed);
        line_obs->SetLineWidth(2);
        line_obs->Draw("SAME");
        output_file->cd();
        cqB->Write();

        delete cqB;
        cqB=nullptr;
        delete line_obs;
        line_obs=nullptr;
        delete hqB_check;
        hqB_check = nullptr;

        // Build distribution of q under S+B and count how many times the TOY has a higher q than the PSEUDOEXPERIMENT.
        // This number is CLSB (as defined in https://cds.cern.ch/record/1379837/files/NOTE2011_005.pdf)
        Double_t n_higher_sb = 0.0;
        TH1D *hqSB_check = new TH1D(Form("hqSB_check_%d", i), Form("q distr under sb, pseudoexp %d", i), 50, 0., 20.);
        for (Int_t i_toy = 0; i_toy < nToys; ++i_toy)
        {
            const auto startToy = std::chrono::high_resolution_clock::now();
            // OPTION 1. set the signal and bkg yield to the one obtained from the PSEUDOEXPERIMENT fit.
            sb_params.assign(*sb_params_default_vals);

            // sb_full.generate(b_params);
            // OPTION 2. Sample the signal and bkg yield before generating the PSEUDOEXPERIMENTS.
            // NOT a frequentist construction, but we should not fear this.
            RooDataSet *ds_global_toy = sb_full.generate(globals, 1);
            globals.assign(*ds_global_toy->get(0));
            mu.setVal(1.0);
            RooDataSet *ds_toy = sb_full.generate({mass}, RooFit::Extended());
            ds_toy->setGlobalObservables(globals);
            mu.setVal(1.0);
            mu.setConstant(true);
            RooFitResult *result_mu_toy = sb_full.fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                        RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                        RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
            Double_t nll_mu_toy = result_mu_toy->minNll();
            mu.setVal(1E-5);
            mu.setConstant(false);
            RooFitResult *result_mu_hat_toy = sb_full.fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                            RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                            RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
            Double_t nll_mu_hat_toy = result_mu_hat_toy->minNll();

            Double_t q_toy = 2. * (nll_mu_toy - nll_mu_hat_toy);

            delete result_mu_toy;
            result_mu_toy = nullptr;
            delete result_mu_hat_toy;
            result_mu_hat_toy=nullptr;
            delete ds_toy;
            ds_toy = nullptr;
            delete ds_global_toy;
            ds_global_toy = nullptr;
            

            if (q_toy >= q_obs)
            {
                n_higher_sb += 1.0;
            }
            hqSB_check->Fill(q_toy);

            // remove the line below if you want to do debug printing.
            const auto stopToy = std::chrono::high_resolution_clock::now();
            auto durationToy = std::chrono::duration_cast<std::chrono::microseconds>(stopToy - startToy);
            auto durationPE = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startPseudoExp);
            auto durationNow = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startFull);

            Double_t durationToyNice = durationToy.count() / 1000.; // ms
            Double_t durationPENice = durationNow.count() / 1000.;  // s
            Double_t durationNowNice = durationNow.count() / 1000.; // s
            if ((i_toy+1) % (nToys / 2) == 0)
            {
                std::cout << "\r Ran pseudo-experiment " << i+1 << "/" << nPseudoExps
                          << ": s+b toy " << i_toy+1 << "/" << nToys
                          << " Duration: " << durationNowNice << "[s]" << std::flush;
            }
            hTimeSBToy->Fill(durationToyNice);
        }

        TCanvas *cqSB = new TCanvas(Form("cqSB_%i", i),Form("cqSB_%i", i), 800, 600);
        cqSB->cd();
        hqSB_check->Draw();
        line_obs = new TLine(q_obs, 0., q_obs, hqSB_check->GetMaximum());
        line_obs->SetLineColor(kRed);
        line_obs->SetLineWidth(2);
        line_obs->Draw("SAME");
        output_file->cd();
        cqSB->Write();

        delete cqSB;
        cqSB=nullptr;
        delete line_obs;
        line_obs=nullptr;
        delete hqSB_check;
        hqSB_check = nullptr;

        // Computation of CLs. Exactly as described in the notes
        Double_t clSB = (1.0 * n_higher_sb) / (1.0 * nToys);
        Double_t clB = (1.0 * n_higher_bkg) / (1.0 * nToys);
        Double_t CLS = clB > 0. ? clSB / clB : 0.0;
        hCLS->Fill(CLS);
        hCLSB->Fill(clSB);
        hCLB->Fill(clB);
        // Assume 95% confindence level
        if (CLS < 0.05)
        {
            nTimesExcluded += 1.0;
        }
        nTotalSB += n_higher_sb;
        
        const auto stopPE = std::chrono::high_resolution_clock::now();
        auto durationPE = std::chrono::duration_cast<std::chrono::milliseconds>(stopPE - startPseudoExp);
        Double_t durationPENice = durationPE.count() / 1000.; // s
        hTimePseudoExp->Fill(durationPENice);
    }

    */

    std::cout << std::endl;
    Double_t ratioExcluded = (1.0*nTimesExcluded)/(1.0*nPseudoExps);
    std::cout << "Point excluded " << 100*ratioExcluded << "\% of the time.";
    if (ratioExcluded < 0.0228){
        std::cout << " It is outside the +-2sigma band; below the median";
    }
    else if (ratioExcluded < 0.1587){
        std::cout << "It is in the +-2sigma band, below the median";
    }
    else if (ratioExcluded < 0.5){
        std::cout << "It is in the +-1sigma band, below the median";
    }
    else if (ratioExcluded < 0.843){
        std::cout << "It is in the +-1sigma band, above the median";
    }
    else if (ratioExcluded < 0.9772){
        std::cout << "It is in the +- 2 sigma band, above the median";
    }
    else{
        std::cout << "It is outside the +-2sigma band, above the median";
    }
    std::cout << std::endl;
    
    prob_file << std::fixed << std::setprecision(2) << point.m_s << "," << std::setprecision(5) << std::scientific << nTotalSB/(nPseudoExps*nToys) << '\n';
    std::cout << std::fixed << std::setprecision(2);

    TCanvas *c_cls = new TCanvas("c_cls", "c_cls");
    c_cls->Divide(3, 1);
    c_cls->cd(1);
    hCLS->DrawCopy();
    c_cls->cd(2);
    hCLSB->DrawCopy();
    c_cls->cd(3);
    hCLB->DrawCopy();

    TCanvas *c_time = new TCanvas("c_time", "c_time", 1600, 600);
    c_time->Divide(3, 1);
    c_time->cd(1);
    hTimeBkgToy->DrawCopy();
    c_time->cd(2);
    hTimeSBToy->DrawCopy();
    c_time->cd(3);
    hTimePseudoExp->DrawCopy();

    output_file->Write();
    output_file->Close();
    std::cout << "Done\n\n" << std::endl;

    // Cleanup 
    delete cSigConstraint;
    cSigConstraint = nullptr;
    delete cBkgConstraint;
    cBkgConstraint = nullptr;
    delete c_cls;
    c_cls = nullptr;
    delete c_time;
    c_time = nullptr;

}


/*
    Main analysis function
*/
int main(int argc, char *argv[]) {
    
    // Make ROOT run in batch mode
    gROOT->SetBatch(1);

    ROOT::EnableThreadSafety();

    TApplication app("app", &argc, argv);

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // NB: change this to get different results every time.
    RooRandom::randomGenerator()->SetSeed(42);
    
    try {
        // Load the config file 
        const char* configFile = app.Argv(1);
        std::map<std::string, std::string> config = load_Config(configFile);

        // Check if the user entered the runType
        if(config["runType"] == "") {
            throw std::runtime_error("Error: Please specify the work mode of the program\n");
        }

        // Read the file containing the event yields
        discriminator = int(std::stod(config["discriminator"])*1000);
        const char* dataFile = Form("results/mChi2/signal_yields/sig_bkg_D%d.csv", discriminator);
        std::vector<DataPoint> data = read_CSV(dataFile);
        std::ofstream prob_file(Form("results/mChi2/roofit_results/out_D%d/p_values_S825.csv", discriminator));
        prob_file << "M_S,p_value\n";

        // Initialize the number of pseudo-experiments and toys
        nPseudoExps = std::stoi(config["nPseudoExps"]);
        nToys = std::stoi(config["nToys"]);

        // Check the runType for this instance
        if(config["runType"] == "full") {
            // Perform the analysis over the entire file
            std::cout << "Running full analysis over all data points\n\n";
            for(auto point : data)    
                analysisRun(point, prob_file);
        } 
        else if(config["runType"] == "point") {
            if(config["mass"] == "") {
                // If the mass point is not specified for analysis, exit with an error
                throw std::runtime_error("Error: Please specify a mass point to be analysed\n");
            }
            else {
                // Search for the index of the point with the specified mass
                double mass = std::stod(config["mass"]);
                auto it = std::find_if(data.begin(), data.end(),
                                    [mass](const DataPoint& p) {return p.m_s == mass;});
                if(it != data.end()) {
                    size_t idx = std::distance(data.begin(), it);
                    // Perform the analysis only on the specified point
                    std::cout << "Running single point analysis\n\n";
                    analysisRun(data[idx], prob_file);
                }
                else {
                    throw std::runtime_error("Error: Please input a valid mass point\n");
                }
            }
        }
        else throw std::runtime_error("Error: Please specify a valid work mode for the program\n");

        prob_file.close();

    }
    catch(const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        return EXIT_FAILURE;
    }


    // app.Run();
    app.Terminate();

    return 0;
}
