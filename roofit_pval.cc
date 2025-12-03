#include <chrono>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <filesystem>
#include <omp.h>

#include "nlohmann/json.hpp"
#include "ROOT/TProcessExecutor.hxx"
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

#include "data_types.hxx"


/*
======================== RUNNING FROM SCRIPT ================================
    It is enough to run the command:
    ./run_pval.sh
    and the program will be compiled and run with the parameters given in stats.ini
*/


// Value of the discriminator used in ML
int discriminator;
// Number of PSEUDOEXPERIMENTS. Recommend 100 at least, maybe 1000
int nPseudoExps = 100;
// Number of TOYS per PSEUDOEXPERIMENTS. Recommend 1000 at least, maybe 10000
int nToys = 1000;
// Path to the process 
std::string path;



/*
    Function taking as argument the path of a .csv datafile and
    returns a vector of DataPoint type values.
*/
std::vector<DataPoint> read_CSV(std::string inputFile) {
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



PseudoExperimentResult runPseudoExp(PseudoExperimentInput input) {
    DataPoint point = input.point;
    int index = input.experimentIndex;

    // Use experiment index as RNG seed
    RooRandom::randomGenerator()->SetSeed(index);

    const auto startPseudoExp = std::chrono::high_resolution_clock::now();

    /*
        Define a dummy "discriminator" variable, for example reconstructed mass.
        We will assume it is uniformly distributed for both signal and background.
        In this case it will simply not matter that we have it, but the construction
        is easier to understand when all ingredients are there.
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

    // OPTION 1. Set the XX_true parameters to the default values in order to generate PSEUDOEXPERIMENTS the same way every-time.

    sb_params.assign(*sb_params_default_vals);

    // sb_full->Print("v");

    // 1. Global observables
    // We regenerate those because if we were to "repeat" the experiment, you would redo the simulations as well, obtaining different
    // S0_obs and B0_obs.
    RooDataSet* ds_global = sb_full.generate(globals, 1);

    globals.assign(*ds_global->get(0));

    // 2. Main data
    mu.setVal(1.0); // we're generating s+b PSEUDOEXPERIMENT
    mu.setConstant(true);
    RooDataSet *ds = sb_full.generate({mass}, RooFit::Extended());
    ds->setGlobalObservables(globals); // set global observables so that RooFit does not fit S0_obs and B0_obs


    // fit to b and record parameters and nll
    // set and fix mu to 0.0, meaning background only
    sb_params.assign(*sb_params_default_vals);
    mu.setVal(0.0);
    mu.setConstant(true);

    RooFitResult *result_0;
    result_0 = sb_full.fitTo(*ds, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                            RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                            RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));

    // std::cout << "================ mu=1 fit ==============================" << std::endl;
    // result_mu->Print("V");
    Double_t nll_0 = result_0->minNll();
    RooArgSet *params_fit_sb = sb_params.snapshot();
    // params_fit_sb->Print("V");


    // fit to s+b and record parameters and nll
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
    // Q = -2 ln (L(mu=0)/L(mu_hat)) = -2 * (ln(L(mu=0)) - ln(L(mu_hat))) = 2*(nll(mu=0) - nll(mu_hat))
    Double_t q_obs = 2. * (nll_0 - nll_mu_hat);


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
    Double_t n_higher_bkg = 0.0;

    std::cout << "Starting pseudo-experiment #" << index << "\n" << std::flush;
    const auto startFull = std::chrono::high_resolution_clock::now();
    for (Int_t i_toy = 0; i_toy < nToys; ++i_toy)
    {
        const auto startToy = std::chrono::high_resolution_clock::now();
        // OPTION 1. Set the bkg yield to the one obtained from the PSEUDOEXPERIMENT fit.
        sb_params.assign(*params_fit_b);


        // generate the "XX_obs" variables
        RooDataSet *ds_global_toy = sb_full.generate(globals, 1);
        globals.assign(*ds_global_toy->get(0));

        mu.setVal(0.0);
        RooDataSet *ds_toy = sb_full.generate({mass}, RooFit::Extended());
        ds_toy->setGlobalObservables(globals);

        // Fit under the mu = 0 hypothesis
        mu.setVal(0.0);
        mu.setConstant(true);
        RooFitResult *result_0_toy = sb_full.fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                    RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                    RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
        Double_t nll_0_toy = result_0_toy->minNll();

        // Fit under the mu_hat hypothesis
        // mu.setVal(1E-5);
        mu.setConstant(false);
        RooFitResult *result_mu_hat_toy = sb_full.fitTo(*ds_toy, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::PrintLevel(-1),
                                                        RooFit::PrintEvalErrors(-1), RooFit::Warnings(false),
                                                        RooFit::Verbose(false), RooFit::Save(), RooFit::GlobalObservables(globals));
        Double_t nll_mu_hat_toy = result_mu_hat_toy->minNll();

        Double_t q_toy = 2. * (nll_0_toy - nll_mu_hat_toy);
        delete result_0_toy;
        result_0_toy = nullptr;
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

        const auto stopToy = std::chrono::high_resolution_clock::now();
        auto durationToy = std::chrono::duration_cast<std::chrono::microseconds>(stopToy - startToy);
        auto durationPE = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startPseudoExp);
        auto durationNow = std::chrono::duration_cast<std::chrono::milliseconds>(stopToy - startFull);

        Double_t durationToyNice = durationToy.count() / 1000.; // ms
        Double_t durationPENice = durationNow.count() / 1000.;  // s
        Double_t durationNowNice = durationNow.count() / 1000.; // s
        if ((i_toy+1) == nToys)
        {
            std::cout << " Duration: " << durationNowNice << "[s]\n" << std::flush;
        }
    }

    

    PseudoExperimentResult result;

    result.nTotalB = n_higher_bkg;
    result.experimentIndex = index;



    return result;
}



/*
    Function that performs the statistical analysis and checks if the points are excluded or not
    As an argument, pass the DataPoint associated to a given mass of Suu
*/
void analysisRun(DataPoint point, std::ofstream &prob_file) {

    std::cout << Form("=== Running analysis for M_S = %.2f TeV ===\n", point.m_s);


    // Create a pool of processes (size = number of CPU cores by default)
    ROOT::TProcessExecutor pool(96);

    // Create a vector to store the inputs for each pseudo-experiment
    std::vector<PseudoExperimentInput> inputs(nPseudoExps);
    for (int i = 0; i < nPseudoExps; ++i) {
        inputs[i] = {point, i+1};
    }

    // Run in parallel: each process executes runPseudoExperiment
    PseudoExperimentResult results = pool.MapReduce(runPseudoExp, inputs, [] (const std::vector<PseudoExperimentResult>& results) {
        PseudoExperimentResult combined;
        combined.nTotalB = 0.0;

        for (const auto& res : results) {
            combined.nTotalB += res.nTotalB;
        }

        return combined;
    });

    // Counter to check how often we exclude the parameter point (in which the signal yield is what it is above).
    Double_t nTotalB = results.nTotalB;


    prob_file << std::fixed << std::setprecision(2) << point.m_s << "," << std::setprecision(5) << std::scientific << nTotalB/(nPseudoExps*nToys) << '\n';
    std::cout << std::fixed << std::setprecision(2);


    std::cout << "Done\n\n" << std::endl;
}


/*
    Main analysis function
*/
int main(int argc, char *argv[]) {

    // Make ROOT run in batch mode
    gROOT->SetBatch(1);

    TApplication app("app", &argc, argv);

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // NB: change this to get different results every time.
    RooRandom::randomGenerator()->SetSeed(42);

    // Let Cling know about our data structures
    gInterpreter->LoadFile("data_types.hxx");

    try {
        
        // Load the file containing the paths
        std::ifstream pathFile("analysis_paths.json");
        nlohmann::json paths = nlohmann::json::parse(pathFile);

        // Load the config file
        const char* configFilePath = app.Argv(1);
        std::ifstream configFile(configFilePath);
        nlohmann::json config = nlohmann::json::parse(configFile);
    
        // Check if the user entered the runType
        if(config["runType"] == "") {
            throw std::runtime_error("Error: Please specify the work mode of the program\n");
        }

        // Read the file containing the event yields
        discriminator = int(config["discriminator"].get<float>()*1000);
        auto process = config["process"].get<std::string>();
        path = paths[process].get<std::string>();
        std::string dataFile = path + Form("/signal_yields/sig_bkg_D%d.csv", discriminator);
        std::vector<DataPoint> data = read_CSV(dataFile);
        std::ofstream prob_file(path + Form("/roofit_results/out_D%d/p_values.csv", discriminator));
        prob_file << "M_S,p_value\n";

        // Initialize the number of pseudo-experiments and toys
        nPseudoExps = config["nPseudoExps"].get<int>();
        nToys = config["nToys"].get<int>();
        auto runType = config["runType"].get<std::string>();
        double mass = config["mass"].get<double>();

        // Check the runType for this instance
        if(runType == "full") {
            // Perform the analysis over the entire file
            std::cout << "Running full analysis over all data points\n\n";
            for(auto point : data)
                analysisRun(point, prob_file);
        }
        else if(runType == "point") {
            if(!mass) {
                // If the mass point is not specified for analysis, exit with an error
                throw std::runtime_error("Error: Please specify a mass point to be analysed\n");
            }
            else {
                // Search for the index of the point with the specified mass
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

        pathFile.close();
        configFile.close();
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
