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
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooRandom.h"
#include "RooMinimizer.h"
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

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
vector<double> point_exclusion(DataPoint point) {
    vector<double> exclusion_limits;

    RooRandom::randomGenerator()->SetSeed(69);

    // Create workspace
    RooWorkspace* wspace = new RooWorkspace("wspace");
    
    // Define sig constraints
    RooRealVar theta_S("theta_S", "theta_S", -5., 5.);
    wspace->import(theta_S);
    RooGaussian gamma_S("gamma_S", "gamma_S", theta_S, RooConst(0.), RooConst(1.));

    // Define bkg constraints
    RooRealVar theta_B("theta_B", "theta_B", -5., 5.);
    wspace->import(theta_B);
    RooGaussian gamma_B("gamma_B", "gamma_B", theta_B, RooConst(0.), RooConst(1.));

    // Define expected sig counts
    RooRealVar mu("mu", "mu", 1., 0., 5.);
    wspace->import(mu);
    RooRealVar S0("S0", "S0", point.sig);
    wspace->import(S0);
    RooRealVar sigma_S_frac("sigma_S", "sigma_S", point.sigma_sig/point.sig);
    wspace->import(sigma_S_frac);
    RooFormulaVar S_exp("S_exp", 
                        "@0 * @1 * (1 + @2 * @3)",
                        {mu, S0, theta_S, sigma_S_frac});

    // Define expected bkg counts
    RooRealVar B0("B0", "B0", point.bkg);
    wspace->import(B0);
    RooRealVar sigma_B_frac("sigma_B", "sigma_B", point.sigma_bkg/point.bkg);
    wspace->import(sigma_B_frac);
    RooFormulaVar B_exp("B_exp", 
                        "@0 * (1 + @1 * @2)",
                        {B0, theta_B, sigma_B_frac});

    // Define total expected counts
    RooFormulaVar N_exp("N_exp",
                        "@0 + @1",
                        {S_exp, B_exp});

    // Define Poisson distribution of events
    RooRealVar N_obs("N_obs", "N_obs", 0., 70.);
    wspace->import(N_obs);
    RooPoisson count("count", "count", N_obs, N_exp);

    // Define full model
    RooProdPdf model("model", "model", {count, gamma_S, gamma_B});
    wspace->import(model);

    // Define set of observables
    wspace->defineSet("obs", "N_obs");
    // Define set of parameters of interest
    wspace->defineSet("poi", "mu");
    // Define set of nuisance parameters
    wspace->defineSet("np", "theta_S,theta_B");
 

    // Create s+b model configuration
    ModelConfig* sbModel = new ModelConfig("sbModel", wspace);
    // Add parameters to the model configuration
    sbModel->SetPdf(model);
    sbModel->SetParametersOfInterest(mu);
    sbModel->SetNuisanceParameters({theta_S, theta_B});
    sbModel->SetObservables(N_obs);
    sbModel->SetSnapshot(mu);


    // Create b-only model configuration 
    ModelConfig* bModel = new ModelConfig("bModel", wspace);
    bModel->SetPdf(*wspace->pdf("model"));
    RooRealVar bPoi = mu;
    bPoi.setVal(0.);
    bPoi.setConstant(kTRUE);
    bModel->SetParametersOfInterest(bPoi);
    bModel->SetNuisanceParameters({theta_S, theta_B});
    bModel->SetObservables(N_obs);
    bModel->SetSnapshot(bPoi);



    for(int i = 0; i < 100; i++) {

        // Create dataset       
        wspace->var("mu")->setVal(0.);
        wspace->var("theta_S")->setVal(0.);
        wspace->var("theta_B")->setVal(0.);
        RooDataSet* toyData = wspace->pdf("model")->generate(N_obs, 1);
        toyData->Print("V");
        // wspace->pdf("model")->fitTo(*toyData);

        ProfileLikelihoodCalculator* plrCalc = new ProfileLikelihoodCalculator(*toyData, *sbModel);
        plrCalc->SetConfidenceLevel(0.9);

        LikelihoodInterval* interval = plrCalc->GetInterval();
        RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
        double lower = interval->LowerLimit(*poi);
        double upper = interval->UpperLimit(*poi);
        cout << "RESULT: " << 100*plrCalc->ConfidenceLevel() << "% interval is : ["<< lower << ", "<< upper <<"] "<< endl ;

        // Clear stack
        toyData->Delete();
        interval->Delete();

    }
    
    // Create asymptotic calculator 
    // AsymptoticCalculator* asympCalc = new AsymptoticCalculator(*toyData, *sbModel, *bModel);
    // asympCalc->SetOneSided(true);


    // HypoTestInverter* inverter = new HypoTestInverter(*asympCalc);
    // inverter->SetConfidenceLevel(0.95);
    // inverter->UseCLs(true);
    // inverter->SetVerbose(true);
    // inverter->SetFixedScan(50, 0, 5); 
    // inverter->SetAutoScan();

    // HypoTestInverterResult* result = inverter->GetInterval();
    
    // exclusion_limits.push_back(result->UpperLimit());
    // exclusion_limits.push_back(result->GetExpectedUpperLimit(-2));
    // exclusion_limits.push_back(result->GetExpectedUpperLimit(-1));
    // exclusion_limits.push_back(result->GetExpectedUpperLimit(0));
    // exclusion_limits.push_back(result->GetExpectedUpperLimit(1));
    // exclusion_limits.push_back(result->GetExpectedUpperLimit(2));

    // // Plot results
    // TCanvas* c = new TCanvas();
    // HypoTestInverterPlot* plot = new HypoTestInverterPlot("HTI_Result_Plot", "HypoTest Scan Results", result);
    // plot->Draw("CLB 2CL");
    // c->Draw();      
    
    
    // Clear the stack
    wspace->Delete();
    sbModel->Delete();
    bModel->Delete();

    return exclusion_limits;

}


/*
    Main analysis function
*/
void stat_analysis(const char* inputFile = nullptr) {
    vector<DataPoint> data = read_CSV(inputFile);
    vector<vector <double>> limits;

    limits.push_back(point_exclusion(data[3]));

    // for(auto point : data) 
    //     limits.push_back(point_exclusion(point));


    // // Print results
    // cout << "95% upper limit: " << limits[0][0] << '\n';
    // cout << "Expected upper limits, in the B (alternative) model:\n";
    // cout << "Expected limit (median): " << limits[0][3] << '\n';
    // cout << "Expected limit (-1 sig): " << limits[0][2] << '\n';
    // cout << "Expected limit (+1 sig): " << limits[0][4] << '\n';
    // cout << "Expected limit (-2 sig): " << limits[0][1] << '\n';
    // cout << "Expected limit (+2 sig): " << limits[0][5] << '\n';


}