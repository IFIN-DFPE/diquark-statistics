// parallel_generate.cxx
#include "ROOT/TProcessExecutor.hxx"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooRandom.h"

#include <iostream>
#include <fstream>
#include <memory>

using namespace ROOT;

std::string runToy(int seed) {
    // Create a workspace and PDF *inside each process*
    RooWorkspace w("w", "workspace");

    RooRealVar x("x", "observable", -10, 10);
    RooRealVar mean("mean", "mean", 0);
    RooRealVar sigma("sigma", "sigma", 2);

    RooGaussian gauss("g", "gaussian PDF", x, mean, sigma);
    w.import(gauss);

    RooAbsPdf *pdf = w.pdf("g");
    if (!pdf) throw std::runtime_error("PDF not found in workspace");

    RooRandom::randomGenerator()->SetSeed(seed);

    std::unique_ptr<RooDataSet> data(pdf->generate(RooArgSet(x), 1000));

    // Write out a simple text file with the mean of the sample
    double avg = 0;
    for (int i = 0; i < data->numEntries(); ++i) {
        const RooArgSet *row = data->get(i);
        avg += ((RooRealVar*)row->find("x"))->getVal();
    }
    avg /= data->numEntries();

    std::string fname = "toy_" + std::to_string(seed) + ".txt";
    std::ofstream out(fname);
    out << "Seed " << seed << " average = " << avg << "\n";
    out.close();

    return fname; // return filename to main process
}

int main() {
    // Seeds for each job
    std::vector<int> seeds = {101, 202, 303, 404};

    // Create a pool of processes (size = number of CPU cores by default)
    TProcessExecutor pool;

    // Run in parallel: each process executes runToy(seed)
    auto outputs = pool.Map(runToy, seeds);

    // Collect results
    for (auto &f : outputs) {
        std::cout << "Output written to: " << f << std::endl;
    }

    return 0;
}
