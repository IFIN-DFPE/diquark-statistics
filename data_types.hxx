#pragma once

// We use this file to let ROOT/Cling know what kinds of data types we use for inter-process communication.

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

struct PseudoExperimentInput {
    DataPoint point;
    int experimentIndex;
};

struct PseudoExperimentResult {
    int experimentIndex;
    TH1D hCLS;
    TH1D hCLSB;
    TH1D hCLB;
    double nTimesExcluded;
    double nTotalSB;
};

