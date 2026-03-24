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

struct PseudoExperimentInput {
    DataPoint point;
    SignalUncertainties uncrt;
    int experimentIndex;
};

struct PseudoExperimentResult {
    int experimentIndex;
    TH1D hCLS;
    TH1D hCLSB;
    TH1D hCLB;
    double nTimesExcluded;
    double nTotalSB;
    double nTotalB;
};

