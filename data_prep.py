"""
Run this file as:

python data_prep.py discriminator_value

where discriminator_value = D * 1000. For example, for D = 0.925:
    python data_prep.py 925
"""


import numpy as np
import pandas as pd
import sys

# Take the discriminator value given in the arguments
D = int(sys.argv[1]) / 1000

# Mass value array
M_s = np.array([7.00, 7.25, 7.50, 7.75, 8.00, 8.25, 8.50, 8.75])
# 2D to keep the errors as well
sig = np.empty(shape=(len(M_s), 2))
bkg = np.empty(shape=(len(M_s), 2))

# Read the files for the different masses and build the 2D arrays
df = pd.read_csv("results/ATLAS_136_S700_B650_6j_5f/ATLAS_136_S700_B650_6j_5f_random_forest_counts_summary.csv")
sig[0, 0], sig[0, 1] = df.loc[29, str(D)].split("±")
bkg[0, 0], bkg[0, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S725_B675_6j_5f/ATLAS_136_S725_B675_6j_5f_random_forest_counts_summary.csv")
sig[1, 0], sig[1, 1] = df.loc[29, str(D)].split("±")
bkg[1, 0], bkg[1, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S750_B700_6j_5f/ATLAS_136_S750_B700_6j_5f_random_forest_counts_summary.csv")
sig[2, 0], sig[2, 1] = df.loc[29, str(D)].split("±")
bkg[2, 0], bkg[2, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S775_B725_6j_5f/ATLAS_136_S775_B725_6j_5f_random_forest_counts_summary.csv")
sig[3, 0], sig[3, 1] = df.loc[29, str(D)].split("±")
bkg[3, 0], bkg[3, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S800_B750_6j_5f/ATLAS_136_S800_B750_6j_5f_random_forest_counts_summary.csv")
sig[4, 0], sig[4, 1] = df.loc[29, str(D)].split("±")
bkg[4, 0], bkg[4, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S825_B775_6j_5f/ATLAS_136_S825_B775_6j_5f_random_forest_counts_summary.csv")
sig[5, 0], sig[5, 1] = df.loc[29, str(D)].split("±")
bkg[5, 0], bkg[5, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S850_B800_6j_5f/ATLAS_136_S850_B800_6j_5f_random_forest_counts_summary.csv")
sig[6, 0], sig[6, 1] = df.loc[29, str(D)].split("±")
bkg[6, 0], bkg[6, 1] = df.loc[30, str(D)].split("±")
df = pd.read_csv("results/ATLAS_136_S875_B825_6j_5f/ATLAS_136_S875_B825_6j_5f_random_forest_counts_summary.csv")
sig[7, 0], sig[7, 1] = df.loc[29, str(D)].split("±")
bkg[7, 0], bkg[7, 1] = df.loc[30, str(D)].split("±")


# Write the theoretical yields in a .csv file
df_out = pd.DataFrame({
    'M_S': M_s,
    'SIG': sig[:,0],
    "SIGMA_SIG": sig[:, 1],
    'BKG': bkg[:,0],
    "SIGMA_BKG": bkg[:, 1],
})
df_out.to_csv("results/signal_yields/sig_bkg_D{a}.csv".format(a = int(D*1000)), index=False)