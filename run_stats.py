import stats_models as mod
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
import argparse
from pathlib import Path
import sys
from datetime import datetime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.json")
    parser.add_argument("--paths", default="analysis_paths.json")
    args = parser.parse_args()

    config, path_file = mod.load_config(args.config, args.paths)

    yields_file = "/data/dcostach/diquark/"+path_file[config["process"]]+"/signal_yields/sig_bkg_D"+str(int(config["discriminator"]*1000))+".csv"
    df_yields = pd.read_csv(yields_file)
    uncrt_file = "/data/dcostach/diquark/"+path_file[config["process"]]+"/signal_yields/uncrt.csv"
    df_uncrt = pd.read_csv(uncrt_file)
    output_path = "/data/dcostach/diquark/"+path_file[config["process"]]+"/pyhf_output/out_D"+str(int(config["discriminator"]*1000))+"/"

    handlers = {
        ("cls", "asymptotic"): mod.run_cls_asymptotics,
        ("cls", "toybased"): lambda yld, unc, out: mod.run_cls_toybased(yld, unc, out, config["nToys"]),
        ("limits", "asymptotic"): mod.run_limits_asymptotics,
        ("limits", "toybased"): lambda yld, unc, out: mod.run_limits_toybased(yld, unc, out, config["nToys"]),
        ("pval", "asymptotic"): mod.run_pval_asymptotics,
        ("pval", "toybased"): lambda yld, unc, out: mod.run_pval_toybased(yld, unc, out, config["nToys"]),
    }

    key = (config["runType"], config["calcType"])
    runner = handlers[key]

    results = runner(df_yields, df_uncrt, output_path)

    # print(results)



if __name__ == "__main__":
    main()