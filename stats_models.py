import pyhf
pyhf.set_backend("jax")
import numpy as np
import scipy.stats as stats
import pandas as pd
import warnings
import json
from dataclasses import dataclass
from pathlib import Path
warnings.filterwarnings("ignore", message="Values in x were outside bounds during a minimize step")


@dataclass
class CLsResult:
    observed: float
    expected: np.ndarray
@dataclass
class LimitResult:
    observed: float
    expected: np.ndarray



def load_config(config_file, path_file):
    """
        Function that loads the configuration file.

        Returns: 
            - a dictionary containing the configuration parameters.
            - a dictionary containing the paths to the input files.
    """
    with open(config_file, 'r') as f:
        config = json.load(f)
    with open(path_file, 'r') as f:
        path_file = json.load(f)
    return config, path_file



def build_pyhf_model(yields, uncrt, obs_yield):
    """
        Function that builds a pyhf json model starting from the 
        yields and uncertainties provided as arguments.

        Returns: 
            - a pyhf model specification.
    """
    uncrt_lumi = {"name": "lumi", "type": "normsys",
                  "data": {"hi": 1 + uncrt['lumi_uncrt']/100, "lo": 1 - uncrt['lumi_uncrt']/100}}
    
    uncrt_jes_jer = {"name": "jes_jer", "type": "normsys",
                "data": {"hi": 1 + uncrt['JES_JER_uncrt']/100, "lo": 1 - uncrt['JES_JER_uncrt']/100}}
    
    uncrt_cone = {"name": "cone", "type": "normsys",
                  "data": {"hi": 1 + uncrt['cone_uncrt']/100, "lo": 1 - uncrt['cone_uncrt']/100}}
    
    uncrt_pileup = {"name": "pileup", "type": "normsys",
                    "data": {"hi": 1 + uncrt['pileup_uncrt']/100, "lo": 1 - uncrt['pileup_uncrt']/100}}
    
    model_spec = {
        "channels": [
            {
                "name": "region1",
                "samples": [
                    {
                        "name": "signal",
                        "data": [yields['SIG']],
                        "modifiers": [
                            {"name": "mu", "type": "normfactor", "data": None},

                            {"name": "sig_ml", "type": "normsys", "data": {"hi": 1 + yields['SIGMA_SIG']/yields['SIG'], "lo": 1 - yields['SIGMA_SIG']/yields['SIG']}},

                            {"name": "sig_pdf", "type": "normsys", "data": {"hi": 1 + uncrt['sig_PDF']/100, "lo": 1 - uncrt['sig_PDF']/100}},

                            {"name": "sig_scale", "type": "normsys", "data": {"hi": 1 + uncrt['sig_scale_hi']/100, "lo": 1 - uncrt['sig_scale_lo']/100}},

                            uncrt_lumi, uncrt_jes_jer, uncrt_cone, uncrt_pileup
                        ]
                    },

                    {
                        "name": "hjj",
                        "data": [yields['BKG_JJH']],
                        "modifiers": [
                            {"name": "hjj_ml", "type": "normsys", "data": {"hi": 1 + yields['SIGMA_BKG_JJH']/yields['BKG_JJH'], "lo": 1 - yields['SIGMA_BKG_JJH']/yields['BKG_JJH']}},

                            {"name": "hjj_pdf", "type": "normsys", "data": {"hi": 1 + uncrt['hjj_PDF']/100, "lo": 1 - uncrt['hjj_PDF']/100}},

                            {"name": "hjj_scale", "type": "normsys", "data": {"hi": 1 + uncrt['hjj_scale_hi']/100, "lo": 1 - uncrt['hjj_scale_lo']/100}},

                            uncrt_lumi, uncrt_jes_jer, uncrt_cone, uncrt_pileup
                        ]
                    },

                    {
                        "name": "wj",
                        "data": [yields['BKG_WJ']],
                        "modifiers": [
                            {"name": "wj_ml", "type": "normsys", "data": {"hi": 1 + yields['SIGMA_BKG_WJ']/yields['BKG_WJ'], "lo": 1 - yields['SIGMA_BKG_WJ']/yields['BKG_WJ']}},

                            {"name": "wj_pdf", "type": "normsys", "data": {"hi": 1 + uncrt['Wj_PDF']/100, "lo": 1 - uncrt['Wj_PDF']/100}},

                            {"name": "wj_scale", "type": "normsys", "data": {"hi": 1 + uncrt['Wj_scale_hi']/100, "lo": 1 - uncrt['Wj_scale_lo']/100}},

                            uncrt_lumi, uncrt_jes_jer, uncrt_cone, uncrt_pileup
                        ]
                    },

                    {
                        "name": "qq2gg",
                        "data": [yields['BKG_QQ2GG']],
                        "modifiers": [
                            {"name": "qq2gg_ml", "type": "normsys", "data": {"hi": 1 + yields['SIGMA_BKG_QQ2GG']/yields['BKG_QQ2GG'], "lo": 1 - yields['SIGMA_BKG_QQ2GG']/yields['BKG_QQ2GG']}},

                            {"name": "qq2gg_pdf", "type": "normsys", "data": {"hi": 1 + uncrt['qq2gg_PDF']/100, "lo": 1 - uncrt['qq2gg_PDF']/100}},

                            {"name": "qq2gg_scale", "type": "normsys", "data": {"hi": 1 + uncrt['qq2gg_scale_hi']/100, "lo": 1 - uncrt['qq2gg_scale_lo']/100}},

                            uncrt_lumi, uncrt_jes_jer, uncrt_cone, uncrt_pileup
                        ]
                    }
                ]
            }
        ],

        "measurements": [
            {
                "name": "measurement1",
                "config": {
                    "poi": "mu",
                    "parameters": [
                        {
                            "name": "mu",
                            "bounds": [[0.0, 10.0]]
                         }
                    ]
                }
            }
        ],

        "observations": [
            {
                "name": "region1",
                "data": [obs_yield]
            }
        ],

        "version": "1.0.0"

    }

    return model_spec



def run_cls_asymptotics(df_yields, df_uncrt, output_path):
    """
        Function that runs asymptotic CLs scan.

        Returns:    
            - Observed CLs
            - Expected CLs band array (-2s, -1s, median, +1s, +2s)
    """
    rows = []
    results = []

    for i in range(len(df_yields)):
        yields = df_yields.iloc[i].to_dict()
        uncrt = df_uncrt.iloc[i].to_dict()

        expected_bkg_total = sum([yields['BKG_JJH'], yields['BKG_WJ'], yields['BKG_QQ2GG']])

        model_dict = build_pyhf_model(yields, uncrt, expected_bkg_total)
        workspace = pyhf.Workspace(model_dict)
        model = workspace.model()

        obs_data = workspace.data(model)
        
        print(f"\n===============================================================")
        print(f"Running asymptotic CLs computation for:")
        print(f"M_S: {yields['M_S']} TeV")

        CLs_obs, CLs_bands = pyhf.infer.hypotest(
            poi_test = 1.0,
            data = obs_data,
            pdf = model,
            test_stat = 'qtilde',
            calctype = 'asymptotics',
            return_expected_set = True
        )

        print(f"Observed CLs (Median)         : {CLs_obs}")
        print(f"Expected -1 Sigma             : {CLs_bands[1]}")
        print(f"Expected CLs (Median)         : {CLs_bands[2]}")
        print(f"Expected +1 Sigma             : {CLs_bands[3]}")

        print(f"Done.")
        print(f"===============================================================")

        rows.append({
            "M_S": yields["M_S"],
            "obs": CLs_obs,
            "exp-2s": CLs_bands[0],
            "exp-1": CLs_bands[1],
            "exp_med": CLs_bands[2],
            "exp+1s": CLs_bands[3],
            "exp+2s": CLs_bands[4],
        })

        results.append(CLsResult(
            observed = float(np.asarray(CLs_obs).reshape(-1)[0]),
            expected = np.asarray(CLs_bands)
        ))

    out = Path(output_path+"cls_asymptotic_results.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    return results



def run_cls_toybased(df_yields, df_uncrt, output_path, toys):
    """
        Function that runs toy-based CLs scan.

        Returns:    
            - Observed CLs
            - Expected CLs band array (-2s, -1s, median, +1s, +2s)
    """
    rows = []

    for i in range(len(df_yields)):
        yields = df_yields.iloc[i].to_dict()
        uncrt = df_uncrt.iloc[i].to_dict()
        expected_bkg_total = sum([yields['BKG_JJH'], yields['BKG_WJ'], yields['BKG_QQ2GG']])

        model_dict = build_pyhf_model(yields, uncrt, expected_bkg_total)
        workspace = pyhf.Workspace(model_dict)
        model = workspace.model()

        obs_data = workspace.data(model)
        
        print(f"\n===============================================================")
        print(f"Running toy-based CLs computation for:")
        print(f"M_S: {yields['M_S']} TeV")

        CLs_obs, CLs_bands = pyhf.infer.hypotest(
            poi_test = 1.0,
            data = obs_data,
            pdf = model,
            test_stat = 'qtilde',
            calctype = 'toybased',
            ntoys = toys,
            return_expected_set = True
        )

        print(f"Observed CLs (Median)         : {CLs_obs}")
        print(f"Expected -1 Sigma             : {CLs_bands[1]}")
        print(f"Expected CLs (Median)         : {CLs_bands[2]}")
        print(f"Expected +1 Sigma             : {CLs_bands[3]}")

        print(f"Done.")
        print(f"===============================================================")

        rows.append({
                "M_S": yields["M_S"],
                "obs": CLs_obs,
                "exp-2s": CLs_bands[0],
                "exp-1": CLs_bands[1],
                "exp_med": CLs_bands[2],
                "exp+1s": CLs_bands[3],
                "exp+2s": CLs_bands[4],
            })

    out = Path(output_path+"cls_toybased_results.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    return CLsResult(
        observed = float(CLs_obs),
        expected = np.asarray(CLs_bands)
    )



def run_limits_asymptotics(df_yields, df_uncrt, output_path):
    """
        Function that runs asymptotic 95% Upper Limit scan.

        Returns:    
            - Observed Upper Limit
            - Expected Upper Limit band array (-2s, -1s, median, +1s, +2s)
    """
    rows = []

    for i in range(len(df_yields)):
        yields = df_yields.iloc[i].to_dict()
        uncrt = df_uncrt.iloc[i].to_dict()
        expected_bkg_total = sum([yields['BKG_JJH'], yields['BKG_WJ'], yields['BKG_QQ2GG']])

        model_dict = build_pyhf_model(yields, uncrt, expected_bkg_total)
        workspace = pyhf.Workspace(model_dict)
        model = workspace.model()

        obs_data = workspace.data(model)

        print(f"\n===============================================================")
        print(f"Running asymptotic Upper Limit computation for:")
        print(f"M_S: {yields['M_S']} TeV")

        obs_limit, exp_limits = pyhf.infer.intervals.upper_limits.upper_limit(
            obs_data,
            model,
            scan = np.linspace(0, 10, 50),
            test_stat = 'qtilde',
            calctype = 'asymptotics',
            level = 0.05,
        )

        print(f"Observed Upper Limit (Median) : {obs_limit}")
        print(f"Expected -1 Sigma             : {exp_limits[1]*yields['SIG']}")
        print(f"Expected Upper Limit (Median) : {exp_limits[2]*yields['SIG']}")
        print(f"Expected +1 Sigma             : {exp_limits[3]*yields['SIG']}")

        print(f"Done.")
        print(f"===============================================================")

        rows.append({
            "M_S": yields["M_S"],
            "obs": obs_limit,
            "exp-2s": exp_limits[0]*yields['SIG'],
            "exp-1s": exp_limits[1]*yields['SIG'],
            "exp_med": exp_limits[2]*yields['SIG'],
            "exp+1s": exp_limits[3]*yields['SIG'],
            "exp+2s": exp_limits[4]*yields['SIG'],
        })

    out = Path(output_path+"limits_asymptotic_results.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    return LimitResult(
        observed = float(obs_limit),
        expected = np.asarray(exp_limits)*yields['SIG']
    )



def run_limits_toybased(df_yields, df_uncrt, output_path, ntoys):
    """
        Function that runs toy-based 95% Upper Limit scan.

        Returns:    
            - Observed Upper Limit
            - Expected Upper Limit band array (-2s, -1s, median, +1s, +2s)
    """
    rows = []

    for i in range(len(df_yields)):
        yields = df_yields.iloc[i].to_dict()
        uncrt = df_uncrt.iloc[i].to_dict()
        expected_bkg_total = sum([yields['BKG_JJH'], yields['BKG_WJ'], yields['BKG_QQ2GG']])

        model_dict = build_pyhf_model(yields, uncrt, expected_bkg_total)
        workspace = pyhf.Workspace(model_dict)
        model = workspace.model()

        obs_data = workspace.data(model)

        print(f"\n===============================================================")
        print(f"Running toy-based Upper Limit computation for:")
        print(f"M_S: {yields['M_S']} TeV")

        obs_limit, exp_limits = pyhf.infer.intervals.upper_limits.upper_limit(
            obs_data,
            model,
            scan = np.linspace(0, 10, 50),
            test_stat = 'qtilde',
            calctype = 'toybased',
            ntoys = ntoys,
            level = 0.05,
        )

        print(f"Observed Upper Limit (Median) : {obs_limit}")
        print(f"Expected -1 Sigma             : {exp_limits[1]*yields['SIG']}")
        print(f"Expected Upper Limit (Median) : {exp_limits[2]*yields['SIG']}")
        print(f"Expected +1 Sigma             : {exp_limits[3]*yields['SIG']}")

        print(f"Done.")
        print(f"===============================================================")

        rows.append({
                "M_S": yields["M_S"],
                "obs": obs_limit,
                "exp-2s": exp_limits[0]*yields['SIG'],
                "exp-1s": exp_limits[1]*yields['SIG'],
                "exp_med": exp_limits[2]*yields['SIG'],
                "exp+1s": exp_limits[3]*yields['SIG'],
                "exp+2s": exp_limits[4]*yields['SIG'],
            })

    out = Path(output_path+"limits_toybased_results.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    return LimitResult(
        observed = float(obs_limit),
        expected = np.asarray(exp_limits)*yields['SIG']
    )



def run_pval_asymptotics(df_yields, df_uncrt, output_path):
    """
        Function that runs asymptotic p-value computation.

        Returns:    
            - Expected p-value
    """
    rows = []

    for i in range(len(df_yields)):
        yields = df_yields.iloc[i].to_dict()
        uncrt = df_uncrt.iloc[i].to_dict()
        expected_bkg_total = sum([yields['BKG_JJH'], yields['BKG_WJ'], yields['BKG_QQ2GG']])

        model_dict = build_pyhf_model(yields, uncrt, expected_bkg_total)
        workspace = pyhf.Workspace(model_dict)
        model = workspace.model()

        init_pars = model.config.suggested_init()
        poi_index = model.config.poi_index
        init_pars[poi_index] = 1.0

        obs_data = model.expected_data(init_pars)

        print(f"\n===============================================================")
        print(f"Running asymptotic p-value computation for:")
        print(f"M_S: {yields['M_S']} TeV")

        p0_exp_tmp = pyhf.infer.hypotest(
            0.0,
            obs_data,
            model,
            test_stat='q0',
            calctype='asymptotics'
        )

        p0_exp = float(p0_exp_tmp)

        Z_exp = stats.norm.isf(p0_exp)

        # print(f"Observed p-value: {p0_obs:.3e} (Z = {Z_obs:.2f} sigma)")
        print(f"Expected p-value: {p0_exp:.3e} (Z = {Z_exp:.2f} sigma)")

        print(f"Done.")
        print(f"===============================================================")

        rows.append({
            "M_S": yields["M_S"],
            "p0_exp": p0_exp,
            "Z_exp": Z_exp
        })

    out = Path(output_path+"pval_asymptotics_results.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    return p0_exp



def run_pval_toybased(df_yields, df_uncrt, output_path, ntoys):
    """
        Function that runs toy-based p-value computation.

        Returns:    
            - Expected p-value
    """
    rows = []

    for i in range(len(df_yields)):
        yields = df_yields.iloc[i].to_dict()
        uncrt = df_uncrt.iloc[i].to_dict()
        expected_bkg_total = sum([yields['BKG_JJH'], yields['BKG_WJ'], yields['BKG_QQ2GG']])

        model_dict = build_pyhf_model(yields, uncrt, expected_bkg_total)
        workspace = pyhf.Workspace(model_dict)
        model = workspace.model()

        init_pars = model.config.suggested_init()
        poi_index = model.config.poi_index
        init_pars[poi_index] = 1.0

        obs_data = model.expected_data(init_pars)

        print(f"\n===============================================================")
        print(f"Running toy-based p-value computation for:")
        print(f"M_S: {yields['M_S']} TeV")

        p0_exp_tmp = pyhf.infer.hypotest(
            0.0,
            obs_data,
            model,
            test_stat='q0',
            calctype='toybased',
            ntoys = ntoys
        )

        p0_exp = float(p0_exp_tmp)

        Z_exp = stats.norm.isf(p0_exp)

        print(f"Expected p-value: {p0_exp:.3e} (Z = {Z_exp:.2f} sigma)")

        print(f"Done.")
        print(f"===============================================================")

        rows.append({
            "M_S": yields["M_S"],
            "p0_exp": p0_exp,
            "Z_exp": Z_exp
        })

    out = Path(output_path+"pval_toybased_results.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    return p0_exp