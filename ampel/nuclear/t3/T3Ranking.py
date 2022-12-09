#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t3/T3Ranking.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu
# Date              : 08.06.2020
# Last Modified Date: 29.11.2022
# Last Modified By  : simeon.reusch@desy.de

import datetime, json, io, time, os

from typing import Any, Optional, Union
from collections.abc import Generator

import matplotlib as mpl  # type: ignore
from matplotlib import pyplot as plt  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import astropy.time  # type: ignore
import corner  # type: ignore

from ampel.types import UBson, T3Send
from ampel.view.TransientView import TransientView
from ampel.view.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.nuclear.t3.dropboxIO import DropboxUnit
import ampel.nuclear.t3.classifyme as classifyme

flex_keys = [
    "dtime_peak",
    "band",
    "mag_peak",
    "sigma_rise",
    "sigma_fade",
    "mean_color",
    "color_slope",
    "e_dtime_peak",
    "e_mag_peak",
    "e_mean_color",
    "e_color_slope",
    "e_sigma_rise",
    "mad",
    "rms",
    "chi2",
    "n_rise",
    "n_fade",
    "classification",
    "name",
    "mjd_peak",
    "photoclass",
    "blabber",
]

simple_keys = [
    "name",
    "classification",
    "ra",
    "dec",
    "ndetections",
    "ndetections_negative",
    "peak_mag",
    "latest_mag",
    "peak_jd",
    "latest_jd",
    "age",
    "ndetections_g",
    "offset_rms_g",
    "offset_unc_g",
    "offset_med_g",
    "offset_weighted_g",
    "offset_sig_g",
    "ndetections_r",
    "offset_rms_r",
    "offset_unc_r",
    "offset_med_r",
    "offset_weighted_r",
    "offset_sig_r",
    "offset_sig",
    "lastest_obs_is_detection",
    "latest_maglim",
    "ref_mag_r",
    "ref_mag_g",
    "peak_mag_g",
    "latest_mag_g",
    "peak_diff_mag_g",
    "latest_diff_mag_g",
    "peak_mag_r",
    "latest_mag_r",
    "peak_diff_mag_r",
    "latest_diff_mag_r",
    "ref_mag",
    "latest_diff_mag",
    "peak_diff_mag",
]


class T3Ranking(DropboxUnit):

    plotHisto: bool = False  # create histrograms for of the magnitude of the transients in the different summary plots

    def post_init(self):
        # super().post_init()

        # make empty dataframes
        self.metrics = pd.DataFrame(columns=simple_keys)
        self.metrics_flex = pd.DataFrame(columns=flex_keys)

        # do folder I/O
        self.ranking_location = self.base_location + "/ranking"
        self.sum_location = self.base_location + "/sum_plots"

        if "ranking" not in self.get_files(self.base_location):
            self.create_folder(self.ranking_location)

        if "sum_plots" not in self.get_files(self.base_location):
            self.create_folder(self.sum_location)

        # if needed, add the year subfolder (always needed when we are running pytest)
        self.this_year = str(self.night.datetime.year)
        if self.this_year not in self.get_files(self.ranking_location):
            self.create_folder(self.ranking_location + "/" + self.this_year)

        sum_location_year = self.sum_location + "/" + self.this_year
        if self.this_year not in self.get_files(self.sum_location):
            self.create_folder(sum_location_year)

        # finally, add the dd-mm subfolder for sum_plot
        self.this_mmdd = self.night.datetime.strftime("%m-%d")
        self.plot_dir = sum_location_year + "/" + self.this_mmdd
        if self.this_mmdd not in self.get_files(sum_location_year):
            self.create_folder(self.plot_dir)

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        """ """
        transients = [t for t in gen]

        metrics, metrics_flex = self.collect_metrics(transients)

        if len(metrics):
            self.metrics = pd.concat([self.metrics, metrics], ignore_index=True)
            self.metrics_flex = pd.concat(
                [self.metrics_flex, metrics_flex], ignore_index=True
            )
            self.logger.info(
                "done adding transients, total collected: " + str(len(self.metrics))
            )
            self.make_sumplots(self.metrics, self.metrics_flex)
        else:
            self.metrics = []
            self.metrics_flex = []
            self.logger.info("metrics is empty, skipping summary plots")

        if len(self.metrics):
            self.apply_ranking(self.metrics, self.metrics_flex)

        super().done()

        return None

    def collect_metrics(self, transients):
        """ """

        if transients:

            local_sources = np.array(transients)
            newdetection_count = len(local_sources)
            self.logger.info(
                "Collecting metrics for {} transient(s)".format(newdetection_count)
            )

            simple_res_list = []
            flexfit_res_list = []

            for i, tran_view in enumerate(local_sources):

                # get the classification and some extra info
                # writing of the classifcation to dropbox is done by T3SummaryPlots
                classification, extra_info = classifyme.get_class(
                    tran_view, self, self.logger, write=False
                )

                simple_res = self.get_simple_results(
                    tran_view, classification, extra_info
                )
                flexfit_res = self.get_fit_params(tran_view, classification)

                simple_res_list.append(simple_res)
                flexfit_res_list.append(flexfit_res)

            metrics = pd.DataFrame.from_records(simple_res_list)
            metrics_flex = pd.DataFrame.from_records(flexfit_res_list)

            self.logger.info(
                "successfully collected metrics: "
                + str(all(metrics_flex["name"] == metrics["name"]))
            )
        else:
            self.logger.info("no transients to rank!")
            return [], []

        # Cast object columns to float64
        for k in metrics_flex.columns:
            if k not in {
                "name",
                "classification",
                "photoclass",
                "blabber",
                "n_rise",
                "n_fade",
                "band",
            }:
                metrics_flex[k] = metrics_flex[k].astype(float)

        # Returns results from T2s
        return metrics, metrics_flex

    def get_fit_params(self, tran_view, classification):
        """ """
        try:
            fit_params = dict(
                (tran_view.get_latest_t2_body("T2FlexFit") or {})["fit_params"]
            )
            fit_params["classification"] = classification
            return fit_params

        except KeyError:
            self.logger.warn(
                "WARNING: No Flex fit for source {}".format(tran_view.stock["name"][0])
            )
            emptyrow = dict(zip(np.full(len(flex_keys), None), flex_keys))
            emptyrow["name"], emptyrow["classification"] = (
                tran_view.stock["name"][0],
                classification,
            )
            return emptyrow

    def get_simple_results(self, tran_view, classification, extra_info):
        """ """
        simple_results = dict(
            tran_view.get_latest_t2_body("T2SimpleMetrics" or {}).get("metrics", {})
        )

        simple_results["ra"] = np.mean(simple_results["ra"])
        simple_results["dec"] = np.mean(simple_results["dec"])

        simple_results["classification"] = classification
        simple_results["extra_info"] = extra_info

        # get all the jds for the lightcurve photopoints -> might provide useful
        # lc_photopoints = tran_view.get_lightcurves()[0].get_photopoints()
        # lc_jds = [p.get("body").get("jd") for p in lc_pps]

        simple_results["age"] = (
            self.night.jd
            - astropy.time.Time(tran_view.get_time_updated(output="datetime")).jd
        )

        # to do: read this from the T2 catalog match records?
        alert_base = "/mampel/alerts/" + "20" + simple_results["name"][3:5]
        try:
            # if matched to PS1, read offset and our estimate of the offset significance
            # this could also be computed from the T2 match info, but this is fast enough
            results_infile = self.read_file(
                alert_base
                + "/{}/{}_simple.json".format(
                    simple_results["name"], simple_results["name"]
                )
            )
            if isinstance(results_infile, str):
                results_json = json.loads(results_infile)
            else:
                results_json = results_infile.json()

            simple_results["offset_med_ps1"] = results_json["offset_med_ps1"]
            simple_results["offset_rms_ps1"] = results_json["offset_rms_ps1"]
            simple_results["offset_sig_ps1"] = results_json["offset_sig_ps1"]
            self.logger.info(
                "PS1 offsets calculated for {}".format(simple_results["name"])
            )
        except (FileNotFoundError, KeyError) as e:
            self.logger.warn(
                "WARNING: No PS1 offsets calculated for source {}".format(
                    tran_view.stock["name"][0]
                )
            )
            simple_results["offset_med_ps1"] = 0.0
            simple_results["offset_rms_ps1"] = 0.0
            simple_results["offset_sig_ps1"] = 0.0

        return simple_results

    def apply_ranking(self, metrics, metrics_flex):
        """ """

        ranking = self.get_ranking(metrics, metrics_flex)
        self.metrics["ranking"] = ranking
        metrics = self.metrics
        metrics_flex = self.metrics_flex

        if self.dryRunDir:
            metrics.to_csv(os.path.join(self.dryRunDir, "metrics.csv"))
            metrics_flex.to_csv(os.path.join(self.dryRunDir, "metrics_flex.csv"))

        base_name = (
            self.ranking_location + "/" + self.this_year + "/" + self.this_mmdd
        )  # base of filename

        metrics["classification"] = [
            "None" if (c == "") or (c == np.nan) else c
            for c in metrics["classification"]
        ]

        # disable for now
        iage = metrics["age"] > -9999

        # Select things saved for as part of Ampel filter
        # iampel = ['ZTFBH Nuclear' in [x['comment'] for x in aa] for aa in [ls['autoannotations'] for ls in marshal_data['metadata']]]

        iagn = (
            (metrics["classification"] == "AGN")
            | (metrics["classification"] == "blazar")
            | (metrics["classification"] == "QSO")
            | (metrics["classification"] == "Blazar")
        )  # This would use AGN detection...
        icv = np.array(["CV" in str(x) for x in metrics["classification"]])
        istar = np.array(["st" in str(x) for x in metrics["classification"]]) | icv
        ibogus = np.array(["bogus" in str(x) for x in metrics["classification"]])
        iflexok = (
            (metrics_flex["chi2"] < 10)
            & (metrics_flex["photoclass"] != "no fit")
            & (metrics_flex["photoclass"] != "fast/weird")
        )

        iok = (
            (metrics["peak_mag"] < 20)
            & (metrics["peak_diff_mag"] < 2)
            & (metrics["peak_mag"] > 0)
        )

        ipublic = (
            (metrics["latest_mag"] < 19.5)
            & (metrics["latest_diff_mag"] < 1)
            & (
                metrics["ndetections"]
                > 4 & metrics["lastest_obs_is_detection"].astype(bool)
            )
        )

        inone = metrics["classification"] == "None"  # Marshal classification
        itde = (metrics["classification"] == "TDE") | (
            metrics["classification"] == "TDE?"
        )  # Marshal classification
        imaybe = ["?" in str(x["classification"]) for ind, x in metrics.iterrows()]
        inotsecure = inone | imaybe

        # to do: get this info again from Fritz (and Growth?)
        # ihasspec = np.array([bool(ls['n_spectra']) for ls in marshal_data['metadata']])
        # self.logger.info ('# of sources with spectra: '+  str(sum(ihasspec)))
        # inospec = ihasspec==False

        # try a public stream
        iselect = (
            ipublic
            & iage
            & iflexok
            & (iagn == False)
            & (ibogus == False)
            & (istar == False)
        )
        self.logger.info("# sources in public list: " + str(sum(iselect)))
        filename = base_name + "_public_testing.txt"
        self.print_ranking(filename, metrics[iselect], metrics_flex[iselect])

        # spectrum without secure marshal class
        # iselect =  iok & iflexok & inotsecure # & ihasspec
        # filename = base_name.split('_')[0]+'_spectrum_but_no_classification.txt'
        # self.print_ranking(filename, metrics[iselect], metrics_flex[iselect])

        # everything (including bad fit etc)
        iselect = iage
        self.logger.info("# sources in everything list: " + str(sum(iselect)))
        filename = base_name + "_everything.txt"

        self.print_ranking(filename, metrics[iselect], metrics_flex[iselect])

        # SEDM select bright sources
        # to do: filter out sources wih a spectrum
        ibrightnow = metrics["latest_mag"] < 19.5
        iselect = ibrightnow & iflexok & iage & iok & inotsecure  # & inospec
        self.logger.info("# sources in SEDM list: " + str(sum(iselect)))

        filename = base_name + "_unclassified_SEDM.txt"
        self.print_ranking(filename, metrics[iselect], metrics_flex[iselect])

        # select TDE candidates, plus known TDEs
        iphototde = (metrics_flex["photoclass"] == "TDE") | (
            metrics_flex["photoclass"] == "TDE?"
        )
        iselect = iage & iok & iflexok & iphototde & (inotsecure | itde)
        self.logger.info("# sources in TDE list: " + str(sum(iselect)))

        filename = base_name + "_TDEs.txt"
        self.print_ranking(filename, metrics[iselect], metrics_flex[iselect])

        # select source for rapid Swift follow-up:
        # blue, slow and pre-peak
        ipretde = (
            (metrics_flex["mean_color"] < 0.0)
            & (np.abs(metrics_flex["e_mean_color"]) < 0.35)
            & (metrics_flex["sigma_rise"] > 9)
            & (
                metrics_flex["sigma_rise"]
                / np.clip(metrics_flex["e_sigma_rise"], 1, 1e99)
                > 2
            )
            & (metrics_flex["n_rise"] > 2)
            & (metrics_flex["n_fade"] < 3)
        )
        iselect = iage & iok & iflexok & ipretde & (inotsecure | itde)  # &inospec
        self.logger.info("# sources in pre-peak TDE list: " + str(sum(iselect)))
        filename = base_name + "_prepeak_TDEs.txt"

        self.print_ranking(filename, metrics[iselect], metrics_flex[iselect])

    def get_ranking(self, metrics, metrics_flex):
        """
        attempt at ranking sources for follow-up, using resulting from
        light curve fitting (flexfit) and basic properties  (simplemetrics)

        """

        rank = np.zeros(len(metrics), dtype=float)
        rank[
            metrics["peak_mag"] == 0
        ] = 1000  # if it didnt pass filter, even lower rank
        ipbad = np.array([bool("bad" in str(pc)) for pc in metrics_flex["photoclass"]])
        ipagn = np.array([bool("AGN" in str(pc)) for pc in metrics_flex["photoclass"]])
        ipcv = np.array([bool("CV" in str(pc)) for pc in metrics_flex["photoclass"]])
        iptde = np.array([bool(str(pc) == "TDE") for pc in metrics_flex["photoclass"]])
        iptde_m = np.array(
            [bool(str(pc) == "TDE?") for pc in metrics_flex["photoclass"]]
        )
        ipsn = np.array(
            [bool(str(pc)[0:2] == "SN") for pc in metrics_flex["photoclass"]]
        )
        ipsn_m = np.array(
            [bool(str(pc)[0:3] == "SN?") for pc in metrics_flex["photoclass"]]
        )

        # decrease rank for fainter reference magnitude
        # removed this in Dec 2020 to decrease host bias
        # rank += np.clip(metrics['ref_mag'].values, 17,25)

        # photmetric TDE get higer rank
        rank[iptde] -= 5
        # maybe photmetric TDE gets somewhat higer rank
        rank[iptde_m] -= 4

        # poor fit gets a lower rank
        rank[ipbad] += 1

        # photmetric AGN get lower rank
        rank[ipagn] += 5

        # CV get lower rankp
        rank[ipcv] += 5

        # photometric SN get lower rank
        rank[ipsn] += 5

        # photometric potential SN get somewhat lower rank
        rank[ipsn_m] += 3

        # increase rank value by number of days we didnt see any detections
        rank += np.clip(metrics["age"], 3, 1000) / 10

        # small preference for source with more detections
        rank -= np.clip(metrics["ndetections"], 1, 10) / 10

        # strong anti preference for sources with many negative subtractions
        rank += (
            metrics["ndetections_negative"]
            / np.clip(metrics["ndetections"], 1, 1e99)
            * 10
        )

        # no ranking based on offset for now, because symatics can be large
        # rank += np.clip(metrics['offset_sig'], 2., 20)
        # rank += np.clip(metrics['offset_sig_ps1'], 2., 20)

        # brighter source somewhat higher rank
        rank += np.clip(metrics["peak_mag"] - 19, 0, 3)

        # linear increase with magnitude flare-host mag difference
        rank += np.clip(metrics["latest_diff_mag"], -2, 5) / 2

        # small extra penalty for upper limits in last detection
        # isignon = (metrics['lastest_obs_is_detection']==0) & (metrics['latest_maglim']>20)
        # self.logger.info ('# sources with last obs is nondetection {}'.format(sum(isignon)))
        # rank[isignon] +=0.5

        # number one gets rank=1
        rank -= np.nanmin(rank) - 1

        return rank

    # This print function needs to also show some info about current classification and current followup status
    def print_ranking(self, filename, metrics, metrics_flex):
        """ """

        # header
        ss = "#name           rank   mag     diff     age   neg/pos  off_stat_ztf off_stat_ps1  mean_color  color_slope     rise       fade   photo_class  official_class    more info"

        for i in np.argsort(metrics["ranking"]):

            ss += "\n{0} {1:6.1f}   {2:5.2f}   {3:5.2f}    {4:5.1f}   {5:5.2f}     {6:7.1f}     {7:7.1f}       {8:5.2f}      {9:7.4f}     {10:9.1f}  {11:9.1f}   {12:14}    {13:7}    {14}".format(
                str(metrics_flex.iloc[i]["name"]),
                float(metrics["ranking"].iloc[i]),
                float(metrics.iloc[i]["latest_mag"]),
                float(metrics.iloc[i]["latest_diff_mag"]),
                float(metrics.iloc[i]["age"]),
                float(
                    metrics.iloc[i]["ndetections_negative"]
                    / np.clip(metrics.iloc[i]["ndetections"], 1, 1e99)
                ),
                float(metrics.iloc[i]["offset_sig"]),
                float(metrics.iloc[i]["offset_sig_ps1"]),
                float(metrics_flex.iloc[i]["mean_color"]),
                float(metrics_flex.iloc[i]["color_slope"]),
                np.clip(float(metrics_flex.iloc[i]["sigma_rise"]), 0, 999),
                np.clip(float(metrics_flex.iloc[i]["sigma_fade"]), 0, 999),
                str(metrics_flex.iloc[i]["photoclass"]),
                str(metrics.iloc[i]["classification"]),
                str(metrics.iloc[i]["extra_info"]),
            )

        self.put(filename, bytes(ss + "\n", "utf-8"))

    def make_sumplots(self, metrics, metrics_flex):
        """
        make summary plot
        currently hardcoded to work for ZTFBH Nuclear
        """

        from numpy import log10, sqrt, log

        inotbogus = np.array(
            [
                bool(
                    (rc != "star")
                    & (rc != "bogus")
                    & (rc != "varstar")
                    & (rc != "star?")
                    & (rc != "stellar")
                )
                for rc in metrics["classification"].values
            ]
        )

        results = metrics_flex[inotbogus].copy()
        if len(results) == 0:
            return

        # some cuts on the Fritz/Marshal/TNS/catalog labels (from classifyme)
        iquestion = np.array(
            [
                ("?" in c) and (c != "TDE?")
                for c in results["classification"].astype(str)
            ]
        )
        iagn = (results["classification"] == "AGN") | (
            results["classification"] == "QSO"
        )
        isnIa = np.array(
            [c[0:5] == "SN Ia" for c in results["classification"].astype(str)]
        ) & (iquestion == False)
        isn_all = np.array(
            [c[0:2] == "SN" for c in results["classification"].astype(str)]
        ) & (iquestion == False)
        isn = isn_all & (isnIa == False)
        icv = np.array([c[0:2] == "CV" for c in results["classification"].astype(str)])
        inone = (
            np.array([c == "None" for c in results["classification"].astype(str)])
            | iquestion
        )

        self.logger.info("\n\n----")
        for ii, ttype in zip(
            *([iagn, isnIa, isn, icv, inone], ("AGN", "Ia", "SN other", "CV", "none"))
        ):
            self.logger.info(f"{ttype}: {sum(ii)}")
        self.logger.info("----\n\n")

        # add this back in later:
        itde = np.array(
            [c[0:3] == "TDE" for c in results["classification"].astype(str)]
        )

        # rename for GOT
        # SvV to do: make a version of the summary plot with all the TDEs
        # GoT_dict = self.read_file('/mampel/misc/GOTnames.json')[1].json()

        # for i, source in results.iterrows():
        #   if source['name'] in GoT_dict.keys():
        #       source['name'] = GoT_dict[source['name']]

        # add marshal class to name (for plotting)
        results["name"] = results.apply(
            lambda r: r["name"] + "(" + r["classification"] + ")"
            if (not "(" in r["name"])
            and (r["classification"] != "TDE")
            and (r["classification"] != "None")
            else r["name"],
            axis=1,
        )

        plot_tuple = (
            [itde, iagn, isnIa, isn, isn_all, inone],
            ("TDE", "AGN", "SN Ia", "SN other", "SN all", "unknown"),
            (1, 0.3, 0.5, 0.8, 0.8, 0.5, None, 0.9),
            ("d", ".", "s", "o", None, "x"),
        )
        mew = 0.6
        lw = 0.6
        ms = 2

        # ---

        # ----
        # color maps obsession
        # cmap = mpl.cm.get_cmap('magma')

        cmap = mpl.cm.get_cmap("viridis")
        manual_colors = [cmap(0.1)] + [
            cmap(0.65),
            cmap(0.65),
            cmap(0.65),
            cmap(0.35),
            cmap(0.92),
        ]
        labels = ["unknown", "SN other", "SN all", "SN Ia", "TDE", "AGN"]
        color_dict = dict(zip(labels, manual_colors))

        if self.plotHisto:
            plt.clf()
            ipl = results["mag_peak"] > 0
            plt.title("all sources; N=" + str(sum(ipl)))
            plt.hist(results[ipl]["mag_peak"], range=[17, 22], bins=20)
            plt.xlabel("peak mag")

            plot_fname = self.plot_dir + "/all_maghisto.pdf"
            self.plt2box(plt, plot_fname)

        # ------
        # rise / fade
        finfo = "{0:30} {1:10} {2:10}\n".format("name", "rise_time", "fade_rime")

        ipl = (
            (results["n_rise"] > 2)
            & (results["n_fade"] > 2)
            & (results["e_sigma_fade"] > 0)
            & (results["sigma_rise"] / results["e_sigma_rise"] > 1.5)
            & (results["sigma_rise"] > 1.1)
            & (results["sigma_fade"] / results["e_sigma_fade"] > 1.5)
        )

        if self.plotHisto:
            plt.clf()
            plt.title("rise/fade selection; N=" + str(sum(ipl)))
            plt.hist(results[ipl]["mag_peak"], range=[17, 22], bins=20)
            plt.xlabel("peak mag")

            plot_fname = self.plot_dir + "/rise_fade_maghisto.pdf"
            self.plt2box(plt, plot_fname)

        reference_agn = pd.read_csv(
            io.BytesIO(self.read_file("/mampel/misc/AGN.csv").content)
        )
        reference_sne = pd.read_csv(
            io.BytesIO(self.read_file("/mampel/misc/SN_all.csv").content)
        )
        reference_tde = pd.read_csv(
            io.BytesIO(self.read_file("/mampel/misc/TDE.csv").content)
        )
        plt.clf()
        for idx, labl, alph, fmt in zip(*plot_tuple):

            self.logger.info(f"{labl}, {sum(ipl&idx)}")
            x, y = log10(results[ipl & idx]["sigma_rise"].astype(np.float64)) + log10(
                sqrt(2)
            ), log10(results[ipl & idx]["sigma_fade"].astype(np.float64))

            xerr = (
                1
                / log(10)
                * results[ipl & idx]["e_sigma_rise"].values
                / results[ipl & idx]["sigma_rise"].values
            )
            yerr = (
                1
                / log(10)
                * results[ipl & idx]["e_sigma_fade"].values
                / results[ipl & idx]["sigma_fade"].values
            )
            x, y = x.values, y.values

            if fmt:
                line = plt.errorbar(
                    x,
                    y,
                    xerr=xerr,
                    yerr=yerr,
                    color=color_dict[labl],
                    fmt=fmt,
                    alpha=alph,
                    label=labl,
                    lw=lw,
                    mew=mew,
                    ms=ms,
                    zorder=3 + 1 / (sum(idx) + 1),
                )

            if labl == "TDE":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    # if nm=='NedStark': nm='AT2018zr'
                    finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(nm, x[i], y[i])
                    plt.text(x[i], y[i], nm, fontsize=4, zorder=5)
            elif labl == "unknown":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    if x[i] > 0 or y[i] > 1.0:
                        finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(
                            nm, x[i], y[i]
                        )
                        plt.text(x[i], y[i], nm, fontsize=3, zorder=6)
            elif labl != "AGN":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    if (x[i] > 0 or y[i] > 1.0) or ("ZTF" in nm is False):
                        finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(
                            nm, x[i], y[i]
                        )
                        plt.text(x[i], y[i], nm, fontsize=2, zorder=5)
        ax = plt.gca()
        corner.hist2d(
            log10(reference_agn["sigma_rise"].values) + log10(sqrt(2)),
            log10(reference_agn["sigma_fade"].values),
            bins=[35, 35],
            range=[(0, log10(500)), (0, log10(500))],
            smooth=True,
            levels=([0.66, 0.66, 0.66]),
            ax=ax,
            color=color_dict["AGN"],
            plot_datapoints=False,
            plot_density=False,
            plot_contours=True,
            no_fill_contours=True,
            fill_contours=False,
            contour_kwargs={"zorder": 0},
        )

        corner.hist2d(
            log10(reference_sne["sigma_rise"].values) + log10(sqrt(2)),
            log10(reference_sne["sigma_fade"].values),
            bins=[35, 35],
            range=[(0, log10(500)), (0, log10(500))],
            smooth=True,
            levels=([0.66, 0.66, 0.66]),
            ax=ax,
            color=color_dict["SN all"],
            plot_datapoints=False,
            plot_density=False,
            plot_contours=True,
            no_fill_contours=True,
            fill_contours=False,
            contour_kwargs={"zorder": 0},
        )

        plt.legend(loc=4)
        plt.xlabel(r"rise e-folding time (log day)")
        plt.ylabel(r"fade e-folding time (log day)")

        plt.xlim(0.7, 2.5)
        plt.ylim(1.0, 2.5)

        plot_fname = self.plot_dir + "/rise_fade.pdf"
        self.plt2box(plt, plot_fname)

        # also plot the known TDEs
        plot_fname = plot_fname.replace(".pdf", "_wtdes.pdf")
        line = plt.errorbar(
            log10(reference_tde["sigma_rise"].values) + log10(sqrt(2)),
            log10(reference_tde["sigma_fade"].values),
            color=color_dict["TDE"],
            fmt="d",
            alpha=0.7,
            lw=lw,
            mew=mew,
            ms=ms,
            zorder=1,
        )
        self.plt2box(plt, plot_fname)

        # ------
        # color / cooling
        finfo += "\n\n {0:30} {1:10} {2:10}\n".format(
            "name", "mean_color", "color_change"
        )
        ipl = (
            (results["mean_color"] != 0)
            & (np.abs(results["e_mean_color"]) < 0.35)
            & (results["n_fade"] > 3)
            & (np.abs(results["e_color_slope"]) < 0.03)
        )

        if self.plotHisto:
            plt.clf()
            plt.title("color/cooling selection; N=" + str(sum(ipl)))
            plt.hist(results[ipl]["mag_peak"], range=[17, 22], bins=20)
            plt.xlabel("peak mag")

            plot_fname = self.plot_dir + "/color_change_maghisto.pdf"
            self.plt2box(plt, plot_fname)

        plt.clf()
        for idx, labl, alph, fmt in zip(*plot_tuple):

            self.logger.info(f"{labl}, {sum(ipl&idx)}")
            x, y = results[ipl & idx]["mean_color"], results[ipl & idx]["color_slope"]
            x, y = x.values, y.values  # to np array

            # clip to keep very blue sources in range
            x = np.clip(x, -0.6, 0.9)

            if fmt:
                line = plt.errorbar(
                    x,
                    y,
                    xerr=results[ipl & idx]["e_mean_color"],
                    yerr=results[ipl & idx]["e_color_slope"],
                    color=color_dict[labl],
                    fmt=fmt,
                    alpha=alph,
                    lw=lw,
                    mew=mew,
                    ms=ms,
                    label=labl,
                    zorder=3 + 1 / (sum(idx) + 1),
                )

            if labl == "TDE":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(nm, x[i], y[i])
                    plt.text(x[i] - 0.01, y[i], nm, fontsize=1.5, zorder=5)
            elif (labl != "AGN") and (labl[0:2] != "SN"):
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    if (
                        (y[i] < 0.02)
                        and not ("AGN?" in nm)
                        and not ("CV?" in nm)
                        and not ("varstar?" in nm)
                        or (x[i] <= 0.4)
                    ):

                        if "?" in nm:
                            fs = 1
                        else:
                            fs = 1.5
                            finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(
                                nm, x[i], y[i]
                            )
                            plt.text(x[i], y[i], nm, fontsize=fs, zorder=5)
                        self.logger.info(f"{nm}, {not('AGN?' in nm)}")

        ax = plt.gca()
        corner.hist2d(
            reference_sne["mean_color"].values,
            reference_sne["color_slope"].values,
            bins=[10, 20],
            range=[(-0.4, 0.8), (-0.03, +0.045)],
            smooth=True,
            levels=([0.66, 0.66, 0.66]),
            ax=ax,
            color=color_dict["SN all"],
            plot_datapoints=False,
            plot_density=False,
            plot_contours=True,
            no_fill_contours=True,
            fill_contours=False,
            contour_kwargs={"zorder": 0},
        )
        plt.legend()
        plt.xlabel("mean g-r color")
        plt.ylabel("rate of color change (1/day)")
        plt.xlim(-0.6, 0.8)
        plt.ylim(-0.020, 0.04)
        # plt.draw()
        # plt.pause(0.05)

        # option to add a "TDE box"
        # plt.plot([-1,-0.00],[0.015,0.015], '--k',alpha=0.6, lw=0.5)
        # plt.plot([-0.1,-0.00],[-1,0.015], '--k',alpha=0.6, lw=0.5)
        # txt = plt.text(-0.3, -0.01, 'TDE box',  fontsize='large', verticalalignment='center', horizontalalignment='center')
        # txt.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='w'))

        plt.title("late-time (post-peak) selection")
        plot_fname = self.plot_dir + "/color_change.pdf"
        self.plt2box(plt, plot_fname)

        # also plot the known TDEs
        plot_fname = plot_fname.replace(".pdf", "_wtdes.pdf")
        line = plt.errorbar(
            reference_tde["mean_color"].values,
            reference_tde["color_slope"].values,
            color=color_dict["TDE"],
            fmt="d",
            alpha=0.7,
            lw=lw,
            mew=mew,
            ms=ms,
            zorder=1,
        )
        self.plt2box(plt, plot_fname)

        # ------
        # rise / color
        # print ('\n\n rise/color')
        finfo += "\n\n {0:30} {1:10} {2:10}\n".format("name", "mean_color", "rise_time")
        ipl = (
            (results["mean_color"] != 0)
            & (results["n_rise"] > 2)
            & (results["sigma_rise"] / results["e_sigma_rise"] > 1.5)
            & (abs(results["e_mean_color"]) < 0.35)
        )

        if self.plotHisto:
            plt.clf()
            plt.title("rise-time + color selection selection; N=" + str(sum(ipl)))
            plt.hist(results[ipl]["mag_peak"], range=[17, 22], bins=10)
            plt.xlabel("peak mag")

            plot_fname = self.plot_dir + "/rise_color_maghisto.pdf"
            self.plt2box(plt, plot_fname)

        plt.clf()
        for idx, labl, alph, fmt in zip(*plot_tuple):

            self.logger.info(f"{labl}, {sum(ipl&idx)}")

            x, y = results[ipl & idx]["mean_color"], log10(
                results[ipl & idx]["sigma_rise"]
            ) + log10(sqrt(2))
            yerr = (
                1
                / log(10)
                * results[ipl & idx]["e_sigma_rise"]
                / results[ipl & idx]["sigma_rise"]
            )
            x, y = x.values, y.values

            # clip to keep very blue sources in the plot range
            x = np.clip(x, -0.6, 0.9)

            if fmt:
                line = plt.errorbar(
                    x,
                    y,
                    xerr=results[ipl & idx]["e_mean_color"],
                    yerr=yerr,
                    color=color_dict[labl],
                    fmt=fmt,
                    alpha=alph,
                    label=labl,
                    lw=lw,
                    mew=mew,
                    ms=ms,
                    zorder=3 + 1 / (sum(idx) + 1),
                )

            if labl == "TDE":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(nm, x[i], y[i])
                    plt.text(x[i], y[i], nm, fontsize=2, zorder=5)
            elif labl != "AGN":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    if ((x[i] < 0.3 and y[i] > 0) and not ("SN" in nm)) or (
                        x[i] < -0.1 or y[i] > 0.5
                    ):

                        if ("?" in nm) or ("Ia" in nm):
                            fs = 1
                        else:
                            fs = 1.6
                        finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(
                            nm, x[i], y[i]
                        )
                        plt.text(x[i], y[i], nm, fontsize=fs, zorder=6)

        ax = plt.gca()
        corner.hist2d(
            reference_agn["mean_color"].values,
            log10(reference_agn["sigma_rise"].values) + log10(sqrt(2)),
            bins=[30, 30],
            range=[(-0.6, 0.8), (0, log10(300))],
            smooth=True,
            levels=([0.66, 0.66, 0.66]),
            ax=ax,
            color=color_dict["AGN"],
            plot_datapoints=False,
            plot_density=False,
            plot_contours=True,
            no_fill_contours=True,
            fill_contours=False,
            contour_kwargs={"zorder": 0},
        )  # 'lw':1 keyword not used by corner
        corner.hist2d(
            reference_sne["mean_color"].values,
            log10(reference_sne["sigma_rise"].values) + log10(sqrt(2)),
            bins=[30, 30],
            range=[(-0.6, 0.8), (0, log10(300))],
            smooth=True,
            levels=([0.66, 0.66, 0.66]),
            ax=ax,
            color=color_dict["SN all"],
            plot_datapoints=False,
            plot_density=False,
            plot_contours=True,
            no_fill_contours=True,
            fill_contours=False,
            contour_kwargs={"zorder": 0},
        )

        plt.title("real-time (pre-peak) selection")
        plt.legend(loc=1)
        plt.xlabel("mean g-r color")
        plt.ylabel(r"rise e-folding time (log$_{10}$ day)")
        plt.xlim(-0.6, 0.8)
        plt.ylim(log10(4), log10(300))

        plot_fname = self.plot_dir + "/rise_color.pdf"
        self.plt2box(plt, plot_fname)

        # also plot the known TDEs
        plot_fname = plot_fname.replace(".pdf", "_wtdes.pdf")
        line = plt.errorbar(
            reference_tde["mean_color"].values,
            log10(reference_tde["sigma_rise"].values) + log10(sqrt(2)),
            color=color_dict["TDE"],
            fmt="d",
            alpha=0.7,
            lw=lw,
            mew=mew,
            ms=ms,
            zorder=1,
        )
        self.plt2box(plt, plot_fname)

        # as a bonus: write a file with the names of TDE candidates (based on some boxing)
        # for each of the plots (requested by Suvi)
        self.put(self.plot_dir + "/boxed.txt", bytes(finfo, "utf-8"))

    # pyplot to dropbox function
    def plt2box(self, plt, fname):
        """
        write PDF of plot in pyplot to dropbox
        (could be more neat to work with ax objects everywhere, but hey)
        >> plt2box(plt, plot_name)
        """
        buf = io.BytesIO()
        ax = plt.gca()
        ax.figure.savefig(buf, format="pdf")
        buf.seek(0)
        self.put(fname, buf.read())

    # some last helper functions
    def datetomjd(self, d):
        d0 = datetime.datetime(1858, 11, 17, 0, 0, 0)
        dt = d - d0
        # dt is a timedelta object.
        return dt.days + (dt.seconds + dt.microseconds / 1e6) / 86400.0

    def mjdtodate(self, mjd):
        if not np.isfinite(mjd):
            return None
        jd = mjd + 2400000.5
        unixtime = (jd - 2440587.5) * 86400.0  # in seconds
        return datetime.datetime.utcfromtimestamp(unixtime)

    def format_date(self, transient_date):
        return astropy.time.Time(
            transient_date[6:10]
            + "-"
            + transient_date[3:5]
            + "-"
            + transient_date[0:2]
            + transient_date[10:]
        ).jd
