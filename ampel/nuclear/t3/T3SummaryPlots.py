#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t3/T3MetricsPlots.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu
# Date              : 08.06.2020
# Last Modified Date: 29.11.2022
# Last Modified By  : simeon.reusch@desy.de


import os, sys, datetime, io, json, warnings
from collections.abc import Generator
from typing import Any, Optional, Union

from pydantic import BaseModel
import matplotlib as mpl  # type: ignore
import corner  # type: ignore
import astropy  # type: ignore
from numpy import log10, sqrt, log
import numpy as np
import pandas as pd  # type: ignore
from matplotlib import pyplot as plt

from ampel.types import UBson, T3Send
from ampel.view.TransientView import TransientView
from ampel.view.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.nuclear.t3.dropboxIO import DropboxUnit
import ampel.nuclear.t3.classifyme as classifyme

warnings.filterwarnings("ignore", category=RuntimeWarning)

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
    "e_sigma_fade",
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


class T3SummaryPlots(DropboxUnit, AbsPhotoT3Unit):

    plotHisto: bool = False  # create histrograms for of the magnitude of the transients in the different summary plots
    force_date: Optional[str] = None  # force the save date for testing

    def post_init(self):
        super().post_init()
        self.save_location = self.base_location + "/sum_plots"
        if "sum_plots" not in self.get_files(self.base_location):
            self.create_folder(self.save_location)

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

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        """ """
        # DUMMY FUNCTION FOR NOW
        return None

    def add(self, transients):
        """
        print some detection stats
        and make summary plost
        currently hardcoded to work for ZTFBH Nuclear
        date_str='yyyy-mm-dd', optional input to selected what date to print
        days_back=60, days since the last detection  to be included in the plot
        """

        drop_dir = self.save_location + "/" + str(self.night.datetime.year)
        if str(self.night.datetime.year) not in self.get_files(self.save_location):
            self.create_folder(drop_dir)

        plot_dir = drop_dir + "/" + self.night.datetime.strftime("%m-%d")
        if self.night.datetime.strftime("%m-%d") not in self.get_files(drop_dir):
            self.create_folder(plot_dir)

        metrics, metrics_flex = self.collect_metrics(transients)

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

        mjdmay = self.datetomjd(datetime.datetime(2018, 5, 1))
        mjdnov = self.datetomjd(datetime.datetime(2018, 10, 1))
        jdnow = self.datetomjd(self.night.datetime) + 2400000.5
        mjdmax = np.nanmax(metrics_flex["mjd_peak"])

        self.logger.info(
            "date of max mjd_peak from fit: " + str(self.mjdtodate(mjdmax))
        )

        number_of_days = mjdmax - mjdnov
        ipostnov = (
            metrics_flex["mjd_peak"].values > mjdnov
        )  # used for printing detection stat
        irecent = metrics["latest_jd"].values > -999  # select all, deprecated daysback

        results = metrics_flex.iloc[ipostnov & inotbogus].copy()
        self.logger.info("# of days since 2018-Nov-1: " + str(number_of_days))

        mag_lim = np.arange(20, 18, -0.5)
        iokrise = results["n_rise"] > 2

        self.logger.info("\n>2 pre-peak detections & not AGN:")

        for magl in mag_lim:
            ilim = results["mag_peak"] < magl
            rr = results[iokrise & ilim].copy()
            inotagn = (rr["photoclass"] != "AGN?") & (rr["classification"] != "AGN")
            self.logger.info(
                "{0:0.1f} {1:0.1f}/week".format(magl, sum(inotagn) / number_of_days * 7)
            )

        self.logger.info("\n>2 pre-peak detections & not AGN & not photo SNe:")

        for magl in mag_lim:
            ilim = results["mag_peak"] < magl
            rr = results[iokrise & ilim].copy()
            inotagn = (rr["photoclass"] != "AGN?") & (rr["classification"] != "AGN")
            inotsn = rr["photoclass"] == "not Ia?"
            self.logger.info(
                "{0:0.1f} {1:0.1f}/week".format(
                    magl, sum(inotsn & inotagn) / number_of_days * 7
                )
            )

        self.logger.info("\n\n>3 post-peak detections & not AGN:")
        iokfade = results["n_fade"] > 3
        for magl in mag_lim:
            ilim = results["mag_peak"] < magl
            rr = results[iokfade & ilim].copy()
            inotagn = (rr["photoclass"] != "AGN?") & (rr["classification"] != "AGN")
            self.logger.info(
                "{0:0.1f} {1:0.1f}/week".format(magl, sum(inotagn) / number_of_days * 7)
            )

        self.logger.info(
            "\n\n>3 post-peak detections & TDE candidate (blue or not cooling):"
        )
        iokfade = results["n_fade"] > 3
        for magl in mag_lim:
            ilim = results["mag_peak"] < magl
            rr = results[iokfade & ilim].copy()
            itde = [bool("TDE" in pc) for pc in rr["photoclass"] if "bad" not in pc]
            self.logger.info(
                "{0:0.1f} {1:0.1f}/week".format(magl, sum(itde) / number_of_days * 7)
            )

        self.logger.info("\n>10 post-peak detections & TDE candidate & not red")
        iokfade = results["n_fade"] > 3
        for magl in mag_lim:
            ilim = results["mag_peak"] < magl
            rr = results[iokfade & ilim].copy()
            itde = [
                "TDE" in pc
                for pc in rr["photoclass"]
                if "red" not in pc and "bad" not in pc
            ]
            self.logger.info(
                "{0:0.1f} {1:0.1f}/week".format(magl, sum(itde) / number_of_days * 7)
            )

        results = metrics_flex[irecent & inotbogus].copy()
        if len(results) == 0:
            return

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
        GoT_dict = self.read_file("/mampel/misc/GOTnames.json")[1].json()

        for i, source in results.iterrows():
            if source["name"] in GoT_dict.keys():
                source["name"] = GoT_dict[source["name"]]

        # results[results['name']=='ZTF18ablllyw'] = 'ZTF18ablllyw (SN?)' # SN with pre-cursor?!
        # results[results['name']=='ZTF18aayatjf'] = 'ZTF18aayatjf (SN?)' # slow fading SN

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

            plot_fname = plot_dir + "/all_maghisto.pdf"
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

            plot_fname = plot_dir + "/rise_fade_maghisto.pdf"
            self.plt2box(plt, plot_fname)

        reference_agn = pd.read_csv(
            io.BytesIO(self.read_file("/mampel/misc/AGN.csv")[1].content)
        )
        reference_sne = pd.read_csv(
            io.BytesIO(self.read_file("/mampel/misc/SN_all.csv")[1].content)
        )
        reference_tde = pd.read_csv(
            io.BytesIO(self.read_file("/mampel/misc/TDE.csv")[1].content)
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
                    if x[i] > 1.2 or y[i] > 1.6:
                        finfo += "{0:26} {1:10.3f}     {2:10.3f}\n".format(
                            nm, x[i], y[i]
                        )
                        plt.text(x[i], y[i], nm, fontsize=3, zorder=6)
            elif labl != "AGN":
                for i, nm in enumerate(results[ipl & idx]["name"]):
                    if (x[i] > 1.3 or y[i] > 1.5) or ("ZTF" in nm is False):
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

        plot_fname = plot_dir + "/rise_fade.pdf"
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

            plot_fname = plot_dir + "/color_change_maghisto.pdf"
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
        plot_fname = plot_dir + "/color_change.pdf"
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

            plot_fname = plot_dir + "/rise_color_maghisto.pdf"
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

        plot_fname = plot_dir + "/rise_color.pdf"
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
        self.put(plot_dir + "/boxed.txt", bytes(finfo, "utf-8"))

        self.commit()

    def collect_metrics(self, transients):
        """ """
        jd_today = self.night.jd

        alert_name = "/mampel/alerts"
        alert_folders = self.get_files(alert_name)
        if "ps1_offset.json" not in alert_folders:
            buf = io.BytesIO()
            buf.write(bytes(json.dumps({}), "utf-8"))
            buf.seek(0)
            self.put(alert_name + "/ps1_offset.json", buf.read())
            ps1_offset = {}
        else:
            ps1_offset = self.read_file(alert_name + "/ps1_offset.json")[1].json()

        if transients is not None:

            local_sources = np.array(transients)  # deprecated daysback soln
            newdetection_count = len(local_sources)
            self.logger.info(
                "Collecting metrics for {} transient(s)".format(newdetection_count)
            )
            metrics, metrics_flex = [], []
            for i, tran_view in enumerate(local_sources):

                # get classifcation plus some neoWISE info
                classification, extra_info = classifyme.get_class(
                    tran_view, self, self.logger, write=True
                )

                metrics_flex.append(self.get_fit_params(tran_view, classification))
                metrics.append(self.get_simple_results(tran_view, classification))

            metrics_flex = pd.DataFrame(metrics_flex, columns=flex_keys)
            metrics = pd.DataFrame(metrics, columns=simple_keys)

            self.logger.info(
                "successfully collected metrics: "
                + str(metrics_flex["name"].all() == metrics["name"].all())
            )
            if metrics_flex["name"].all() != metrics["name"].all():
                warnings.warn(
                    "WARNING: flex metrics and simple metrics are not aligned, plots may be erroneous."
                )
        else:
            self.logger.info("no transients to rank!")

        return metrics, metrics_flex

    def get_fit_params(self, tran_view, classification):
        try:
            fit_params = (tran_view.get_t2_result("T2FlexFit") or {})["fit_params"]
            return {**fit_params, **{"classification": classification}}

        except KeyError:
            self.logger.warn(
                "WARNING: No Flex fit for source {}".format(tran_view.stock["name"][0])
            )
            return {
                "name": tran_view.stock["name"][0],
                "classification": classification,
            }

    def get_simple_results(self, tran_view, classification):
        try:
            simple_results = (tran_view.get_t2_result("T2SimpleMetrics") or {})[
                "metrics"
            ]
            return {**simple_results, **{"classification": classification}}

        except KeyError:
            self.logger.warn(
                "WARNING: No simple metrics for source {}".format(
                    tran_view.stock["name"][0]
                )
            )
            return {
                "name": tran_view.stock["name"][0],
                "classification": classification,
            }

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
