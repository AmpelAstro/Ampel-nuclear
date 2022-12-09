#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t3/T3MetricsPlots.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu
# Date              : 08.06.2020
# Last Modified Date: 29.11.2022
# Last Modified By  : simeon.reusch@desy.de

import io, time, json

from typing import Tuple, Dict, List, Any, Optional, Union
from collections.abc import Generator

import numpy as np
from pydantic import BaseModel
import astropy  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import tempfile

from ampel.types import UBson, T3Send
from ampel.view.TransientView import TransientView
from ampel.view.T3Store import T3Store
from ampel.view.LightCurve import LightCurve
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.nuclear.flexfit import flexfit
from ampel.nuclear.t3.classifyme import jsonify
from ampel.nuclear.t2.T2FlexFit import get_raw_lc
from ampel.nuclear.t3.dropboxIO import DropboxUnit

dump_keys = [
    "name",
    "ampelID",
    "dropbox_path",  # :"alerts/2019/"" or False if no path
    "ra",
    "dec",
    "magpsf_r",
    "magpsf_g",
    "magnr_r",
    "magnr_g",
    "flexclass",
    "neoWISE_chi2_w1",
    "neoWISE_chi2_w2",
    "allWISE_w1-w2",
    "neoWISE_w1-w2",
    "catalogs_matched",  # ["PS1", "Gaia", ...]
    "Gaia_compactness",
]


def get_t2_result_with_compound(
    view: TransientView, unit_id: str
) -> Tuple[Optional[Dict[str, Any]], Optional[Union[int, bytes, str]]]:

    t2_res = reversed(list(view.get_t2_views(unit=unit_id) or []))

    for t2 in t2_res:
        result = t2.get_payload()

        if result and t2.code >= 0:
            assert isinstance(result, dict)
            return result, t2.link

    return None, None


def get_lightcurve(
    view: TransientView, compound_id: Union[int, bytes, str]
) -> Optional[LightCurve]:
    for lc in reversed(view.get_lightcurves() or []):
        if lc.compound_id == compound_id:
            return lc
    return None


class T3MetricsPlots(DropboxUnit):

    verbose: bool = False  #: give feedback on what is plotted

    def post_init(self):
        """
        Run config keys: save_location, verbose, mode
        """
        # super().post_init()
        self.save_location = self.base_location + "/alerts"
        if "alerts" not in self.get_files(self.base_location):
            self.create_folder(self.save_location)
        self.dump_log = []

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        """
        Loops through transients and plots results of fit, offset plots, and dumps fit info, simple metrics, and raw lightcurve
        """

        all_transients = [t for t in gen]

        if all_transients is not None:
            # select sources with new data
            good_sources = np.array(all_transients)

            tran_years = np.array(
                ["20" + tran.stock["name"][0][3:5] for tran in good_sources]
            )
            year_folders = self.get_files(
                self.save_location
            )  # this should also include neowise log and ps1 log
            sources_years = (
                {}
            )  # dict: {year1: sources first detected in year1, year2:, '' year2, ...}
            for year in np.unique(tran_years):  # create folders for all relevant years
                if str(year) not in year_folders:
                    self.create_folder(self.save_location + f"/{year}")
                thisyear_trans = good_sources[tran_years == year]
                self.logger.info(
                    "Adding to {} sources in year {}".format(
                        str(len(thisyear_trans)), year
                    )
                )
                sources_years[str(year)] = thisyear_trans

            for (
                year,
                sources,
            ) in sources_years.items():  # loop through and create missing folders
                year_path = self.save_location + f"/{year}"
                entries = self.get_files(year_path)
                for i, tran_view in enumerate(sources):
                    tran_name = tran_view.stock["name"][0]
                    self.dump_log.append({})
                    self.dump_log[-1]["name"] = tran_name
                    if tran_name not in entries:
                        self.create_folder(year_path + f"/{tran_name}")

                for i, tran_view in enumerate(sources):
                    tran_name = tran_view.stock["name"][0]
                    filepath = year_path + f"/{tran_name}/{tran_name}"
                    self.dump_log[-1]["dropbox_path"] = filepath

                    flex, compound_id = get_t2_result_with_compound(
                        view=tran_view, unit_id="T2FlexFit"
                    )
                    simple = [
                        dv.get_payload()
                        for dv in tran_view.get_t2_views(unit="T2SimpleMetrics")
                    ][0]

                    ps1_match_dvs = tran_view.get_t2_views(unit="T2CatalogMatch") or {}
                    ps1_match = [dv.get_value("PS1", dict) for dv in ps1_match_dvs][0]

                    if (
                        flex is not None
                        and compound_id is not None
                        and (lc := get_lightcurve(tran_view, compound_id)) is not None
                    ):
                        dict_lc = get_raw_lc(lc)
                        dump_fname = f"{filepath}_lc.json"
                        buf = io.BytesIO()
                        buf.write(
                            bytes(json.dumps(jsonify(dict_lc), indent=3), "utf-8")
                        )
                        buf.seek(0)
                        self.put(dump_fname, buf.read())

                    if flex and "fit_params" in flex.keys():
                        if flex["fit_params"]["photoclass"] != "no fit":
                            self.fit_plot(flex, filepath)
                            if self.verbose:
                                self.logger.info(
                                    "successfully plotted flex fit for source {}".format(
                                        tran_name
                                    )
                                )
                            dump_fname = f"{filepath}_flex.json"
                            buf = io.BytesIO()
                            buf.write(
                                bytes(json.dumps(flex["fit_params"], indent=3), "utf-8")
                            )
                            buf.seek(0)
                            self.put(dump_fname, buf.read())

                        else:
                            self.logger.info(
                                "no flex fit for source {}".format(tran_name)
                            )
                    else:
                        self.logger.warn(
                            "empty flex fit for source {}".format(tran_name)
                        )

                    if ps1_match and simple:
                        if "metrics" in simple.keys():
                            (
                                dra_ps1,
                                ddec_ps1,
                                offset_med_ps1,
                                offset_rms_ps1,
                                offset_sig_ps1,
                            ) = self.get_ps1(simple, ps1_match, filepath)
                            simple = dict(simple)
                            simple["plot_info"] = {
                                "dra_ps1": dra_ps1,
                                "ddec_ps1": ddec_ps1,
                                **simple["plot_info"],
                            }
                            simple["metrics"] = {
                                "offset_med_ps1": offset_med_ps1,
                                "offset_rms_ps1": offset_rms_ps1,
                                "offset_sig_ps1": offset_sig_ps1,
                                **simple["metrics"],
                            }

                    if simple:
                        if "metrics" in simple.keys():
                            self.offset_plots(simple, filepath)
                            if self.verbose:
                                self.logger.info(
                                    "successfully plotted offset plots for source {}".format(
                                        tran_name
                                    )
                                )

                            dump_fname = f"{filepath}_simple.json"
                            buf = io.BytesIO()
                            buf.write(
                                bytes(json.dumps(simple["metrics"], indent=3), "utf-8")
                            )
                            buf.seek(0)
                            self.put(dump_fname, buf.read())

                        else:
                            self.logger.info(
                                "no simple metrics for source {}".format(tran_name)
                            )
                    else:
                        self.logger.warn(
                            "empty simple metrics for source {}".format(tran_name)
                        )

        self.commit()
        return None

    def get_ps1(self, simple, ps1_match, filepath):
        metrics, plot_info = simple["metrics"], simple["plot_info"]
        ps1_ra, ps1_dec = ps1_match["RA"], ps1_match["Dec"]
        dra_ps1 = (np.array(metrics["ra"]) - ps1_ra) * 3600  # mean dra
        ddec_ps1 = (np.array(metrics["dec"]) - ps1_dec) * 3600  # mean dec
        dstd_ps1 = np.clip(
            np.std(np.append(dra_ps1, ddec_ps1)) / np.sqrt(len(dra_ps1)), 0.1, 1
        )  # double-check this!
        offset_med_ps1 = np.sqrt(np.median(dra_ps1) ** 2 + np.median(ddec_ps1) ** 2)
        offset_rms_ps1 = np.std(np.append(dra_ps1, ddec_ps1))
        offset_sig_ps1 = offset_med_ps1 / dstd_ps1

        # self.logger.info('simple_metrics: distance to PS1 (median, rms, significance) : {0:0.3f} {1:0.3f} {2:0.3f}'.format(
        #   offset_med_ps1, offset_rms_ps1, offset_sig_ps1))

        return dra_ps1, ddec_ps1, offset_med_ps1, offset_rms_ps1, offset_sig_ps1

    def fit_plot(self, result, filepath):
        fit_params, lc_data, plot_info = (
            result["fit_params"],
            result["lightcurve_data"],
            result["plot_info"],
        )

        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)

        filters = ["g", "r"]
        ax.clear()

        if plot_info:

            dtime_out = plot_info["dtime_out"]
            if (dtime_out != None) and (plot_info["flux_out"] != None):
                ax.plot(dtime_out, plot_info["flux_out"], "sk", zorder=1e99)

            dtime, flux, flux_err, fid, isdetect = (
                np.array(lc_data["dtime"]),
                np.array(lc_data["flux"]),
                np.array(lc_data["flux_err"]),
                np.array(lc_data["fid"]),
                np.array(lc_data["isdetect"]),
            )
            lsq_out = plot_info["lsq_out"]
            dt_fix = plot_info["dt_fix"]
            photoclass = fit_params["photoclass"]
            name = fit_params["name"]
            try:
                classification = fit_params["classification"]
            except KeyError:
                classification = photoclass

            xx = np.linspace(min(dtime), max(dtime), 500)
            ax.set_ylim(0, max(flux[isdetect]) * 1.2)
            ax.set_xlim(min(xx), max(xx))
            ax.set_xlabel("Days since first detection")
            ax.set_ylabel("Flux")

            islim = isdetect == False
            for k in filters:
                try:  # handling missing color fits w/ exceptions, could also do this in T2
                    iflt = fid == k
                    ll = ax.errorbar(
                        dtime[isdetect * iflt],
                        flux[isdetect * iflt],
                        flux_err[isdetect * iflt],
                        fmt="s",
                        color=k,
                    )
                    ax.errorbar(
                        dtime[islim * iflt],
                        flux[islim * iflt],
                        flux_err[islim * iflt],
                        alpha=0.3,
                        fmt="v",
                        color=ll[0].get_color(),
                    )
                    ax.plot(
                        xx,
                        flexfit.broken_gauss(lsq_out[k][0], xx, dt_fix=dt_fix[k]),
                        "--",
                        color=k,
                        lw=1.3,
                        alpha=0.7,
                    )
                    ax.set_ylim(0, max(flux[isdetect] * 1.2))
                except:
                    self.logger.warn("no fit for {} band for source {}".format(k, name))

            xx_0 = np.linspace(min(dtime), max(dtime), sum(fid == filters[0]) * 10)
            xx_1 = np.linspace(min(dtime), max(dtime), sum(fid == filters[1]) * 10)
            xx = np.append(xx_0, xx_1)
            fid_in = np.append(
                np.repeat(filters[0], len(xx_0)), np.repeat(filters[1], len(xx_1))
            )

            # x_peakc arguement needed to make sure we use the time normaliation for the color function
            if "color" in lsq_out.keys():
                model_out = flexfit.broken_gauss_twocolor(
                    lsq_out["color"][0],
                    xx,
                    dt_fix=dt_fix["color"],
                    fid=fid_in,
                    x_peakc=np.median(dtime),
                )
                for flt in filters:
                    ax.plot(
                        xx[fid_in == flt],
                        model_out[fid_in == flt],
                        "-",
                        color=flt,
                        lw=2.0,
                        alpha=0.7,
                    )
            else:
                self.logger.warn("no color calculated for source {}".format(name))

            txt = ax.text(
                0.03,
                0.92,
                plot_info["finfo"],
                transform=ax.transAxes,
                verticalalignment="center",
                fontsize=11,
            )
            txt.set_bbox(dict(facecolor="white", alpha=0.5, edgecolor="k"))

            ax.set_xlim(-8, max(np.array(lc_data["jd"]) - np.array(lc_data["jd0"])) + 2)
            ax.set_title(
                fit_params["name"]
                + " ("
                + str(classification).strip()
                + " - "
                + photoclass
                + ")"
            )
        else:
            txt = ax.text(
                0.03,
                0.92,
                "Fewer than 3 observations for source",
                transform=ax.transAxes,
                verticalalignment="center",
            )
            txt.set_bbox(dict(facecolor="white", alpha=0.5, edgecolor="k"))

            try:
                classification = fit_params["classification"]
            except KeyError:
                classification = fit_params["photoclass"]

            ax.set_title(
                fit_params["name"]
                + " (no classification - "
                + fit_params["photoclass"]
                + ")"
            )

        plot_fname = f"{filepath}_flex.pdf"
        buf = io.BytesIO()
        ax.figure.savefig(buf, format="pdf")
        buf.seek(0)
        self.put(plot_fname, buf.read())
        ax.clear()

    def offset_plots(self, simple, filepath):
        plt.close()
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        metrics, plot_info = simple["metrics"], simple["plot_info"]

        if plot_info["offset_plot"]:
            dra, ddec, snr_weight = (
                np.array(plot_info["dra"]),
                np.array(plot_info["ddec"]),
                np.array(plot_info["snr_weight"]),
            )
            ax.plot(
                dra, ddec, "s", markerfacecolor="none", markeredgecolor="k", label=""
            )

            for band in "gr":
                if f"offset_weighted_{band}" in metrics:
                    mask = np.asarray(plot_info[f"i{band}"])
                    weighted_dra = self.wmean(dra[mask], snr_weight[mask])
                    weighted_ddec = self.wmean(ddec[mask], snr_weight[mask])
                    ax.plot(
                        dra[mask],
                        ddec[mask],
                        "s",
                        markerfacecolor=band,
                        markeredgecolor="k",
                        label="IPAC {0} (mean={1:0.3f}, rms={2:0.3f}, sig={3:0.2f})".format(
                            band,
                            metrics[f"offset_weighted_{band}"],
                            metrics[f"offset_rms_{band}"],
                            metrics[f"offset_sig_{band}"],
                        ),
                    )
                    ax.plot(
                        weighted_dra, weighted_ddec, f"*{band}", markeredgecolor="k"
                    )

            if "offset_med_ps1" in metrics:
                dra_ps1, ddec_ps1 = np.array(plot_info["dra_ps1"]), np.array(
                    plot_info["ddec_ps1"]
                )
                ax.plot(
                    dra_ps1,
                    ddec_ps1,
                    "o",
                    markerfacecolor="orange",
                    markeredgecolor="k",
                    label="PS1     (med ={0:0.3f}, rms={1:0.3f}, sig={2:0.2f})".format(
                        metrics["offset_med_ps1"],
                        metrics["offset_rms_ps1"],
                        metrics["offset_sig_ps1"],
                    ),
                )
                ax.plot(
                    np.median(dra_ps1),
                    np.median(ddec_ps1),
                    "*",
                    markerfacecolor="orange",
                    markeredgecolor="k",
                )

            ax.set_xlabel("RA - <RA>")
            ax.set_ylabel("Dec - <Dec>")
            ax.set_xlim(-1.3, 1.3)
            ax.set_ylim(-1.3, 1.3)
            ax.legend(loc=0)

            name = metrics["name"]
            # ax.set_title(name+' ('+str(self.classification).strip()+')') #to-do: when I have marshal
            ax.set_title(name)

            plot_fname = f"{filepath}_offset.pdf"
            buf = io.BytesIO()
            ax.figure.savefig(buf, format="pdf")
            buf.seek(0)
            self.put(plot_fname, buf.read())
            ax.clear()

    def wmean(self, d, ivar):
        """
        >> mean = wmean(d, ivar)
        inverse-variance weighted mean
        """
        d = np.array(d)
        ivar = np.array(ivar)
        return np.sum(d * ivar) / np.sum(ivar)

    def format_date(self, transient_date):
        return astropy.time.Time(
            transient_date[6:10]
            + "-"
            + transient_date[3:5]
            + "-"
            + transient_date[0:2]
            + transient_date[10:]
        ).jd
