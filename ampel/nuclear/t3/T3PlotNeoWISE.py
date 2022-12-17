#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t3/T3PlotNeoWISE.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu
# Date              : 08.06.2020
# Last Modified Date: 29.11.2022
# Last Modified By  : simeon.reusch@desy.de

import io, json, warnings
from collections.abc import Generator
from typing import Any, Optional, Union

import astropy.io.ascii  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from pydantic import BaseModel
import requests
import tempfile

from ampel.types import UBson, T3Send
from ampel.view.TransientView import TransientView
from ampel.view.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.nuclear.t3.dropboxIO import DropboxUnit
from ampel.secret.NamedSecret import NamedSecret

warnings.filterwarnings("ignore", category=RuntimeWarning)


class T3PlotNeoWISE(DropboxUnit):

    # dropbox_token: None | NamedSecret[str]
    apply_qcuts: bool  #: whether to apply quality cuts to data based on n detections
    plot_allWISE: bool  #: whether or not to plot allWISE data
    verbose: bool = False

    def post_init(self):
        # super().post_init()
        self.save_location = self.base_location + "/alerts"
        if "alerts" not in self.get_files(self.base_location):
            self.create_folder(self.save_location)

    def jsonify(self, obj):
        """
        Recursively cast to JSON native types
        """
        if hasattr(obj, "tolist"):
            return obj.tolist()
        elif hasattr(obj, "keys"):
            return {self.jsonify(k): self.jsonify(obj[k]) for k in obj.keys()}
        elif hasattr(obj, "__len__") and not isinstance(obj, str):
            return [self.jsonify(v) for v in obj]
        else:
            return obj

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        """
        downloads for neoWISE data at IPAC
        this function is based on plot_neoWISE.py from the sZTF repo by sjoertvv.

        plot neo-wise photometry to alert dropbox: {name}-neoWISE.pdf
        back up raw neoWISE data for future matches: {name}-neoWISE.txt
        save in json: info_list, wise_class, out_dict: {name}_neoWISE.json

        time axis is normalized to MJD=58119, Jan-01-2018.

        in the future, we could make this more stand-alone to allow usage by other Ampel teams
        however the need to download the data makes it hard to make this portable.
        """
        selcols_str = "selcols=ra,dec,sigra,sigdec,mjd,w1mpro,w1sigmpro,w1flux,w1sigflux,w1sky,w2flux,w2sigflux,w2mpro,w2sigmpro,w2sky,cc_flags,nb,saa_sep,qual_frame,w1mpro_allwise,w1sigmpro_allwise,w2mpro_allwise,w2sigmpro_allwise,w3mpro_allwise,w3sigmpro_allwise"

        url_to_fill = (
            "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=neowiser_p1bs_psd&objstr={0}&radius=2.5&spatial=Cone&radunits=arcsec&"
            + selcols_str
        )
        t0 = 58119

        all_transients = list(gen)
        all_transient_names = [
            str(t.stock["name"][0]) for t in all_transients if t.stock is not None
        ]

        if all_transient_names is not None:
            try:
                # neowise_log.json tracks which sources we have local data for.  to-do: update log for sufficiently old entries
                data_log = self.read_file(self.save_location + "/neowise_log.json").json()
                self.logger.info(
                    "found {} out of {} sources in data log".format(
                        sum(
                            [
                                True if name in data_log.keys() else False
                                for name in all_transient_names
                            ]
                        ),
                        len(all_transient_names),
                    )
                )
            except Exception as e:
                self.logger.warn(str(e))
                self.logger.warn(
                    "New or no data log (check exception above), downloading all sources from ipac"
                )
                data_log = {}
                buf = io.BytesIO()
                buf.write(bytes(json.dumps(data_log, indent=3), "utf-8"))
                buf.seek(0)
                self.put(self.save_location + "/neowise_log.json", buf.read())
                self.logger.info("created new (blank) neowise_log")

            tran_years = np.array(["20" + name[3:5] for name in all_transient_names])

            year_folders = self.get_files(
                self.save_location
            )  # this should also include neowise log and ps1 log
            sources_years = (
                {}
            )  # dict: {year1: sources first detected in year1, year2:, '' year2, ...}
            for year in np.unique(tran_years):  # create folders for all relevant years
                if str(year) not in year_folders:
                    self.create_folder(self.save_location + f"/{year}")
                thisyear_trans = np.array(all_transients)[tran_years == year]

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
                    print(tran_view.stock["name"])
                    tran_name = tran_view.stock["name"][0]
                    if tran_name not in entries:
                        self.create_folder(year_path + f"/{tran_name}")

            entries = {}  # track entries to know which folder to create
            years = self.get_files(self.save_location)

            # if "neowise_log.json" in years:
            #     years.remove("neowise_log.json")

            for entry in years:
                if entry in ["neowise_log.json", ".DS_Store"]:
                    years.remove(entry)

            for year in years:
                try:
                    entries[year] = self.get_files(
                        self.save_location + "/{}".format(year)
                    )  # [entry.name for entry in self.dbx.files_list_folder( self.save_location + f'/{year}').entries]
                except NotADirectoryError:
                    ...

            for i, tran_view in enumerate(all_transients):
                tran_name = tran_view.stock["name"][0]
                tran_year = "20" + tran_name[3:5]
                need_data = True

                ra, dec = tran_view.get_lightcurves()[-1].get_values(
                    "ra"
                ), tran_view.get_lightcurves()[-1].get_values("dec")
                filebase = self.save_location + f"/{tran_year}/{tran_name}/{tran_name}"
                out_dict: Optional[dict] = {}

                needs_plot = not self.exists(
                    f"{self.save_location}/{tran_year}/{tran_name}/{filebase}-neowise.pdf"
                )

                if (tran_name in data_log.keys()) and (needs_plot):
                    need_data = False
                    self.logger.info(f"reading: {tran_name}_neoWISE.txt")
                    try:
                        astro_tab = astropy.io.ascii.read(
                            self.read_file(filebase + "_neoWISE.txt")[1].text
                        )
                    except Exception as e:  # This occurs for anachronisms, usually due to testing.
                        self.logger.warn(str(e))
                        self.logger.warn(
                            f"Could not find neoWISE log for source {tran_name}, this should not happen.  Check to make sure the path is correct. Re-downloading..."
                        )
                        need_data = True
                if need_data:
                    if tran_name not in entries[tran_year]:
                        self.create_folder(
                            self.save_location + f"/{tran_year}/{tran_name}"
                        )
                        self.logger.info(f"Creating folder for transient {tran_name}")

                    obs_str = "{0:0.6f}+{1:+0.6f}".format(np.median(ra), np.median(dec))

                    url = url_to_fill.format(obs_str)
                    if self.verbose:
                        self.logger.info(f"getting neoWISE data from: \n{url}")
                    try:
                        r = requests.get(url, timeout=120)

                        content = str(r.content)

                        urlindex = content.find("/workspace/")
                        urlindex2 = content.find(".tbl")

                        table = requests.get(
                            "https://irsa.ipac.caltech.edu"
                            + content[urlindex:urlindex2]
                            + ".tbl"
                        )

                    except Exception as e:
                        self.logger.warn(
                            "plot_neoWISE: no connection to IRSA; url was:\n" + url
                        )
                        self.logger.warn(str(e))
                        data_log[tran_name] = False
                        info_list = []
                        wise_class = None
                        out_dict = None

                        continue

                    astro_tab = astropy.io.ascii.read(table.text)

                    # save the data
                    save_fname = f"{filebase}_neoWISE.txt"
                    # buf = io.BytesIO()
                    # buf.write(bytes(json.dumps(table.text, indent = 3), 'utf-8'))
                    # buf.seek(0)
                    self.put(save_fname, bytes(table.text, "utf-8"))
                    data_log[tran_name] = True

                dd = np.array(astro_tab)

                if len(dd) == 0:
                    info_list, wise_class, out_dict = ["neoWISE: no match"], None, None

                # apply some quality cuts
                iqual = (
                    (dd["saa_sep"] > 5)
                    * (dd["qual_frame"] > 0)
                    * (dd["cc_flags"] == "0000")
                )
                # if self.verbose:
                # print ('# neoWISE observations', len(dd))
                # print ('# data point to be removed by quality cuts', sum(iqual==False))

                # print ('apply_qcuts=', self.apply_qcuts)

                if sum(iqual) < 3:
                    if self.verbose & (sum(iqual) > 0):
                        self.logger.info(
                            "not enough good data."
                        )  # + np.array(astro_tab[('qual_frame', 'saa_sep', 'cc_flags','w1mpro', 'w1sigmpro')]))
                    elif self.verbose:
                        self.logger.info("no data for this source.")
                    if self.apply_qcuts:
                        info_list, wise_class, out_dict = (
                            [
                                "neoWISE: not enough detections (N={0})".format(
                                    sum(iqual)
                                )
                            ],
                            None,
                            None,
                        )
                        continue

                if self.apply_qcuts:
                    dd = dd[iqual]
                else:
                    self.logger.info("not applying quality cuts")
                dd = dd[np.argsort(dd["mjd"])]

                xx = dd["mjd"] - t0

                # attempt to make smart bins that are catch the entire light curve
                roundtime = np.unique(np.round(xx, decimals=-2))  # bins in 100 days
                neartime = [
                    xx[np.argmin(np.abs(xx - rt))] for rt in roundtime
                ]  # catch nearest data
                bins = []  # find the edges of the bins
                for nt in neartime:
                    ii = abs(xx - nt) < 30
                    if sum(ii):
                        nb = [min(xx[ii]) - 0.1, max(xx[ii]) + 0.1]  # some padding
                        # avoid duplicates
                        if not nb in bins:
                            bins.append(nb)

                bins = list(np.array(bins).flatten())

                # compute: W1
                ii1 = dd["w1sigmpro"] > 0
                xbin1, ybin1 = self.binthem(
                    xx[ii1],
                    dd["w1mpro"][ii1],
                    dd["w1sigmpro"][ii1],
                    use_wmean=True,
                    std=True,
                    sqrtN=True,
                    bins=bins,
                    silent=True,
                )
                inz1 = ybin1[0, :] > 0

                if sum(ii1):
                    w1_med = np.median(dd["w1mpro"][ii1])
                    w1_N = sum(ybin1[3, :])
                    # w1_chi2 = sum((dd['w1mpro'][ii1]-w1_med)**2 / dd['w1sigmpro'][ii1]**2) / sum(ii1)  # use all data with reported uncertainty
                    inz1_chi = (
                        ybin1[3, :] > 2
                    )  # need enough detection to compute scatter
                    w1_chi2 = sum(
                        (ybin1[0, inz1_chi] - w1_med) ** 2 / ybin1[1, inz1_chi] ** 2
                    ) / sum(
                        inz1_chi
                    )  # use uncertainty estimated from scatter in each bin

                    if out_dict is not None:
                        out_dict["w1_mag_mean"] = ybin1[0, inz1]
                        out_dict["w1_mag_med"] = w1_med
                        out_dict["w1_mag_rms"] = ybin1[1, inz1]
                        out_dict["w1_N_obs"] = ybin1[3, inz1]
                        out_dict["w1_time"] = xbin1[inz1]

                else:
                    w1_N, w1_med, w1_chi2 = 0, np.nan, np.nan

                # plot and compute: W2

                ii2_mag = dd["w2sigmpro"] > 0  # significant detections
                ii2 = dd["w2sigflux"] > 0  # all detections
                if sum(ii2) and sum(ii2_mag):

                    zp2 = np.median(
                        dd["w2mpro"][ii2_mag] + 2.5 * np.log10(dd["w2flux"][ii2_mag])
                    )

                    xbin2, ybin2 = self.binthem(
                        xx[ii2],
                        dd["w2flux"][ii2],
                        dd["w2sigflux"][ii2],
                        use_wmean=True,
                        std=True,
                        sqrtN=True,
                        bins=bins,
                        silent=True,
                    )

                    inz2 = (
                        ybin2[0, :] > 0
                    )  # note, this removed epoch where the mean flux is negative

                    # convert mean flux to mag
                    w2binflux = ybin2[0, :].copy()
                    ybin2[0, :] = -2.5 * np.log10(w2binflux) + zp2
                    ybin2[1, :] = 2.5 / np.log(10) * ybin2[1, :] / w2binflux

                    w2_med = (
                        -2.5
                        * np.log10(self.wmean(dd["w2flux"][ii2], dd["w2sigflux"][ii2]))
                        + zp2
                    )  # convert to mag
                    w2_N = sum(ybin2[3, :])

                    # w2_chi2 = sum((dd['w2mpro'][ii2]-w2_med)**2 / dd['w2sigmpro'][ii2]**2) / sum(ii2) # use all data with reported uncertainty
                    inz2_chi = (
                        ybin2[3, :] > 2
                    )  # need enough detection to compute scatter
                    w2_chi2 = sum(
                        (ybin2[0, inz2_chi] - w2_med) ** 2 / ybin2[1, inz2_chi] ** 2
                    ) / sum(
                        inz2_chi
                    )  # use uncertainty estimated from scatter in each bin

                    if out_dict is not None:
                        out_dict["w2_mag_mean"] = ybin2[0, inz2]
                        out_dict["w2_mag_med"] = w2_med
                        out_dict["w2_mag_rms"] = ybin2[1, inz2]
                        out_dict["w2_N_obs"] = ybin2[3, inz2]
                        out_dict["w2_time"] = xbin2[inz2]

                else:
                    w2_N, w2_med, w2_chi2 = 0, np.nan, np.nan

                # compute W1-W2 in each bin
                w1minw2 = []
                if (sum(ii1) > 10) and (sum(ii2) > 10) and (sum(ii2_mag) > 5):
                    for x in xbin2[inz2]:
                        # check that we have data in the same bins
                        i1_near = abs(xbin1[inz1] - x) < 20
                        i2_near = abs(xbin2[inz2] - x) < 20
                        n1_near = sum(ybin1[3, :][inz1][i1_near])
                        n2_near = sum(ybin2[3, :][inz2][i2_near])
                        if (n1_near > 2) and (n2_near > 2):
                            w1_near = np.mean(ybin1[0, :][inz1][i1_near])
                            w2_near = np.mean(ybin2[0, :][inz2][i2_near])
                            w1minw2.append(w1_near - w2_near)

                # do classification
                wise_class = None

                # and make a nice string that we can return
                if w2_N == 0 and w1_N == 0:
                    wise_info = "neoWISE: no match"
                else:
                    wise_info = "neoWISE: N(w1)={0:0.0f}".format(w1_N)

                # keep the string small
                # if w1_N>0:
                #   wise_info += '; <w1>={0:0.2f}'.format(w1_med)
                # if w2_N>0:
                #   wise_info += '; <w2>={0:0.2f}'.format(w2_med)

                if out_dict is not None:

                    if not np.isnan(np.median(w1minw2)):
                        wise_info += "; <w1-w2>={0:0.2f}".format(np.median(w1minw2))
                        out_dict["<w1-w2>"] = np.median(w1minw2)

                    if not np.isnan(w1_chi2):
                        wise_info += "; chi2(w1)={0:0.1f}".format(w1_chi2)
                        out_dict["chi2(w1)"] = w1_chi2

                    if not np.isnan(w2_chi2):
                        wise_info += "; chi2(w2)={0:0.1f}".format(w2_chi2)
                        out_dict["chi2(w2)"] = w2_chi2

                # we want at least 10 W2 detections over 3 epochs
                if w2_N > 10:
                    if sum(ybin2[3, inz2]) > 10 and (len(xbin2[inz2]) > 3):

                        # Stern12
                        # if (w2_med<15) and (np.median(w1minw2)>0.8):
                        #   wise_class = 'AGN'

                        # Assaf+13
                        if np.median(w1minw2) > 0.662 * np.exp(
                            0.232 * (np.clip(w2_med - 13.97, 0, 1e99)) ** 2
                        ):
                            wise_class = "AGN"
                            # wise_info+='; AGN classification based on color (Assef+12)'

                    # significant variability
                    if w1_chi2 > 10:
                        wise_class = "AGN?"
                        # wise_info+='; significant variability'
                    # evidence for variability
                    elif w1_chi2 > 5:
                        wise_class = None
                        # wise_info+='; some evidence for variability'

                info_list = [wise_info]

                # also add allWISE info
                if max(dd["w2mpro_allwise"]) and max(dd["w1mpro_allwise"]):
                    w1w2_allwise = max(dd["w1mpro_allwise"]) - max(dd["w2mpro_allwise"])
                    w1w2_allwise_err = np.sqrt(
                        max(dd["w1sigmpro_allwise"]) ** 2
                        + max(dd["w2sigmpro_allwise"]) ** 2
                    )
                    all_info = "allWISE: w1-w2={0:0.2f}+/-{1:0.2f}".format(
                        w1w2_allwise, w1w2_allwise_err
                    )

                    if np.median(w1w2_allwise) > 0.662 * np.exp(
                        0.232 * (np.clip(dd[0]["w2mpro_allwise"] - 13.97, 0, 1e99)) ** 2
                    ):
                        wise_class = "AGN"
                        all_info += "; AGN classification based on color (Assef+12)"

                    info_list += [all_info]

                save_fname = filebase + "_neoWISE.json"
                buf = io.BytesIO()
                buf.write(
                    bytes(
                        json.dumps(
                            self.jsonify(
                                {
                                    "info_list": info_list,
                                    "wise_class": wise_class,
                                    "out_dict": out_dict,
                                }
                            ),
                            indent=3,
                        ),
                        "utf-8",
                    )
                )
                buf.seek(0)
                self.put(save_fname, buf.read())

                if self.verbose:
                    # print ('allWISE: median(w1)-median(w2) = {0:0.3f}'.format(w1_med-w2_med))
                    # [print(x) for x in info_list]
                    self.logger.info("neoWISE classification: " + str(wise_class))

                # make a plot
                plt.close()
                fig = plt.figure()
                ax = fig.add_subplot(111)

                if w1_N:
                    line = ax.errorbar(
                        xx[ii1],
                        dd["w1mpro"][ii1],
                        dd["w1sigmpro"][ii1],
                        fmt="s",
                        alpha=0.5,
                        label="W1",
                    )
                    ax.plot(
                        xbin1[inz1], ybin1[0, inz1], "d--", color=line[0].get_color()
                    )
                    ax.fill_between(
                        xbin1[inz1],
                        ybin1[0, inz1] + ybin1[1, inz1],
                        ybin1[0, inz1] - ybin1[1, inz1],
                        color=line[0].get_color(),
                        alpha=0.5,
                    )
                    if self.plot_allWISE and max(dd["w1mpro_allwise"]):
                        ax.plot(
                            xx,
                            np.repeat(max(dd["w1mpro_allwise"]), len(xx)),
                            ":",
                            lw=2,
                            color=line[0].get_color(),
                            alpha=0.6,
                        )

                if w2_N:
                    line = ax.errorbar(
                        xx[ii2],
                        dd["w2mpro"][ii2],
                        dd["w2sigmpro"][ii2],
                        fmt="o",
                        alpha=0.5,
                        label="W2",
                    )
                    ax.plot(
                        xbin2[inz2], ybin2[0, inz2], "d--", color=line[0].get_color()
                    )
                    ax.fill_between(
                        xbin2[inz2],
                        ybin2[0, inz2] + ybin2[1, inz2],
                        ybin2[0, inz2] - ybin2[1, inz2],
                        color=line[0].get_color(),
                        alpha=0.5,
                    )
                    if self.plot_allWISE and max(dd["w2mpro_allwise"]):
                        ax.plot(
                            xx,
                            np.repeat(max(dd["w2mpro_allwise"]), len(xx)),
                            ":",
                            lw=2,
                            color=line[0].get_color(),
                            alpha=0.6,
                        )

                ax.set_ylim(ax.get_ylim()[::-1])
                ax.set_xlabel("MJD - {0:0.0f}".format(t0))
                ax.set_ylabel("mag (Vega)")
                titstr = "N(w1)={0:0.0f}; N(w2)={1:0.0f}; chi2(w1)={2:0.1f}; chi2(w2)={3:0.1f} <w1-w2>={4:0.2f}".format(
                    w1_N, w2_N, w1_chi2, w2_chi2, np.median(w1minw2)
                )
                ax.set_title(titstr, fontsize=11)
                ax.legend()
                fig.tight_layout()

                plot_fname = f"{filebase}-neoWISE.pdf"
                buf = io.BytesIO()
                ax.figure.savefig(buf, format="pdf")
                buf.seek(0)
                self.put(plot_fname, buf.read())
                self.logger.debug(f"plotting {tran_name}")
                ax.clear()

        buf = io.BytesIO()
        buf.write(bytes(json.dumps(data_log, indent=3), "utf-8"))
        buf.seek(0)
        self.put(self.save_location + "/neowise_log.json", buf.read())

        self.commit()
        return None

    def binthem(
        self,
        x,
        y,
        yerr=None,
        bins=10,
        range=[],
        use_mean=False,
        use_sum=False,
        use_wmean=False,
        cl=0.9,
        poisson=False,
        std=False,
        sqrtN=True,
        smart_expand=False,
        silent=False,
    ):
        """
        >> xmid, ymid = binthem(x, y)
        bin parameter y using bins in x-direction
        output:
         - xmid  (N) the mean value of x for each bin
         - ymid  (4,N) containing the [median/mean/total value, -/+ of uncertainty/dispersion, number in bin]
        input:
         - x,y  equal length arrays
         - bins  number of xbins or array with bins (default is 10)
        optional input:
         - range=[xmin, xmax] range for the bins (default is [min(x),max(x)])
         - cl=0.9  confidence level for computing the uncertainty
         - poisson=False use Poisson stat to find uncertainty on number in bin
         - std=False  use standard deviation / sqrt(N) to compute uncertainty (ie, symtric errorbars)
         - sqrtN=True  set to False to use std only
         - use_mean=False  use mean (median is default)
         - use_sum=False sum bin content
         - use_wmean=False  compute weighted mean, requires yerr= input
         - smart_expand=False, if set, we try to increase bins to catch clusters of data (warning: experimentenal!)
         - silent  shutup
        """

        if use_wmean:
            if yerr is None:
                self.logger.warn("ERROR, please give yerr= input when wmean=True")
                return

        if np.isscalar(bins):
            if len(range) == 2:
                x_bins = np.linspace(range[0], range[1], bins)
            else:
                x_bins = np.linspace(np.min(x), np.max(x), bins)
        else:
            x_bins = bins.copy()

        if smart_expand:

            for i in np.arange(1, len(x_bins)):
                dbin = x_bins[i] - x_bins[i - 1]
                iright = x > x_bins[i]
                if sum(iright):
                    # points closer to the current bins than next bins
                    if min(x[iright] - x_bins[i]) < dbin / 2.0:
                        x_bins[i] = x_bins[i] + dbin / 2.0
        xmid = np.zeros(len(x_bins) - 1)
        ymid = np.zeros((5, len(x_bins) - 1))

        for i in np.arange(len(xmid)):

            ibin = np.where((x >= x_bins[i]) & (x < x_bins[i + 1]))[0]
            if i == (len(xmid) - 1):  # close last bin
                ibin = np.where((x >= x_bins[i]) & (x <= x_bins[i + 1]))[0]

            if len(ibin) > 0:

                xmid[i] = np.mean(x[ibin])
                ymid[4, i] = np.std(x[ibin])

                y_arr = y[ibin]
                ymid[3, i] = len(ibin)

                if use_mean:
                    ymid[0, i] = np.mean(y_arr)
                elif use_sum:
                    ymid[0, i] = np.sum(y_arr)
                elif use_wmean:
                    y_arr_err = yerr[ibin]
                    ymid[0, i] = self.wmean(y_arr, 1 / y_arr_err * 2)
                else:
                    ymid[0, i] = np.median(y_arr)

                if std:
                    ymid[[1, 2], i] = np.std(y_arr)
                    if sqrtN:
                        ymid[[1, 2], i] = ymid[[1, 2], i] / np.sqrt(len(ibin))
                elif poisson and (len(ibin) > 1):
                    ymid[[1, 2], i] = np.abs(
                        poisson_limits(np.sum(y[ibin]), cl) / len(ibin) - ymid[0, i]
                    )
                elif use_wmean:
                    ymid[[1, 2], i] = 1 / np.sqrt(sum(1 / y_arr_err**2))
                else:
                    y_arrs = np.sort(y_arr)
                    ii = np.arange(len(ibin))
                    ymid[1, i] = np.abs(
                        ymid[0, i] - np.interp((cl / 2.0 - 0.5) * len(ibin), ii, y_arrs)
                    )
                    ymid[2, i] = np.abs(
                        ymid[0, i] - np.interp((cl / 2.0 + 0.5) * len(ibin), ii, y_arrs)
                    )
            else:
                xmid[i] = (x_bins[i] + x_bins[i + 1]) / 2.0

        if sum(ymid[3, :]) > len(x):
            self.logger.warn(
                "binthem: WARNING: more points in bins ({0}) compared to lenght of input ({1}), please check your bins".format(
                    sum(ymid[3, :]), len(x)
                )
            )
            # key = input()

        return xmid, ymid

    def wmean(self, d, ivar):
        """
        >> mean = wmean(d, ivar)
        inverse-variance weighted mean
        """
        d = np.array(d)
        ivar = np.array(ivar)
        return np.sum(d * ivar) / np.sum(ivar)
