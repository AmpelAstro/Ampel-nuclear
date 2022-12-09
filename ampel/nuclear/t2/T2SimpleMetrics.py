#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t2/T2SimpleMetrics.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu
# Date              : 11.05.2020
# Last Modified Date: 08.06.2020
# Last Modified By  : mitchell@nyu.edu

import numpy as np
import os, sys
import numpy as np
from astropy import units as u  # type: ignore
from astropy.coordinates import SkyCoord  # type: ignore
import astropy.time  # type: ignore
import pandas as pd  # type: ignore

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.ztf.util.ZTFIdMapper import to_ztf_id
from ampel.view import LightCurve


class T2SimpleMetrics(AbsLightCurveT2Unit):
    """
    Computes simple lightcurve metrics
    """

    def process(self, light_curve) -> UBson | UnitResult:
        """
        Parameters
        ----------
        light_curve: `ampel.view.LightCurve` instance.
        See original docstring & photopoint docstring.
        Note: arguments of '.get_values()' method are strings and are in the following table: https://zwickytransientfacility.github.io/ztf-avro-alert/schema.html

        run_config: 'dict'
        'local_data': plotting folder
        'name': ZTF name of object,

        Returns
        -------
        Dict of two dicts:
        {
        metrics: Dict of simple metrics and their values
        plot_info: Info used for plotting t3

        }
        """
        lc = _lc_to_dict(light_curve)
        name = to_ztf_id(light_curve.get_photopoints()[0]["stock"][0])
        meta_data = {"name": name}
        result = jsonify(compute(meta_data, lc, self.logger))
        return result


def _lc_to_dict(lc):
    # Turns Ampel Lightcurve object into a dict
    keys = (
        "jd",
        "fid",
        "pid",
        "diffmaglim",
        "pdiffimfilename",
        "programpi",
        "programid",
        "candid",
        "isdiffpos",
        "tblid",
        "nid",
        "rcid",
        "field",
        "xpos",
        "ypos",
        "ra",
        "dec",
        "magpsf",
        "sigmapsf",
        "chipsf",
        "magap",
        "sigmagap",
        "distnr",
        "magnr",
        "sigmagnr",
        "chinr",
        "sharpnr",
        "sky",
        "magdiff",
        "fwhm",
        "classtar",
        "mindtoedge",
        "magfromlim",
        "seeratio",
        "aimage",
        "bimage",
        "aimagerat",
        "bimagerat",
        "elong",
        "nneg",
        "nbad",
        "rb",
        "rbversion",
        "drb",
        "drbversion",
        "ssdistnr",
        "ssmagnr",
        "ssnamenr",
        "sumrat",
        "magapbig",
        "sigmagapbig",
        "ranr",
        "decnr",
        "ndehist",
        "ncovhist",
        "jdstarhist",
        "jdendhist",
        "scorr",
        "tooflag",
        "objectidps1",
        "sgmag1",
        "srmag1",
        "simag1",
        "szmag1",
        "sgscore1",
        "distpsnr1",
        "objectidps2",
        "sgmag2",
        "srmag2",
        "simag2",
        "szmag2",
        "sgscore2",
        "distpsnr2",
        "objectidps3",
        "sgmag3",
        "srmag3",
        "simag3",
        "szmag3",
        "sgscore3",
        "distpsnr3",
        "nmtchps",
        "rfid",
        "jdstartref",
        "jdendref",
        "nframesref",
        "dsnrms",
        "ssnrms",
        "dsdiff",
        "magzpsci",
        "magzpsciunc",
        "magzpscirms",
        "nmatches",
        "clrcoeff",
        "clrcounc",
        "zpclrcov",
        "zpmed",
        "clrmed",
        "clrrms",
        "neargaia",
        "neargaiabright",
        "maggaia",
        "maggaiabright",
        "exptime",
    )

    vals = []

    for key in keys:
        vals.append(np.array(lc.get_values(key)))

    return dict(zip(keys, vals))


def wmean(d, ivar):
    """
    >> mean = wmean(d, ivar)
    inverse-variance weighted mean
    """
    d = np.array(d)
    ivar = np.array(ivar)
    return np.sum(d * ivar) / np.sum(ivar)


def jsonify(obj):
    """
    Recursively cast to JSON native types
    """
    if hasattr(obj, "tolist"):
        return obj.tolist()
    elif hasattr(obj, "keys"):
        return {jsonify(k): jsonify(obj[k]) for k in obj.keys()}
    elif hasattr(obj, "__len__") and not isinstance(obj, str):
        return [jsonify(v) for v in obj]
    else:
        return obj


def compute(meta_data, lc, logger, verbose=False):

    # make the rec
    result = {}
    plot_info = {}

    ztfname = meta_data["name"]

    # can be useful to remake offset plots, including ps1_data
    plot_info["plot_fname"] = ztfname + "_offset.pdf"
    plot_info["offset_plot"] = False

    # force remake plot if source is still young
    result["name"] = ztfname

    result["ra"], result["dec"] = lc["ra"], lc["dec"]

    data_dict1 = {
        key: value for key, value in lc.items() if len(value) == len(lc["jd"])
    }
    photoarr0 = pd.DataFrame(data_dict1)

    isdetect = photoarr0["magpsf"] < 99  # should be everything
    idiffspos = photoarr0["isdiffpos"] == "t"

    result["ndetections"] = int(sum(idiffspos & isdetect))
    result["ndetections_negative"] = int(sum((idiffspos == False) & isdetect))

    source_age = astropy.time.Time.now().jd - max(lc["jd"])
    # logging.info('simple_metrics: {0} ({1});  age={2:4.1f}d; ra,dec = {3:0.5f} {4:0.5f}; ndetection={5}'.format(ztfname, str(meta_data['classification']), source_age, result['ra'],result['dec'], result['ndetections'] ))

    # we only work with positive detections
    if result["ndetections"] < 3:
        logger.warn("WARNING. not enough positive detection for this source!")

        return {"metrics": result, "plot_info": plot_info}

    photoarr = photoarr0[idiffspos & isdetect]
    photoarr = photoarr.iloc[np.argsort(photoarr["jd"])]

    # here we fill in the metric
    ipeak = photoarr["magpsf"] == min(
        photoarr["magpsf"]
    )  # not very robust but OK for now
    ilatest = photoarr["jd"] == max(photoarr["jd"])

    result["peak_mag"] = min(photoarr["magpsf"])
    result["latest_mag"] = max(photoarr[ilatest]["magpsf"])
    result["peak_jd"] = min(photoarr[ipeak]["jd"])
    result["latest_jd"] = max(photoarr[ilatest]["jd"])
    result["age"] = astropy.time.Time.now().jd - result["latest_jd"]

    for flt in [1, 2]:

        photoarr_flt = photoarr[photoarr["fid"] == flt]
        flt_name = "g" if flt == 1 else "r" if flt == 2 else None
        if "ranr" in photoarr.columns:
            result["ndetections_{0}".format(flt_name)] = len(photoarr_flt)
        else:
            result["ndetections_{0}".format(flt_name)] = 0

        if result["ndetections_{0}".format(flt_name)] > 2:

            dra = (
                photoarr_flt["ra"] - photoarr_flt["ranr"]
            ) * 3600  # offset to nearest ZTF PSF catalog source
            ddec = (photoarr_flt["dec"] - photoarr_flt["decnr"]) * 3600
            drms = np.clip(np.std(np.append(dra, ddec)), 0.1, 1)
            dstd = np.clip(drms / np.sqrt(len(photoarr_flt)), 0.05, 1)

            # offset_sigma = np.clip(0.27 + 1.66 * (photoarr_flt['sigmamagpsf']-0.1) , 0.1, 10)  # equation from the NedStark paper
            sigma_offset = np.clip(
                0.24 + 0.04 * (photoarr_flt["magpsf"] - 20), 0.1, 1
            )  # equation from the NedStark paper
            snr_weight = 1 / sigma_offset**2

            if "magnr" in photoarr.columns:
                sigma_offset_nr = np.clip(
                    0.24 + 0.04 * (np.median(photoarr_flt["magnr"]) - 20), 0.1, 1
                )  # / np.sqrt(5) # assume x images make a ref frame
            else:
                sigma_offset_nr = 0.05

            dunc_theo = np.sqrt(
                1 / (sum(1 / sigma_offset**2)) + sigma_offset_nr**2
            )  # can be too optimistic if large number of observations
            dunc_rms = np.sqrt(
                dstd**2 + sigma_offset_nr**2
            )  # can be too optimisiric if small number of observations
            dunc = np.max([dunc_rms, dunc_theo])

            weighted_dra = wmean(dra, snr_weight)
            weighted_dec = wmean(ddec, snr_weight)
            logger.info(
                "weighted (dra,drec) {0:0.5f} {1:0.5f}".format(
                    weighted_dra, weighted_dec
                )
            )
            offset_weighted = np.sqrt(weighted_dra**2 + weighted_dec**2)

            result["offset_rms_{0}".format(flt_name)] = np.std(np.append(dra, ddec))
            result["offset_unc_{0}".format(flt_name)] = dunc
            result["offset_med_{0}".format(flt_name)] = np.sqrt(
                np.median(dra) ** 2 + np.median(ddec) ** 2
            )
            result["offset_weighted_{0}".format(flt_name)] = offset_weighted
            result["offset_sig_{0}".format(flt_name)] = offset_weighted / dunc

            logger.info(
                str(flt_name)
                + ": "
                + "simple_metrics: distance to ref (median, weighted, significance): {0:0.3f} {1:0.3f} {2:0.3f}".format(
                    result["offset_med_{0}".format(flt_name)],
                    result["offset_weighted_{0}".format(flt_name)],
                    result["offset_sig_{0}".format(flt_name)],
                )
            )
            logger.info(
                str(flt_name)
                + ": "
                + "simple_metrics: distance to ref (rms, unc_rm, unc_theo) : {0:0.3f} {1:0.3f} {2:0.3f}".format(
                    result["offset_rms_{0}".format(flt_name)], dunc_rms, dunc_theo
                )
            )

        else:
            result["offset_sig_{0}".format(flt_name)] = -999

    # select which offset measurement to use based on rm
    if (result["ndetections_g"] > 2) and (result["ndetections_r"] > 2):
        if result["offset_unc_g"] < result["offset_unc_r"]:
            result["offset_sig"] = result["offset_sig_g"]
        else:
            result["offset_sig"] = result["offset_sig_r"]
    elif result["ndetections_r"] > 2:
        result["offset_sig"] = result["offset_sig_r"]
    elif result["ndetections_g"] > 2:
        result["offset_sig"] = result["offset_sig_g"]

    result["lastest_obs_is_detection"] = int(
        max(photoarr0["jd"]) == result["latest_jd"]
    )
    result["latest_maglim"] = min(
        photoarr[ilatest]["diffmaglim"]
    )  # assumed limmag meant diffmaglim

    # by default assume we have no detection in ref image
    result["ref_mag_r"], result["ref_mag_g"] = 23, 23

    # get ref magnitude if detected (note cut on distnr is needed to avoid matches to nearby unrelated sources for non-detections in a single band)
    if "magnr" in photoarr.columns:

        for flt in (1, 2):

            flt_name = "g" if flt == 1 else "r" if flt == 2 else None
            ii = (photoarr["fid"] == flt) & (photoarr["distnr"] < 2)
            if sum(ii):
                result["ref_mag_" + flt_name] = np.median(photoarr[ii]["magnr"])
    else:
        logger.warn(f"warning: magnr not found in photoarr: {photoarr.columns}")

    # now compute magnitude difference
    for flt in (1, 2):

        flt_name = "g" if flt == 1 else "r" if flt == 2 else None
        ii = photoarr["fid"] == flt

        if sum(ii):
            ip = photoarr["magpsf"][ii] == min(photoarr[ii]["magpsf"])
            il = photoarr[ii]["jd"] == max(photoarr[ii]["jd"])

            result["peak_mag_" + flt_name] = (photoarr["magpsf"][ii][ip]).values[0]
            result["latest_mag_" + flt_name] = (photoarr["magpsf"][ii][il]).values[0]
            result["peak_diff_mag_" + flt_name] = (
                (photoarr["magpsf"][ii] - result["ref_mag_" + flt_name])[ip]
            ).values[0]
            result["latest_diff_mag_" + flt_name] = (
                (photoarr["magpsf"][ii] - result["ref_mag_" + flt_name])[il]
            ).values[0]
            logger.info(
                str(flt_name)
                + ": "
                + "(ref_mg, peak_diff_mag)  = ({0:0.2f}, {1:0.2f})".format(
                    result["ref_mag_" + flt_name], result["peak_diff_mag_" + flt_name]
                )
            )
        else:
            result["peak_mag_" + flt_name] = 999
            result["latest_mag_" + flt_name] = 999
            result["peak_diff_mag_" + flt_name] = 999
            result["latest_diff_mag_" + flt_name] = 999

    # finally, store the largest flux increase between the two bands
    result["ref_mag"] = min(result["ref_mag_r"], result["ref_mag_g"])
    result["latest_diff_mag"] = min(
        result["latest_diff_mag_r"], result["latest_diff_mag_g"]
    )
    result["peak_diff_mag"] = min(result["peak_diff_mag_r"], result["peak_diff_mag_g"])

    if "ranr" in photoarr.columns:
        plot_info["offset_plot"] = True
        plot_info["dra"] = (
            photoarr["ra"] - photoarr["ranr"]
        ) * 3600  # offset to nearest ZTF PSF catalog source
        plot_info["ddec"] = (photoarr["dec"] - photoarr["decnr"]) * 3600

        sigma_offset = np.clip(
            0.24 + 0.04 * (photoarr["magpsf"] - 20), 0.1, 1
        )  # equation from the NedStark paper
        plot_info["snr_weight"] = 1 / sigma_offset**2

        plot_info["ig"] = photoarr["fid"] == 1

        plot_info["ir"] = photoarr["fid"] == 2

    return {"metrics": result, "plot_info": plot_info}
