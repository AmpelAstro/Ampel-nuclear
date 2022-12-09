#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t2/T2FlexFit.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu
# Date              : 28.04.2020
# Last Modified Date: 08.06.2020
# Last Modified By  : mitchell@nyu.edu

import sys
import numpy as np
import logging

logging.basicConfig()

from ampel.nuclear.flexfit import flex_metrics
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.view import LightCurve
from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.ztf.util.ZTFIdMapper import to_ztf_id


class T2FlexFit(AbsLightCurveT2Unit):
    """
    Least squares fitting on (ZTF) data with observations in two bands
    We fit for the peak, using a Gaussian rise and an exponential decay
    We measure the mean color pre-peak, and allow for a color-change post-peak

    Goal is photometric typing of alerts into:

     - "SN" or "Ia" : based on rise time
     - "AGN"        : slow rise or fading
     - "CV"         : fast
     - "TDE"        : not the above and constant color or blue
    """

    oldest_upper_limits: int
    max_post_peak: int

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
        Dictionary of three dicts: fit parameters, lightcurve data, and plotting info:
        {'fit_params':
        ['dtime_peak', 'band', 'mag_peak', 'sigma_rise', 'sigma_fade', 'mean_color', 'color_slope', 'e_dtime_peak',
        'e_mag_peak', 'e_mean_color', 'e_color_slope', 'e_sigma_rise', 'mad', 'rms', 'chi2', 'n_rise', 'n_fade',
        'classification', 'name', 'mjd_peak', 'photoclass', 'blabber'],

        'lc_data':
        ['dtime', 'flux', 'flux_err', 'fid', 'isdetect', 'jd', 'jd0'],

        'plot_info':
        ['lsq_out', 'dt_fix', 'dtime_out', 'flux_out', , 'pid', 'programid', 'finfo']
        }

        """
        lc = get_raw_lc(light_curve)
        name = to_ztf_id(light_curve.get_photopoints()[0]["stock"][0])
        meta_data = {
            "name": name,
            "oldest_upper_limits": self.oldest_upper_limits,
            "max_post_peak": self.max_post_peak,
            "logger": self.logger,
        }
        result = flex_metrics.compute(meta_data, lc)
        return jsonify(result)


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


def _append_limits(lc, lcurve):
    # lc is dict, lcurve is LightCurve instance
    # If there's a more efficient way to do this with Ampel let us know!
    if lcurve.get_values("jd", of_upper_limits=True):
        ul_jd = lcurve.get_values("jd", of_upper_limits=True)
        ul_fid = lcurve.get_values("fid", of_upper_limits=True)
        ul_diffmaglim = lcurve.get_values("diffmaglim", of_upper_limits=True)
        ul_isdiffpos = np.full(len(ul_jd), "t")

        lc["jd"] = np.append(lc["jd"], ul_jd)
        lc["fid"] = np.append(lc["fid"], ul_fid)
        lc["diffmaglim"] = np.append(lc["diffmaglim"], ul_diffmaglim)
        lc["isdiffpos"] = np.append(lc["isdiffpos"], ul_isdiffpos)

        for key in lc:
            if key not in ["jd", "fid", "diffmaglim", "isdiffpos"]:
                lc[key] = np.append(lc[key], np.full(len(ul_jd), 999))

    return lc


def get_raw_lc(light_curve):
    lc = _lc_to_dict(light_curve)
    lc = _append_limits(lc, light_curve)
    return lc
