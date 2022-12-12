"""
function to compute simple light curve metrics
"""
import os
import numpy as np
import pandas as pd  # type: ignore
import json
import astropy
import warnings

from . import flexfit


name = "flex"

dt = [
    ("name", "U35"),
    ("classification", "U35"),
    ("photoclass", "U35"),
    ("blabber", "U200"),
    ("band", "U6"),
    ("rms", "f8"),
    ("mad", "f8"),
    ("chi2", "f8"),
    ("mjd_peak", "f8"),
    ("dtime_peak", "f8"),
    ("e_dtime_peak", "f8"),
    ("mag_peak", "f8"),
    ("e_mag_peak", "f8"),
    ("sigma_rise", "f8"),
    ("sigma_fade", "f8"),
    ("e_sigma_rise", "f8"),
    ("e_sigma_fade", "f8"),
    ("n_rise", "f8"),
    ("n_fade", "f8"),
    ("mean_color", "f8"),
    ("color_slope", "f8"),
    ("e_mean_color", "f8"),
    ("e_color_slope", "f8"),
]


# magnitude/flux conversion with some zp
zp = 26


def flux2mag(flux, zp=zp):
    return -2.5 * np.log10(flux) + zp


def mag2flux(mag, zp=26):
    return 10 ** (-0.4 * (mag - zp))


def jdtomjd(jd):
    return jd - 2400000.5


exception_dict = {}
exception_dict["ZTF17aaazdba"] = {
    "jd_min": 2458515
}  # has early bogus data, so we require jd larger than jd_min


def compute(meta_data, lc):
    """
    Returns {'fit_params':
    ['dtime_peak', 'band', 'mag_peak', 'sigma_rise', 'sigma_fade', 'mean_color', 'color_slope', 'e_dtime_peak',
    'e_mag_peak', 'e_mean_color', 'e_color_slope', 'e_sigma_rise', 'mad', 'rms', 'chi2', 'n_rise', 'n_fade',
    'classification', 'name', 'mjd_peak', 'photoclass', 'blabber'],

    'lc_data':
    ['dtime', 'flux', 'flux_err', 'fid', 'isdetect', 'jd', 'jd0', 'pid', 'programid'],

    'plot_info':
    ['lsq_out', 'dt_fix', 'dtime_out', 'flux_out', , 'pid', 'programid', 'finfo']
    }
    """

    # allow some manual trimming
    if meta_data["name"] in exception_dict:
        edict = exception_dict[meta_data["name"]]
    else:
        edict = None

    lc_data, fit_params, plot_info = {}, {}, {}
    logger = meta_data["logger"]
    dtime, flux, flux_err, fid, isdetect, jd0, magnr_dict, mdata = trim_ztf_ampel_data1(
        lc,
        logger,
        exception_dict=edict,
        oldest_upper_limits=meta_data["oldest_upper_limits"],
        max_post_peak=meta_data["max_post_peak"],
    )
    (
        lc_data["dtime"],
        lc_data["flux"],
        lc_data["flux_err"],
        lc_data["fid"],
        lc_data["isdetect"],
    ) = (dtime, flux, flux_err, fid, isdetect)

    if sum(isdetect) < 3:
        lc_data["jd"] = mdata["jd"]
        lc_data["jd0"] = jd0
        lc_data["pid"] = mdata["pid"]
        lc_data["programid"] = mdata["programid"]

        fit_params = {}
        fit_params["name"] = meta_data["name"]
        fit_params["photoclass"] = "not enough data"
        return {
            "fit_params": fit_params,
            "lightcurve_data": lc_data,
            "plot_info": None
            # {'lsq_out' : None,
            # 'dt_fix': None,
            # 'dtime_out': None,
            # 'flux_out': None,
            # 'pid': None,
            # 'programid': None,
            # 'finfo': None}
        }

    # run flexfit, get dictionary with fit info
    fit_params, lc_data, plot_info = flexfit.flex_fit_wclip(
        dtime, flux, flux_err, fid, logger, isdetect, verbose=False
    )
    lc_data["jd"] = mdata["jd"]
    lc_data["jd0"] = jd0
    lc_data["pid"] = mdata["pid"]
    lc_data["programid"] = mdata["programid"]

    # add some meta data
    fit_params["name"] = meta_data["name"]

    # collect some info to show on plots
    if not "e_dtime_peak" in fit_params:
        str_em = "fix"
    elif fit_params["e_dtime_peak"] < 0:
        str_em = "fix"
    else:
        str_em = "{0:0.1f}".format(fit_params["e_dtime_peak"])

    if not "dtime_peak" in fit_params:
        fit_params["dtime_peak"] = -999
    if "mag_peak" in fit_params:
        finfo = "t_peak={0:0.1f} ({1}); m_peak={2:0.1f}; rise={3:0.1e}; fade={4:0.1e}; chi2={5:0.1f}".format(
            fit_params["dtime_peak"],
            str_em,
            fit_params["mag_peak"],
            fit_params["sigma_rise"],
            fit_params["sigma_fade"],
            fit_params["chi2"],
        )
    else:
        finfo = "no fit"

    if "mean_color" in fit_params:
        if not "e_color_slope" in fit_params:
            fit_params["e_color_slope"] = -999
        if not "e_mean_color" in fit_params:
            fit_params["e_mean_color"] = -999
        finfo += "\ncolor={0:0.2f} +/- {1:0.2f}; slope={2:0.3f} +/- {3:0.3f}".format(
            fit_params["mean_color"],
            fit_params["e_mean_color"],
            fit_params["color_slope"],
            fit_params["e_color_slope"],
        )

    logger.info(f"flex_metrics: {finfo}")
    plot_info["finfo"] = finfo

    # ---
    # get the class, plus some extra chatter
    blab, photoclass = flexfit.flex_class(fit_params)

    logger.info(f"flex_metrics: {photoclass}")
    logger.info(f"flex_metrics: {blab}")

    if jd0:
        fit_params["mjd_peak"] = jdtomjd(fit_params["dtime_peak"] + jd0)
    fit_params["photoclass"] = photoclass
    fit_params["blabber"] = blab

    return {
        "fit_params": fit_params,
        "lightcurve_data": lc_data,
        "plot_info": plot_info,
    }


def trim_ztf_ampel_data1(
    data_json,
    logger,
    oldest_upper_limits=14,
    filters=[1, 2],
    only_ztf=True,
    instruments=["ZTF", "SEDM", "IOO"],
    max_post_peak=200,
    remove_neg=True,
    exception_dict=exception_dict,
):
    """
    convenience function for dealing with data from the ampel

    >>> dtime, flux, flux_err, fid, isdetect, jd0, magnr_dict0, data = trim_ztf_ampel_data(input_json)

    input: dict created by reading output of uploaded_photometry

    upper limits (5-sigma) are stored in flux_err and have flux=0.
    for these limits I try to account for the flux of host galaxy (magnr), I we assume the event is nuclear

    options:
    - filter=[g,r], but can be other filters ;
    - only_ZTF=True, set False to include data from intruments;
    - instruments=['ZTF','SEDM', 'IOO']
    """

    dt = [
        ("name", "U35"),
        ("classification", "U35"),
        ("photoclass", "U35"),
        ("blabber", "U200"),
        ("band", "U6"),
        ("rms", "f8"),
        ("mad", "f8"),
        ("chi2", "f8"),
        ("mjd_peak", "f8"),
        ("dtime_peak", "f8"),
        ("e_dtime_peak", "f8"),
        ("mag_peak", "f8"),
        ("e_mag_peak", "f8"),
        ("sigma_rise", "f8"),
        ("sigma_fade", "f8"),
        ("e_sigma_rise", "f8"),
        ("e_sigma_fade", "f8"),
        ("n_rise", "f8"),
        ("n_fade", "f8"),
        ("mean_color", "f8"),
        ("color_slope", "f8"),
        ("e_mean_color", "f8"),
        ("e_color_slope", "f8"),
    ]

    data_dict1 = {
        key: value
        for key, value in data_json.items()
        if len(value) == len(data_json["magnr"])
    }  # everything that isn't a single value

    # I expect to see RuntimeWarnings in this block, median([])=nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        magnr_dict0 = {
            str(k): np.median(
                data_json["magnr"][
                    np.where(
                        (data_json["fid"] == k) & (data_json["magnr"] is not None)
                    )[0]
                ]
            )
            for k in (1, 2, 3)
        }

    data0 = pd.DataFrame(data_dict1)
    isdetect0 = data0["magpsf"] < 23
    islim0 = isdetect0 == False

    # some trimming
    if remove_neg:
        itrim0 = (data0["isdiffpos"].values == "t") | (
            data0["isdiffpos"].values == "1"
        )  # only positive detections
        # print('neg:', itrim0)
    else:
        itrim0 = np.repeat(True, len(data0))

    # don't use upper limits
    # (we can insert them later, after figure out how to account for host contribution adding to the uncertainty)
    itrim0 &= isdetect0

    # selection of bands
    itrim0 &= np.array([x in filters for x in data0["fid"]])

    # option to look only close to peak
    itrim0 &= (
        data0["jd"] - data0["jd"][data0["magpsf"].idxmin()]
    ) < max_post_peak  # upto xx days post guessed peak

    if exception_dict is not None:
        if "jd_min" in exception_dict:
            itrim0 &= data0["jd"] > exception_dict["jd_min"]

    # apply trimming
    data = data0[itrim0]
    # sort
    data = data.sort_values("jd")

    # get the limits and detections
    fid = np.array(
        [
            "g" if d == 1 else "r" if d == 2 else "i" if d == 3 else "check filter?"
            for d in data["fid"]
        ]
    )
    isdetect = data["magpsf"].values < 22
    islim = isdetect == False

    logger.info(
        "# data point used for fitting {0} (with {1} upper limits)".format(
            len(data), sum(islim)
        )
    )

    flux = mag2flux(data["magpsf"].values)
    flux_lim = mag2flux(
        data["diffmaglim"].values
    )  # 5 sigma limit for entire image (ie, not including host noise)
    flux[islim] = 0
    flux_err = np.log(10) / 2.5 * np.clip(data["sigmapsf"].values, 0.1, 2) * flux

    if remove_neg == False:
        ineg = data["isdiffpos"].values == 0
        flux[ineg] *= -1

    # try to guess the real flux upper limit by including the host contribution
    magnr_arr = np.array([magnr_dict0.get(k) for k in fid])
    magnr_arr[magnr_arr == None] = 25  # we dont have the host magnitude, dont use it...
    gain_corr_est = 5.0  # tuned for ZTF

    # flux_err = gain_corr_est * ( (flux_err)**2 + mag2flux(magnr_arr))**0.5 # not clear if this correction is needed, but I have seen residuals for ZTF difference photometry ontop of bright galaxies
    flux_err[islim] = (
        (flux_lim[islim] / 5) ** 2 + gain_corr_est**2 * mag2flux(magnr_arr[islim])
    ) ** 0.5

    if sum(isdetect):
        jd0 = min(data[isdetect]["jd"])
        dtime = data["jd"].values - jd0  # pick time normalization
    else:
        jd0 = None
        dtime = None

    return dtime, flux, flux_err, fid, isdetect, jd0, magnr_dict0, data.to_records()
