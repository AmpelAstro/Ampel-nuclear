"""
least squares fitting on (ZTF) data with observations in two bands
we fit for the peak, using a Gaussian rise and an exponential decay
we measure the mean color pre-peak, and allow for a color-change post-peak

goal is photometric typing of alerts into:

 - "SN" or "Ia" : based on rise time 
 - "AGN"		: slow rise or fading
 - "CV"			: fast
 - "TDE"		: not the above and constant color or blue

"""

import numpy as np

log10, log, exp, sqrt = np.log10, np.log, np.exp, np.sqrt
from scipy.optimize import leastsq
import astropy.stats


# flux/mag conversion without zeropoint (only used to force fit in mag space)
def f2m(flux):
    return -2.5 * log10(flux)


def m2f(mag):
    return 10 ** (-0.4 * (mag))


def line(p, x):
    return p[0] + p[1] * x


def gauss(x, sigma):
    return np.exp(-0.5 * x**2 / (sigma**2))


minfs, maxfs = log(1e-9), log(1e9)  # max fit slope for exponential/gaussian


def broken_gauss(p, x, dt_fix=False, fid=None):

    p = list(p)
    # check if we want to return only a fading/rising function
    if dt_fix is not False:
        x_peak = dt_fix

        if dt_fix > 0:  # rising only
            a1 = m2f(p[0])
            p[1] = np.clip(p[1], minfs, maxfs)
            b1 = exp(p[1])
            return a1 * gauss(x - x_peak, b1)
        else:  # fading only
            a2 = m2f(p[0])
            p[1] = np.clip(p[1], minfs, maxfs)
            b2 = exp(p[1])
            return a2 * np.exp(-(x - x_peak) / b2)

    # gaussian rise
    x_peak = np.clip(p[0], 0, max(x))
    a1 = m2f(p[1])
    p[2] = np.clip(p[2], minfs, maxfs)
    p[3] = np.clip(p[3], minfs, maxfs)
    b1 = exp(p[2])
    leftside = a1 * gauss(x - x_peak, b1)

    # exponential decay
    a2 = a1 * gauss(0, b1)
    b2 = exp(p[3])
    rightside = a2 * np.exp(-(x - x_peak) / b2)

    leftside[x > x_peak] = 0
    rightside[x <= x_peak] = 0

    return leftside + rightside


minmaxcs = 0.1  # min/max colorslope
mincolor, maxcolor = -1, +3


def color_func(x, x_peak, mean_color_in, slope_in):

    mean_color = np.clip(mean_color_in, mincolor, maxcolor)
    slope = np.clip(slope_in, -minmaxcs, minmaxcs)
    color = mean_color + (x - x_peak) * slope
    # ileft = x<=x_peak
    # if sum(ileft) and (x_peak>0):
    # 	color[ileft] = mean_color
    return color


def broken_gauss_twocolor(
    p, x, dt_fix=False, fid=None, flt1="g", flt2="r", x_peakc=None
):

    p = list(p)
    both = broken_gauss(p, x, dt_fix=dt_fix)

    if dt_fix is not False:
        if not x_peakc:
            x_peakc = np.median(x)  # new in 2019 (was 0)
        p[2] = np.clip(p[2], mincolor, maxcolor)
        # p[3] = np.clip(p[3], -minmaxcs, minmaxcs)
        if dt_fix > 0:
            p_col = p[2], 0
        else:
            p_col = p[2], p[3]
    else:
        if not x_peakc:
            x_peakc = np.median(x)  # new in 2019 np.clip(p[0], 0, max(x))
        p[4] = np.clip(p[4], mincolor, maxcolor)
        p[5] = np.clip(p[5], -minmaxcs, minmaxcs)
        p_col = p[4], p[5]

    iflt1 = fid == flt1
    iflt2 = fid == flt2

    color = color_func(x, x_peakc, p_col[0], p_col[1])

    # this can be tricky, when one band has very little data thing can go wrong
    both[iflt1] *= 10 ** (-0.4 * color[iflt1] / 2.0)
    both[iflt2] *= 10 ** (+0.4 * color[iflt2] / 2.0)

    return both


def res(p, x, y, yerr, model_func, dt_fix=False, fid=None, islim=None):

    model_flux = model_func(p, x, dt_fix=dt_fix, fid=fid)

    chi = (model_flux - y) / yerr
    chi = np.clip(chi, -1e5, 1e5)

    # do a rough prior on rise time, increase chi2 if larger/smaller than 100/0.01

    if (model_func is broken_gauss) or (model_func is broken_gauss_twocolor):

        p_col = 0, 0
        # check if we we have a fixed time of peak
        if dt_fix is False:
            logb1 = p[2]
            logb2 = p[3]
            if model_func is broken_gauss_twocolor:
                p_col = p[4], p[5]
        else:
            if dt_fix > 0:  # rising only
                logb1 = p[1]
                logb2 = 0
                if model_func is broken_gauss_twocolor:
                    p_col = p[2], 0
            else:  # fading only
                logb1 = 0
                logb2 = p[1]
                if model_func is broken_gauss_twocolor:
                    p_col = p[2], p[3]

        # log(100) = 4.60517
        rfmax = log(100)
        if abs(logb1) > rfmax:
            chi += (abs(logb1) - rfmax) / len(x)
        if abs(logb2) > 4.605:
            chi += (abs(logb2) - rfmax) / len(x)

        # weak push for constant color (helpful if limited data available?)
        chi += (p_col[1] - 0) / 0.05 / len(x)

    return chi


def flex_fit_wclip(
    dtime,
    flux,
    flux_err,
    fid,
    logger,
    isdetect=None,
    filters=["g", "r"],
    ax=None,
    verbose=True,
    niter=0,
):
    """
    run flex fit a few time for clipping
    also pick the preferred band for reporting the results
    """

    lc_data = {}
    plot_info = {}
    result, lsq, flux_diff = flex_fit(
        dtime,
        flux,
        flux_err,
        fid,
        logger,
        isdetect=isdetect,
        filters=["g", "r"],
        ax=ax,
        verbose=True,
    )
    # print('dtime:', dtime)
    (
        lc_data["dtime"],
        lc_data["flux"],
        lc_data["flux_err"],
        lc_data["fid"],
        lc_data["isdetect"],
    ) = (dtime, flux, flux_err, fid, isdetect)
    plot_info["lsq_out"] = lsq
    plot_info["dt_fix"] = result["dt_fix"]
    outliers = (abs(flux_diff / flux_err) > 7) * isdetect
    logger.info("flex_fit_wclip: # of outliers {0}".format(sum(outliers)))

    # ---
    # pick what band(s) to use store the final results
    final_result = {}
    flt_pick = None
    if ("color" in lsq) and (("r" in lsq) or ("g" in lsq)):
        flt_pick = "color"
    elif "r" in lsq:
        flt_pick = "r"
    elif "g" in lsq:
        flt_pick = "g"
    if flt_pick is not None:
        for item in result[flt_pick]:
            final_result[item] = result[flt_pick][item]
            final_result["band"] = flt_pick
    else:
        logger.info("flex_fit_wclip: no band yielded a succesful fit")

    if (sum(outliers) == 0) or niter > 5:
        plot_info["dtime_out"], plot_info["flux_out"] = None, None
        return final_result, lc_data, plot_info
    else:
        # print (flux_diff/flux_err)
        # print (flux_diff)
        logger.info(
            "flex_fit_wclip: running again after outlier rejection, attempts left={0}".format(
                5 - niter
            )
        )

        # allow one rejection
        maxout = max(abs(flux_diff[isdetect] / flux_err[isdetect]))
        outliers = abs(flux_diff / flux_err) == maxout
        logger.info(f"outliers with flux: {flux[outliers]}")

        plot_info["dtime_out"] = dtime[outliers]
        plot_info["flux_out"] = flux[outliers]

        dtime = dtime[outliers == False]
        flux = flux[outliers == False]
        flux_err = flux_err[outliers == False]
        fid = fid[outliers == False]
        isdetect = isdetect[outliers == False]

        if len(flux[isdetect]) < 4:
            logger.info("flex_fit_wclip: not enough data left, keeping original result")
            return final_result, lc_data, plot_info

        # loop continues
        return flex_fit_wclip(
            dtime,
            flux,
            flux_err,
            fid,
            logger,
            isdetect=isdetect,
            niter=niter + 1,
            filters=["g", "r"],
            ax=ax,
            verbose=verbose,
        )


def flex_fit(
    dtime,
    flux,
    flux_err,
    fid,
    logger,
    isdetect=None,
    filters=["g", "r"],
    ax=None,
    verbose=True,
):
    """
    function to fit a light curve with one or two bands with two Gaussian and get some color info
    returning scipy.optimize.leastq output dicts for each band and for the fit to both bands

    observations that are upper limits should have zero flux and be indicated with the isdetect=[] input (boolians array)
    """

    from .flex_metrics import flux2mag, mag2flux

    # unless otherwise indicated, asume all data is a detection
    if isdetect is None:
        isdetect = np.repeat(True, len(dtime))

    # output dicts
    dt_fix = {}
    result = {k: {} for k in filters + ["color"]}
    lsq_out = {}

    # do initial guess of parameters
    flux_peak_guess = flux[isdetect][
        np.argsort(flux[isdetect])[-1]
    ]  # pick brightest point as peak
    t_peak_guess = dtime[isdetect][np.argmax(flux[isdetect])]
    islim = isdetect == False

    flux_diff = np.zeros(len(flux))

    # irise = dtime<t_peak_guess0:
    # if sum(irise):
    sigma_rise_guess = np.clip(t_peak_guess - min(dtime[isdetect]), 1, 20)
    ipost = dtime > t_peak_guess
    if sum(ipost * isdetect):
        sigma_fade_guess = np.clip(
            np.interp(
                flux_peak_guess / 2.0,
                flux[ipost * isdetect][::-1],
                dtime[ipost * isdetect][::-1],
            )
            - t_peak_guess,
            1,
            40,
        )
    else:
        sigma_fade_guess = 20

    if verbose:
        logger.info(
            "guessed time of peak (wrt first obs) {0:0.2f}".format(t_peak_guess)
        )
        logger.info(
            "guessed peak (flux, mag)  {0:0.2f} {1:0.2f}".format(
                flux_peak_guess, flux2mag(flux_peak_guess)
            )
        )
        logger.info(
            "guessed sigma (rise, fade)  {0:0.2f} {1:0.2f}".format(
                sigma_rise_guess, sigma_fade_guess
            )
        )

    # do single band fit for the two filter in fid
    for k in filters:

        iflt = fid == k
        if verbose:
            logger.info(str(k) + " " + str(sum(iflt * isdetect)) + " " + str(sum(iflt)))

        if (sum(iflt * isdetect) > 2) and (sum(iflt) > 3):

            # do the least-square fit
            p0 = [
                t_peak_guess,
                f2m(flux_peak_guess),
                log(sigma_rise_guess),
                log(sigma_fade_guess),
            ]
            dt_fix[k] = False
            lsq_out[k] = leastsq(
                res,
                p0,
                (dtime[iflt], flux[iflt], flux_err[iflt], broken_gauss, dt_fix[k]),
                full_output=True,
            )

            result[k]["dtime_peak"] = lsq_out[k][0][0]
            result[k]["mag_peak"] = flux2mag(m2f(lsq_out[k][0][1]))
            result[k]["sigma_rise"] = exp(lsq_out[k][0][2])
            result[k]["sigma_fade"] = exp(lsq_out[k][0][3])

            if lsq_out[k][1] is not None:
                result[k]["e_dtime_peak"] = sqrt(lsq_out[k][1][0, 0])
                result[k]["e_mag_peak"] = sqrt(lsq_out[k][1][1, 1])
                result[k]["e_sigma_rise"] = (
                    sqrt(lsq_out[k][1][2, 2]) * result[k]["sigma_rise"]
                )
                result[k]["e_sigma_fade"] = (
                    sqrt(lsq_out[k][1][3, 3]) * result[k]["sigma_fade"]
                )

                if verbose:
                    ss = [
                        "{0:0.2f} ({1:0.2f}) ".format(result[k][x], result[k]["e_" + x])
                        for x in ("dtime_peak", "mag_peak", "sigma_rise", "sigma_fade")
                    ]
                    logger.info(f"rise/fade fit: {ss}")

            # if we dont get a decent fit with errorbars,
            # check that the source is not only risng/fading
            else:

                if verbose:
                    ss = [
                        "{0:0.2f} ".format(result[k][x])
                        for x in ("dtime_peak", "mag_peak", "sigma_rise", "sigma_fade")
                    ]
                    logger.info(f"rise/fade fit: {ss}")

                if lsq_out[k][0][0] > (max(dtime[iflt]) - 3):
                    dt_fix[k] = max(dtime[isdetect])
                    p0 = [f2m(flux_peak_guess), log(sigma_rise_guess)]
                elif lsq_out[k][0][0] < (min(dtime[iflt]) + 3):
                    dt_fix[k] = 0
                    p0 = [f2m(flux_peak_guess), log(sigma_fade_guess)]

            if dt_fix[k] is not False:
                if verbose:
                    logger.info(
                        "retrying fit without break, using x_peak fixed at {0:0.2f}".format(
                            dt_fix[k]
                        )
                    )
                lsq_out[k] = leastsq(
                    res,
                    p0,
                    (dtime[iflt], flux[iflt], flux_err[iflt], broken_gauss, dt_fix[k]),
                    full_output=True,
                )

                if verbose:
                    logger.info(
                        "dt_fix fit: {0:0.2f} {1:0.2f}".format(
                            lsq_out[k][0][0], lsq_out[k][0][1]
                        )
                    )

                # overwrite results with single fit
                result[k]["dtime_peak"] = dt_fix[k]
                result[k]["e_dtime_peak"] = -1.0
                result[k]["mag_peak"] = flux2mag(m2f(lsq_out[k][0][0]))
                if lsq_out[k][1] is not None:
                    result[k]["e_mag_peak"] = sqrt(lsq_out[k][1][0, 0])

                if dt_fix[k] > 0:
                    result[k]["sigma_rise"] = exp(lsq_out[k][0][1])
                    result[k]["sigma_fade"] = 0
                    if lsq_out[k][1] is not None:
                        result[k]["e_sigma_rise"] = (
                            sqrt(lsq_out[k][1][1, 1]) * result[k]["sigma_rise"]
                        )
                else:
                    result[k]["sigma_fade"] = exp(lsq_out[k][0][1])
                    result[k]["sigma_rise"] = 0
                    if lsq_out[k][1] is not None:
                        result[k]["e_sigma_fade"] = (
                            sqrt(lsq_out[k][1][1, 1]) * result[k]["sigma_fade"]
                        )

            result[k]["n_rise"] = int(
                sum(dtime[isdetect * iflt] < result[k]["dtime_peak"])
            )
            result[k]["n_fade"] = int(
                sum(dtime[isdetect * iflt] > result[k]["dtime_peak"])
            )

            flux_model = broken_gauss(lsq_out[k][0], dtime[iflt], dt_fix=dt_fix[k])
            flux_diff[iflt] = flux_model - flux[iflt]

            result[k]["mad"] = astropy.stats.median_absolute_deviation(
                flux_diff[isdetect]
            ) / np.median(flux[isdetect * iflt])
            result[k]["rms"] = np.std(flux_diff[isdetect]) / np.mean(
                flux[isdetect * iflt]
            )
            result[k]["chi2"] = sum(
                (flux_diff[iflt * isdetect] / flux_err[iflt * isdetect]) ** 2
            ) / (sum(isdetect) - len(p0))
            logger.info(
                "chi2={0:0.2f}; RMS={1:0.2f}, MAD={2:0.2f}".format(
                    result[k]["chi2"], result[k]["rms"], result[k]["mad"]
                )
            )

    # ----------
    # fit both colors at the same time
    iflt_0 = fid == filters[0]
    iflt_1 = fid == filters[1]

    # at least 6 detection are needed given the number of free parameters
    # plus we require at least two detections in a given band,
    # otherwise the fit has too much freedom to go wild due to the color evolution term
    if (
        sum(isdetect) > 5
        and (sum(isdetect * iflt_0) > 1)
        and (sum(isdetect * iflt_1) > 1)
    ):

        if verbose:
            logger.info(f"two color {sum(isdetect)}  {len(dtime)}")

        dt_fix["color"] = False
        p0 = [
            t_peak_guess,
            f2m(flux_peak_guess),
            log(sigma_rise_guess),
            log(sigma_fade_guess),
            -0.1,
            0.0,
        ]
        lsq_out["color"] = leastsq(
            res,
            p0,
            (dtime, flux, flux_err, broken_gauss_twocolor, dt_fix["color"], fid),
            full_output=True,
        )  # , maxfev=int(1e6), xtol=1e-99, ftol=1e-99)

        # store results
        result["color"]["dtime_peak"] = lsq_out["color"][0][0]
        result["color"]["mag_peak"] = flux2mag(m2f(lsq_out["color"][0][1]))
        result["color"]["sigma_rise"] = exp(lsq_out["color"][0][2])
        result["color"]["sigma_fade"] = exp(lsq_out["color"][0][3])
        result["color"]["mean_color"] = lsq_out["color"][0][4]
        result["color"]["color_slope"] = lsq_out["color"][0][5]

        if lsq_out["color"][1] is not None:
            result["color"]["e_dtime_peak"] = sqrt(lsq_out["color"][1][0, 0])
            result["color"]["e_mag_peak"] = sqrt(lsq_out["color"][1][1, 1])
            result["color"]["e_sigma_rise"] = (
                sqrt(lsq_out["color"][1][2, 2]) * result["color"]["sigma_rise"]
            )
            result["color"]["e_sigma_fade"] = (
                sqrt(lsq_out["color"][1][3, 3]) * result["color"]["sigma_fade"]
            )
            result["color"]["e_mean_color"] = sqrt(lsq_out["color"][1][4, 4])
            result["color"]["e_color_slope"] = sqrt(lsq_out["color"][1][5, 5])

            if verbose:
                ss = [
                    "{0:0.2f} ({1:0.2f}) ".format(
                        result["color"][x], result["color"]["e_" + x]
                    )
                    for x in (
                        "dtime_peak",
                        "mag_peak",
                        "sigma_rise",
                        "sigma_fade",
                        "mean_color",
                        "color_slope",
                    )
                ]
                logger.info(f"rise/fade fit: {ss}")

        # again, if we dont get a decent fit with errorbars,
        # check that the source is not only rising/fading
        else:

            if verbose:
                ss = [
                    "{0:0.2f}".format(result["color"][x])
                    for x in (
                        "dtime_peak",
                        "mag_peak",
                        "sigma_rise",
                        "sigma_fade",
                        "mean_color",
                        "color_slope",
                    )
                ]
                logger.info(f"rise/fade fit: {ss}")

            if lsq_out["color"][0][0] > (max(dtime) - 3):
                dt_fix["color"] = max(dtime[isdetect])
                p0 = [f2m(flux_peak_guess), log(sigma_rise_guess), 0, -0.02]
            elif lsq_out["color"][0][0] < (min(dtime) + 3):
                dt_fix["color"] = 0
                p0 = [f2m(flux_peak_guess), log(sigma_fade_guess), 0, -0.02]

        if dt_fix["color"] is not False:

            if verbose:
                logger.info(
                    "retrying fit without break, using x_peak fixed at {0:0.2f} ".format(
                        dt_fix["color"]
                    )
                )

            if dt_fix["color"] > 0:
                if verbose:
                    logger.info("no allowing color evolution")
                p0 = [f2m(flux_peak_guess), log(sigma_rise_guess), -0.1]
            else:
                p0 = [f2m(flux_peak_guess), log(sigma_rise_guess), -0.1, 0]

            lsq_out["color"] = leastsq(
                res,
                p0,
                (dtime, flux, flux_err, broken_gauss_twocolor, dt_fix["color"], fid),
                full_output=True,
            )

            if verbose:
                logger.info(
                    "two-color dt_fix: m_peak={0:0.2f} rise={1:0.2f}     | color={2:0.3f}".format(
                        flux2mag(m2f(lsq_out["color"][0][0])),
                        exp(lsq_out["color"][0][1]),
                        lsq_out["color"][0][2],
                    )
                )

            # overwrite results with single shape two color fit
            result["color"]["dtime_peak"] = dt_fix["color"]
            result["color"]["e_dtime_peak"] = -1
            result["color"]["mag_peak"] = flux2mag(m2f(lsq_out["color"][0][0]))
            result["color"]["mean_color"] = lsq_out["color"][0][2]
            if dt_fix["color"] > 0:
                result["color"]["color_slope"] = 0
            else:
                result["color"]["color_slope"] = lsq_out["color"][0][3]

            if lsq_out["color"][1] is not None:
                result["color"]["e_mag_peak"] = sqrt(lsq_out["color"][1][1, 1])
                result["color"]["e_mean_color"] = sqrt(lsq_out["color"][1][2, 2])
                if dt_fix["color"] > 0:
                    result["color"]["e_color_slope"] = -999
                else:
                    result["color"]["e_color_slope"] = sqrt(lsq_out["color"][1][3, 3])

            if dt_fix["color"] > 0:
                result["color"]["sigma_rise"] = exp(lsq_out["color"][0][1])
                result["color"]["sigma_fade"] = 0
                if lsq_out["color"][1] is not None:
                    result["color"]["e_sigma_rise"] = (
                        sqrt(lsq_out["color"][1][1, 1]) * result["color"]["sigma_rise"]
                    )
            else:
                result["color"]["sigma_fade"] = exp(lsq_out["color"][0][1])
                result["color"]["sigma_rise"] = 0
                if lsq_out["color"][1] is not None:
                    result["color"]["e_sigma_fade"] = (
                        sqrt(lsq_out["color"][1][1, 1]) * result["color"]["sigma_fade"]
                    )

        flux_model = broken_gauss_twocolor(
            lsq_out["color"][0], dtime, dt_fix=dt_fix["color"], fid=fid
        )
        flux_diff = flux_model - flux

        result["color"]["mad"] = astropy.stats.median_absolute_deviation(
            flux_diff[isdetect]
        ) / np.median(flux[isdetect])
        result["color"]["rms"] = np.std(flux_diff[isdetect]) / np.mean(flux[isdetect])
        result["color"]["chi2"] = sum(
            (flux_diff[isdetect] / flux_err[isdetect]) ** 2
        ) / (sum(isdetect) - len(p0))
        logger.info(
            "chi2={0:0.2f}; RMS={1:0.2f}, MAD={2:0.2f}".format(
                result["color"]["chi2"], result["color"]["rms"], result["color"]["mad"]
            )
        )

        result["color"]["n_rise"] = int(
            sum(dtime[isdetect] <= result["color"]["dtime_peak"])
        )
        result["color"]["n_fade"] = int(
            sum(dtime[isdetect] >= result["color"]["dtime_peak"])
        )

    result["dt_fix"] = dt_fix

    return result, lsq_out, flux_diff


# make some higher-level statements about the light curve
def flex_class(result):

    photoclass = ""
    chat_bot = ""

    # check that we have any results (as selected by flex_fit_wclip)
    if not "band" in result:
        return "", "no fit"

    # print color information
    if (result["band"] == "color") and (
        (result["n_rise"] > 2) or (result["n_fade"] > 2)
    ):

        if "e_mean_color" in result:
            if abs(result["e_mean_color"]) < 0.1:
                chat_bot += "<g-r>={0:0.2f}+-{1:0.2f}, ".format(
                    result["mean_color"], np.clip(result["e_mean_color"], 0, 999)
                )

    # print info if we are on the rise
    if result["n_rise"] > 1:
        if result["sigma_rise"] < 300:
            if result["n_fade"] == 0:
                chat_bot += "rising, "
            elif result["n_fade"] < 2:
                chat_bot += "at/near peak, "

            # check slope of rise time for SNe or CV like behaviour
            if not ("e_sigma_rise" in result):
                result["e_sigma_rise"] = 999  # fix for printing

            if not ("e_sigma_fade" in result):
                result["e_sigma_fade"] = 999

            # slow rise time is quite rare, can be AGN, TDE, SN II
            if result["sigma_rise"] > 10:
                photoclass = "TDE?"

                chat_bot += "slow rise time (sigma={0:0.1f}+-{1:0.1f} d), ".format(
                    result["sigma_rise"], np.clip(result["e_sigma_rise"], 0, 999)
                )

                # if also blue this could be a TDE
                if (result["band"] == "color") and (result["mean_color"] < 0.0):
                    photoclass = "TDE"
            elif result["sigma_rise"] > 4:
                chat_bot += "SN-like rise time (sigma={0:0.1f}+-{1:0.1f} d), ".format(
                    result["sigma_rise"], np.clip(result["e_sigma_rise"], 0, 999)
                )
            # 	photoclass = 'SN?'

        else:
            chat_bot += "very slow rise, "
            if (result["n_fade"] > 2) and (result["sigma_fade"] > 500):
                photoclass = "AGN?"

    if (
        (result["n_rise"] > 1)
        and (result["sigma_rise"] < 2)
        and ("e_sigma_rise" in result)
    ):
        chat_bot += "warning, steep increase (sigma={0:0.1f}+-{1:0.1f}), ".format(
            result["sigma_rise"], np.clip(result["e_sigma_rise"], 0, 999)
        )
        photoclass = "CV?"

    if (result["n_rise"] < 2) and (result["sigma_fade"] < 500):
        chat_bot += "missed peak(?), "

    # print info for fading or constant sources
    if (result["n_fade"] > 2) and ("e_sigma_fade" in result):

        if result["sigma_fade"] < 1000:

            chat_bot += "fading (tau={0:0.1f} +- {1:0.1f} d), ".format(
                result["sigma_fade"], np.clip(result["e_sigma_fade"], 0, 999)
            )

        elif result["n_fade"] > 3:

            chat_bot += "very long fading time, "

            # if fading time is actually measured properly, do classification
            if "e_sigma_fade" in result:
                if (result["sigma_fade"] / result["e_sigma_fade"]) > 2:
                    if not photoclass:
                        photoclass = "AGN?"
                elif "color_slope" in result:
                    if abs(result["color_slope"]) < 0.015:
                        photoclass = "not SN?"

    # do checks on color and slope
    if (result["band"] == "color") and ("e_mean_color" in result):

        # add label if blue
        if (result["mean_color"] < -0.2) and (abs(result["e_mean_color"]) < 0.1):
            chat_bot += "blue! ({0:0.1f}+-{1:0.1f}), ".format(
                result["mean_color"], np.clip(result["e_mean_color"], 0, 999)
            )

        # type as TDE if constant color and blue (since May 2019 this now includes ~constant post-peak flux, plus cut on e_color_slope 0.015)
        if result["e_color_slope"] < 0.015:

            if result["color_slope"] < -0.015:
                chat_bot += "getting more blue?! ({0:0.3f}+-{1:0.3f}), ".format(
                    result["color_slope"], np.clip(result["e_color_slope"], 0, 999)
                )
                if (result["mean_color"] < 0) and (result["sigma_fade"] > 5):
                    photoclass = "TDE"
                elif result["mean_color"] < 0.1:
                    photoclass = "TDE?"

            elif abs(result["color_slope"]) < 0.015:
                if result["e_color_slope"] > 0:
                    chat_bot += "near-constant color ({0:0.3f}+-{1:0.3f}) ".format(
                        result["color_slope"], np.clip(result["e_color_slope"], 0, 999)
                    )
                if (result["mean_color"] < 0) and (result["sigma_fade"] > 5):
                    photoclass = "TDE"
                elif result["mean_color"] < 0.1:
                    photoclass = "TDE?"

            elif result["color_slope"] < 0.020:
                chat_bot += "cooling? ({0:0.3f}+-{1:0.3f}), ".format(
                    result["color_slope"], np.clip(result["e_color_slope"], 0, 999)
                )
                if (result["mean_color"] < -0.1) and (result["sigma_fade"] > 5):
                    photoclass = "TDE"
                elif result["mean_color"] < 0.0:
                    photoclass = "TDE?"

            elif (result["color_slope"] >= 0.020) and (result["sigma_fade"] < 100):
                chat_bot += "cooling ({0:0.3f}+-{1:0.3f}), ".format(
                    result["color_slope"], np.clip(result["e_color_slope"], 0, 999)
                )
                if result["color_slope"] != minmaxcs:
                    photoclass = "SN"

    # check fast fading
    if (result["n_fade"] > 2) and (result["sigma_fade"] < 5):
        chat_bot += "warning fast fading, "
        photoclass = "fast/weird"

    if not photoclass:
        photoclass = "unknown"

    if result["chi2"] > 10:
        chat_bot += "warning, poor fit (chi2/dof={0:0.1f})".format(result["chi2"])
        # photoclass +='_badfit'

    return chat_bot.strip(", "), photoclass
