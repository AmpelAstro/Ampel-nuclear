#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t0/NuclearFilter.py
# License           : BSD-3-Clause
# Author            : sjoertvv <sjoert@umd.edu>
# Date              : 26.02.2018
# Last Modified Date: 29.11.2022
# Last Modified By  : sr <simeon.reusch@desy.de>

from typing import Optional, Union, List, Dict, Any

import numpy as np
import astropy.coordinates  # need for l, b computation

from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.abstract.AbsAlertFilter import AbsAlertFilter

from ampel.nuclear.t0.GaiaVeto import GaiaVetoMixin

REJECTION_REASON_CODES = {
    999: "Passed",
    10: "Large sgscore",
    11: "Large sgscore2",
    12: "No sgscore available",
    20: "Bad PS1 photometry",
    21: "PS1 too far",
    22: "Bright 1st PS1 match",
    23: "Bright 2nd PS1 match",
    24: "Bright 3rd PS1 match",
    30: "Nothing passed default filter",
    40: "Not enough detections",
    50: "Too bright in ZTF reference image",
    60: "Too many PS1 matches (field too crowded)",
    61: "Too close to Galactic Plane (field too crowded)",
    70: "RB too low",
    80: "Only duplicate detections",
    90: "Potential mover",
    100: "Large gap between detections",
    110: "Too faint",
    120: "Flux increase too small",
    130: "Gaia veto",
    140: "Distance veto",
    141: "No distance could be computed",
    150: "Last detection did not pass default filter (only when lastOnly=True)",
}


class NuclearFilter(AbsAlertFilter, GaiaVetoMixin):
    """
    Filter for a clean sample of nuclear transients in galaxies
    includes a veto of PS1 star-galaxy score using Gaia data
    """

    maxDeltaRad: float = 0.5  #: the max flare-ref distance to be accepted, we try mulitple ways of summing distnr for multiple detection
    minDeltaRad: float = 0.0  #: the min flare-ref distance to be accepted

    maxSgscore: float = 0.3  #: min star-galaxy score

    minDeltaJD: float = 0.01  #: remove movers with min time distance in days between detections that pass default filter

    maxDeltaJD: float = 10.0  #: remove lightcurves with gaps, we require we have at least two detections seperated by this distance

    minDetections: int = 3  #: number of detections in each band

    diffmagLimit: float = 20  #: magnitude limit applied to difference flux (magpsf)

    maxDeltaMag: float = 2.5  #: lower limit to flux increase: PSF magnitude in difference image minus PSF magnitude in the references

    minRealBogusScore: float = 0.3  #: min RealBogus score of *any* observation

    maxDistPS1source: float = 1  #: min distance to nearest PS1 source (useful for removing SNe with detections in the ref image, bright stars, and ghostly things)

    closePS1dist: float = 0.5  #: distance when we consider the second PS1 match into our checks for the sgscore and photometry (0.5" recommended)

    brightRefMag: float = 12.0  #: bright star removal, used for both ZTF filters

    brightPS1gMag: float = 12.0  #: PS1 bright star removal
    brightPS1rMag: float = 12.0  #: PS1 bright star removal
    brightPS1iMag: float = 12.0  #: PS1 bright star removal
    brightPS1zMag: float = 12.0  #: PS1 bright star removal

    removeNegative: int = 1  #: remove negative detections from lightcurve (or not)

    minGalLat: float = 5.0  #: you may stay away from the Galactic plane

    doGaiaVeto: int = 1  #: triggers the use of the Gaia-based filter for the sgsore, designed to remove star that get a high score

    maxPS1matches: int = 100  #: remove source in dense fields (but some galaxies also get many matches, so leave this large!)

    lastOnly: int = (
        0  #: apply all cuts only to the most recent detection in the light curve
    )

    """
    Return values:

    True: Passed all cuts
    
    -10: Large sgscore
    -11: Large sgscore2
    -12: No sgscore available

    -20: Bad PS1 photometry
    -21: PS1 too far
    -22: Bright 1st PS1 match
    -23: Bright 2nd PS1 match
    -24: Birght 3rd PS1 match

    -30: Nothing passed default filter

    -40: Not enough detections

    -50: Too bright in ZTF reference image

    -60: Too many PS1 matches (field too crowded)
    -61: Too close to Galactic Plane (field too crowded)

    -70: RB too low

    -80: Only duplicate detections

    -90: Potential mover

    -100: Large gap between detections

    -110: Too faint

    -120: Flux increase too small

    -130: Gaia veto

    -140: Distance veto
    -141: No distance could be computed

    -150: Last detection did not pass default filter (only when lastOnly=True)
    """

    def post_init(self) -> None:

        # collect bright cuts into dict
        self.brightPS1Mag = {}
        self.brightPS1Mag["g"] = self.brightPS1gMag
        self.brightPS1Mag["r"] = self.brightPS1rMag
        self.brightPS1Mag["i"] = self.brightPS1iMag
        self.brightPS1Mag["z"] = self.brightPS1zMag

        # Instance dictionaries later used in method apply

        # remove upper limits
        isdetect_flt = {"attribute": "candid", "operator": "is not", "value": None}

        # no bad pixels
        nbad_flt = {"attribute": "nbad", "operator": "<=", "value": 0}

        # similar to Growth Marshal mandatory filter (Full Width Half Max assuming a Gaussian core, from SExtractor)
        fwhm_flt = {"attribute": "fwhm", "operator": "<", "value": 5}

        # similar to Growth Marshal mandatory filter (Ratio: aimage / bimage)
        elong_flt = {"attribute": "elong", "operator": "<", "value": 1.3}

        # similar to Growth Marshal mandatory filter (Difference: magap - magpsf )
        magdiff_flt_lower = {"attribute": "magdiff", "operator": ">", "value": -0.3}

        magdiff_flt_upper = {"attribute": "magdiff", "operator": "<", "value": 0.3}

        # has host galaxy detected (removes orphans)
        distnr_flt = {"attribute": "distnr", "operator": ">", "value": 0}

        # combine the current filters
        self._default_filters: List[Dict[str, Any]] = [
            isdetect_flt,
            distnr_flt,
            nbad_flt,
            fwhm_flt,
            elong_flt,
            magdiff_flt_lower,
            magdiff_flt_upper,
        ]

        # remove negative detection (or not)
        if self.removeNegative:

            # is not ref-science
            isdiff_flt1 = {"attribute": "isdiffpos", "operator": "!=", "value": "0"}

            # is not ref-science (again!)
            isdiff_flt2 = {"attribute": "isdiffpos", "operator": "!=", "value": "f"}

            self._default_filters += [isdiff_flt1, isdiff_flt2]

        self.logger.info("NuclearFilter initialized")

    def wmean(self, x, w):
        """
        compute a weighted mean, used for offset statistics
        """
        return np.sum(x * w) / np.sum(w)

    def process(self, alert: AmpelAlertProtocol) -> Optional[Union[bool, int]]:
        """
        Mandatory implementation.
        To exclude the alert, return *None*
        To accept it, either return
        * self.on_match_t2_units
        * or a custom combination of T2 unit names

        Returns None (reject) or True (accept)

        Make a selection on:
        - the distance between the transient and host in reference image
        - the Real/Bogus sore
        - the distance to a bright star

        """

        # first check we have an extended source (note this can remove flares from faint galaxies that missclassified in PS1)
        # these will have to be dealt with in the orphan/faint filter

        if len(alert.get_values("sgscore1")) == 1:

            alert_pps = alert.datapoints[0]
            sgscore1, sgscore2, sgscore3 = (
                alert_pps["sgscore1"],
                alert_pps["sgscore2"],
                alert_pps["sgscore3"],
            )
            distpsnr1, distpsnr2, distpsnr3 = (
                alert_pps["distpsnr1"],
                alert_pps["distpsnr2"],
                alert_pps["distpsnr3"],
            )

        # check that we dont have very old (pre v1.8) schema
        else:

            self.why = "this schema doesnt have sgscore1, rejecting this alert"
            self.logger.info(None, extra={"sgscore1": None})
            return -12

        # print a warning if we have PS1 match issues
        if distpsnr2 < self.closePS1dist:
            self.logger.info(
                "Potential problem with PS1, we have close matches",
                extra={"distpsnr1": distpsnr1, "distpsnr2": distpsnr2},
            )

        # check sgscore, allowing for close PS1 matches with a safescore to overwrite the nearest score
        if (sgscore1 > self.maxSgscore) and not (
            (sgscore2 < self.minSafeSgscore) and (abs(distpsnr2) < self.closePS1dist)
        ):

            self.why = "sgscore1 = {0:0.2f}, which is > {1:0.2f}".format(
                sgscore1, self.maxSgscore
            )
            logextra = {"sgscore1": sgscore1, "maxSgscore": self.maxSgscore}
            if distpsnr2 < self.closePS1dist:
                self.why += "; sgscore2 = {0:0.2f} (dist={1:0.2f})".format(
                    sgscore2, distpsnr2
                )
                logextra["sgscore2"] = sgscore2
                logextra["distpsnr2"] = distpsnr2

            self.logger.info("Large sgscore", extra=logextra)
            return -10

        # also check the sgscore of the close match
        if (sgscore2 > self.maxSgscore) and (abs(distpsnr2) < self.closePS1dist):
            self.why = (
                "sgscore2 = {0:0.2f}, which is > {1:0.2f}  (distpsnr2={2:0.2f})".format(
                    sgscore2, self.maxSgscore, distpsnr2
                )
            )
            self.logger.info(
                "Large sgscore2", extra={"sgscore2": sgscore2, "distpsnr2": distpsnr2}
            )
            return -11

        # check problems with photometry (likely for saturated stars)

        # collect all PS1 magnitude into a dict
        psmag1 = {k: alert_pps["s" + k + "mag1"] for k in ("g", "r", "i", "z")}
        psmag2 = {k: alert_pps["s" + k + "mag2"] for k in ("g", "r", "i", "z")}
        psmag3 = {k: alert_pps["s" + k + "mag3"] for k in ("g", "r", "i", "z")}

        if (psmag1["r"] < 0) and not (
            (psmag2["r"] > 0) and (abs(distpsnr2) < self.closePS1dist)
        ):
            self.why = "1st PS1 match r-band photometry is faulty: mag={0:0.2f} (dist={1:0.2f})".format(
                psmag1["r"], distpsnr1
            )
            if distpsnr2 < self.closePS1dist:
                self.why += "; 2nd PS match r-band also faulty: mag={0:0.2f} (dist={1:0.2f})".format(
                    psmag1["r"], distpsnr1
                )
            self.logger.info(
                "Bad PS1 photometry",
                extra={
                    "psmag1r": psmag1["r"],
                    "distpsnr1": distpsnr1,
                    "psmag2r": psmag2["r"],
                    "distpsnr2": distpsnr2,
                },
            )
            return -20

        # check that the nearest PS1 source is not too far
        if abs(distpsnr1) > self.maxDistPS1source:
            self.why = (
                "distance to 1st PS match is {0:0.2f}, which is > {1:0.2f}".format(
                    distpsnr1, self.maxDistPS1source
                )
            )
            self.logger.info("PS too far", extra={"distpsnr1": distpsnr1})
            return -21

        # check for nearby bright stars
        for k in ("g", "r", "i", "z"):

            if abs(psmag1[k]) < self.brightPS1Mag[k]:
                self.why = "1st PS match in {0}-band, mag={1:0.2f}, which is < {2:0.2f} (dist={3:0.2f})".format(
                    k, psmag1[k], self.brightPS1Mag[k], distpsnr1
                )
                self.logger.info(
                    "Bright 1st PS match",
                    extra={"filter": k, "psmag1": psmag1[k], "distpsnr1": distpsnr1},
                )
                return -22

            if psmag2 is not None:  # old schema check
                if (abs(psmag2[k]) < self.brightPS1Mag[k]) and (
                    abs(distpsnr2) < self.brightObjDist
                ):
                    self.why = "2nd PS match in {0}-band, mag={1:0.2f}, which is < {2:0.2f} (dist={3:0.2f})".format(
                        k, psmag2[k], self.brightPS1Mag[k], distpsnr2
                    )
                    self.logger.info(
                        "Bright 2nd PS match",
                        extra={
                            "filter": k,
                            "psmag2": psmag2[k],
                            "distpsnr2": distpsnr2,
                        },
                    )
                    return -23

            if psmag3 is not None:
                if (abs(psmag3[k]) < self.brightPS1Mag[k]) and (
                    abs(distpsnr3) < self.brightObjDist
                ):
                    self.why = "3rd  PS1 match in {0}-band, mag={1:0.2f}, which is < {2:0.2f} (dist={3:0.2f})".format(
                        k, psmag3[k], self.brightPS1Mag[k], distpsnr3
                    )
                    self.logger.info(
                        "Bright 3rd PS match",
                        extra={
                            "filter": k,
                            "psmag3": psmag3[k],
                            "distpsnr3": distpsnr3,
                        },
                    )
                    return -24

            # don't use the code below because it will remove sources next to objects
            # that were detected in just one pan-starrs band and thus have srmag=-999
            #
            # if ((psmag2['r']<0) or (psmag['g']<0)) and (abs(distpsnr2)< self.brightObjDist):
            #   self.why = "2nd PS1 match saturated(?) sgmag={0:0.2f} srmag={1:0.2f} (dist={2:0.2f})".format(psmag2['g'], psmag2['r'], distpsnr2)
            #   self.logger.info(self.why)
            #   return 0

        # get RealBogus scores for observations, check number of bad pixels
        res_tuples = alert.get_ntuples(
            ["rb", "jd", "magnr"], filters=self._default_filters
        )

        # check that we have anything
        if len(res_tuples) == 0:
            self.why = "nothing passed default filter"
            self.logger.info(self.why)
            return -30

        res = list(map(list, zip(*res_tuples)))
        rb_arr, jd_arr, magnr_arr = np.asarray(res)

        # check number of detections in all bands
        if len(jd_arr) < self.minDetections:
            self.why = "only {0} detections pass default filter, but we require at least {1}".format(
                len(jd_arr), self.minDetections
            )
            self.logger.info("not enough detections", extra={"detections": len(jd_arr)})

            return -40

        minmag = np.min(magnr_arr).item()
        # check that source is not too bright in ZTF ref img
        if self.brightRefMag > minmag:
            self.why = f"min(magnr)={minmag:0.2f}, which is < {self.brightRefMag:0.1f}"
            self.logger.info("too bright in ZTF ref img", extra={"minmag": minmag})

            return -50

        # check source density --> note, we also check Gaia source density, because this doesn't work well for big galaxies
        nmtchps = alert.get_values("nmtchps")[0]
        if nmtchps > self.maxPS1matches:
            self.why = "nmtchps={0}, which is > {1}".format(nmtchps, self.maxPS1matches)
            self.logger.info("nmtchps too large", extra={"nmtchps": nmtchps})
            return -60

        # check Galactic Latitude
        ra, dec = alert.get_values("ra")[0], alert.get_values("dec")[0]

        if ra > -99:
            gc = astropy.coordinates.SkyCoord(
                ra=ra * astropy.units.deg, dec=dec * astropy.units.deg
            ).galactic
            gal_l, gal_b = gc.galactic.l.deg, gc.galactic.b.deg
            if abs(gal_b) < self.minGalLat:
                self.logger.info(
                    "too close to Galactic plane",
                    extra={"gal_l": gal_l, "gal_b": gal_b},
                )

                return -61
        else:
            self.logger.warn(
                "problems with alert coordindates", extra={"ra": ra, "dec": dec}
            )

        # if we want, only check last observation
        if self.lastOnly:

            lastcheck = alert.get_values(
                "jd",
                filters=self._default_filters
                + [
                    {
                        "attribute": "jd",
                        "operator": "==",
                        "value": max(alert.get_values("jd")),
                    }
                ],
            )

            if len(lastcheck) == 0:
                self.why = "last detection did not pass default filter"
                self.logger.info(self.why)
                return -150

            rb_arr = [
                rb_arr[np.argmax(jd_arr)]
            ]  # make sure rb check below is only for last detection

        # if no detections pass real bogus, remove
        if max(rb_arr) < self.minRealBogusScore:
            self.why = "max(rb)={0:0.2f}, which is  < {1:0.2f}".format(
                max(rb_arr), self.minRealBogusScore
            )
            self.logger.info("max(rb) too low", extra={"max_rb": max(rb_arr)})

            return -70

        # do cut on moving sources (with all good detections)
        dt = abs(np.sort(jd_arr)[1:] - np.sort(jd_arr)[0:-1])

        # first check that we dont have bunch of duplicates
        if not sum(dt > 0):
            self.why = "number of detections={0}, but time difference between all is zero".format(
                len(dt)
            )
            self.logger.info(
                "only duplicate detections", extra={"len(dt)": int(len(dt))}
            )

            return -80

        # check that time between detection is not too short (as expected for moving object)
        dt = dt[dt > 0]
        if np.max(dt) < self.minDeltaJD:
            self.why = "potential mover, number of good detections={0}; \
                        max(time diff)={1:1.3f} h, which is < {2:0.3f} h".format(
                len(jd_arr), max(dt) * 24, self.minDeltaJD * 24
            )
            self.logger.info(
                "potential mover",
                extra={"len_jd_arr": len(jd_arr), "max_dt_hour": max(dt) * 24},
            )

            return -90

        # Check time between detecetions (this could be made more fancy using upper limits etc.)
        if np.min(dt) > self.maxDeltaJD:
            self.why = "number of good detections={0}; min(time diff)={1:1.1f} h, \
                        which is > {2:0.1f} h".format(
                len(jd_arr), min(dt) * 24, self.maxDeltaJD * 24
            )
            self.logger.info(
                "large gap between detection",
                extra={"len_jd_arr": len(jd_arr), "min_dt_hour": min(dt) * 24},
            )

            return -100

        # get some more arrays
        res = list(
            map(
                list,
                zip(
                    *alert.get_ntuples(
                        ["distnr", "magpsf", "sigmapsf", "rb", "fwhm", "fid"],
                        filters=self._default_filters,
                    )
                ),
            )
        )

        distnr_arr, magpsf_arr, sigmapsf_arr, rb_arr, fwhm_arr, fid_arr = np.asarray(
            res
        )

        # get indices to bands
        idx_g = fid_arr == 1
        idx_r = fid_arr == 2

        # removed check on detections per filter in Feb 2019; it's too strict, we can deal with this at the light curve fitting stage
        # if (sum(idx_g)<self.minDetections) and (sum(idx_r)<self.minDetections):
        #   self.why = "number of good detections in (g, r)=({0}, {1}), which is < {2:0}".format(sum(idx_g), sum(idx_r), self.minDetections)
        #   self.logger.info("not enough good detection",
        #       extra={"sum_idx_g":int(sum(idx_g)), "sum_idx_r":int(sum(idx_r))}
        #   )

        #   return 0

        # apply magnitude cut to brightest detection
        if np.min(magpsf_arr) > self.diffmagLimit:
            self.why = "too faint, min(mag)={0:0.2f}, which is \
                        > {1:0.1f}".format(
                np.min(magpsf_arr), self.diffmagLimit
            )
            self.logger.info("too faint", extra={"min_magpsf_arr": np.min(magpsf_arr)})

            return -110

        # apply lower limit on flux increase
        if np.min(magpsf_arr - magnr_arr) > self.maxDeltaMag:
            self.why = "flux increase too small, min(magpsf-magnr)={0:0.2f}, \
                        which is > {1:0.2f}".format(
                np.min(magpsf_arr - magnr_arr), self.maxDeltaMag
            )
            self.logger.info(
                "flux increase too small",
                extra={"min_magpsfnr_diff": np.min(magpsf_arr - magnr_arr)},
            )
            return -120

        # allow Gaia to veto
        if self.doGaiaVeto:
            G_veto = self.gaia_veto(alert)
            if G_veto:
                return -130

        # if we make it this far, compute the host-flare distance, using only (decent-enough) detections
        res = list(
            map(
                list,
                zip(
                    *alert.get_ntuples(
                        ["ra", "dec", "ranr", "decnr"], filters=self._default_filters
                    )
                ),
            )
        )

        ra_arr, dec_arr, ranr_arr, decnr_arr = np.asarray(res)

        # compute a few different measures of the distance
        # we compute these for each band seperately
        sigma_offset = np.clip(
            0.24 + 0.04 * (magpsf_arr - 20), 0.1, 1
        )  # equation from the NedStark paper
        snr_weight = 1 / sigma_offset**2

        offset_info = {}

        for idx, bnd in zip([idx_g, idx_r], ["g", "r"]):

            if sum(idx) >= self.minDetections:

                mean_dist = 3600 * np.sqrt(
                    np.mean((ra_arr - ranr_arr)[idx]) ** 2
                    + np.mean((dec_arr - decnr_arr)[idx]) ** 2
                )
                median_dist = 3600 * np.sqrt(
                    np.median((ra_arr - ranr_arr)[idx]) ** 2
                    + np.median((dec_arr - decnr_arr)[idx]) ** 2
                )
                weighted_dist = 3600 * np.sqrt(
                    self.wmean((ra_arr - ranr_arr)[idx], snr_weight[idx]) ** 2
                    + self.wmean((dec_arr - decnr_arr)[idx], snr_weight[idx]) ** 2
                )

                offset_info[bnd] = (mean_dist, median_dist, weighted_dist)

                if (mean_dist < self.maxDeltaRad) and (mean_dist > self.minDeltaRad):
                    self.why = "pass on mean dist={0:0.2f}, band={1}; detections used={2}".format(
                        mean_dist, bnd, sum(idx)
                    )

                    self.logger.info(
                        "pass on mean dist",
                        extra={
                            "mean_dist": mean_dist,
                            "band": bnd,
                            "detections": int(sum(idx)),
                        },
                    )
                    return True

                if (median_dist < self.maxDeltaRad) and (
                    median_dist > self.minDeltaRad
                ):
                    self.why = "pass on median dist={0:0.2f}, band={1}; detections used={2}".format(
                        median_dist, bnd, sum(idx)
                    )

                    self.logger.info(
                        "pass on median dist",
                        extra={
                            "medianDist": median_dist,
                            "band": bnd,
                            "detections": int(sum(idx)),
                        },
                    )

                    return True

                if (weighted_dist < self.maxDeltaRad) and (
                    weighted_dist > self.minDeltaRad
                ):
                    self.why = "pass on weighted dist={0:0.2f}, band={1}; detections used={2}".format(
                        weighted_dist, bnd, sum(idx)
                    )
                    self.logger.info(
                        "pass on weighted dist",
                        extra={
                            "weightedDist": weighted_dist,
                            "band": bnd,
                            "detections": int(sum(idx)),
                        },
                    )

                    return True

        # if none of the measures of the host-flare pass the cut, reject this alert
        if len(offset_info.keys()):
            self.why = ""
            for bnd in offset_info.keys():
                self.why += "mean/median/weighted dist in {0}-band = ({1:0.2f}/{2:0.2f}/{3:0.2f}), which is > {4:0.2f}".format(
                    bnd,
                    offset_info[bnd][0],
                    offset_info[bnd][1],
                    offset_info[bnd][2],
                    self.maxDeltaRad,
                )

                self.logger.info(
                    None,
                    extra={
                        "meanDist-" + bnd: mean_dist,
                        "weightedDist-" + bnd: weighted_dist,
                        "medianDist-" + bnd: median_dist,
                        "band": bnd,
                    },
                )

            return -140
        else:
            self.why = "offset distance not created, likely not enough data passed default filter"
            self.logger.info(self.why, extra={"sumIdx": int(sum(idx))})
            return -141
