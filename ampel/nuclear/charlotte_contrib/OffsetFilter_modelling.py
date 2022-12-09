#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/ztfbh/t0/NuclearFilter.py
# License           : BSD-3-Clause
# Author            : sjoertvv <sjoert@umd.edu>
# Date              : 26.02.2018
# Last Modified Date: 02.07.2018
# Last Modified By  : svv <sjoert@umd.edu>

import numpy as np
import astropy.coordinates  # need for l, b computation
import catsHTM
from math import radians
from urllib.parse import urlparse
from qso_fit import qso_fit
from ampel.abstract.AbsAlertFilter import AbsAlertFilter
from ampel.logging.AmpelLogger import AmpelLogger


class OffsetFilter(AbsAlertFilter):
    """ """

    # Static version info
    version = 1.0
    resources = ("catsHTM.default",)

    def __init__(
        self, on_match_t2_units, base_config=None, run_config=None, logger=None
    ):
        """ """
        if run_config is None or len(run_config) == 0:
            raise ValueError("Please check you run configuration")

        self.logger = AmpelLogger.get_logger() if logger is None else logger
        self.on_match_t2_units = on_match_t2_units

        # the max flare-ref distance (distnr) to be accepted, we try mulitple ways of summing distnr for multiple detection
        self.maxDeltaRad = run_config["MaxDeltaRad"]

        # the min flare-ref distance (distnr) to be accepted
        self.minDeltaRad = run_config["MinDeltaRad"]

        # star-galaxy score
        self.minSgscore = run_config["MinSgscore"]
        # star-galaxy score where we dont allow gaia veto based on photometry
        self.minSafeSgscore = run_config["MinSafeSgscore"]

        # remove movers with min time distance in days between detections that pass default filter
        self.minDeltaJD = run_config["MinDeltaJD"]

        # option to require we have (at least one epoch with) two detections within a night (or even shorter timescale)
        self.maxDeltaJD = run_config["MaxDeltaJD"]

        # number of detections in each band
        self.minDetections = run_config["MinDetections"]

        # magnitude limit applied to difference flux (magpsf)
        self.diffmagLimit = run_config["DiffmagLimit"]

        # lower limit to flux increase: PSF magnitude in difference image minus PSF magnitude in the references
        self.maxDeltaMag = run_config["MaxDeltaMag"]

        # min RealBogus score of *any* observation
        self.minRealBogusScore = run_config["MinRealBogusScore"]

        # min distance to nearest PS1 source
        # (useful for removing SNe with detections in the ref image, bright stars, and ghostly things)
        self.minDistPS1source = run_config["MinDistPS1source"]

        # distance when we consider the second PS1 match into our checks for the sgscore and photometry (not a free parameter)
        self.closePS1dist = 0.5

        # bright star removal: used for both ZTF filters
        self.brightRefMag = run_config["BrightRefMag"]

        # max distance for checking nearby bright stars in Gaia or PS1
        self.brightObjDist = run_config["BrightObjDist"]

        # PS1 bright star removal in mulitple bands
        self.brightPS1Mag = {
            k: run_config["BrightPS1" + k + "Mag"] for k in ["g", "r", "i", "z"]
        }

        # Gaia bright star removal in mulitple bands
        self.brightGaiaMag = run_config["BrightGaiaMag"]

        # remove source in dens fields (but some galaxies also get many matches, so leave this large!)
        self.maxPS1matches = run_config["MaxPS1matches"]

        # remove source in dens fields (almost certainly star)
        self.maxGaiamatches = run_config["MaxGaiamatches"]

        # Gaia parallax cut
        self.maxParallaxSNR = run_config["MaxParallaxSNR"]

        # Distance to check multiple matches in the Gaia (ie, double stars)
        self.singleGaiaDist = run_config["SingleGaiaDist"]

        # Difference between Gaia mag and PS1_converted_Gaia (large values imply more flux in PS1, ie extended sources)
        # only applied to source with sgscore1 >= minSgscore
        self.minDeltaG = run_config["MinDeltaG"]

        # stay away from the Galactic plane
        self.minGalLat = run_config["MinGalLat"]

        # use only most recent detection (attempt to simulate real time)
        self.LastOnly = run_config["LastOnly"]

        # set path for catsHTM
        # self.catsHTM_path = urlparse(base_config['catsHTM.default']).path
        # self.logger.info("using catsHTM files in %s"%self.catsHTM_path)

        # Instance dictionaries later used in method apply
        # remove upper limits
        isdetect_flt = {"attribute": "candid", "operator": "is not", "value": None}

        # check bad pixels
        nbad_flt = {"attribute": "nbad", "operator": "<=", "value": 0}

        # similar to Marshal mandatory filter (Full Width Half Max assuming a Gaussian core, from SExtractor)
        fwhm_flt = {"attribute": "fwhm", "operator": "<", "value": 3.5}

        # similar to Marshal mandatory filter (Ratio: aimage / bimage)
        elong_flt = {"attribute": "elong", "operator": "<", "value": 1.4}

        # similar to Marshal mandatory filter (Difference: magap - magpsf )
        magdiff_flt = {"attribute": "magdiff", "operator": "<", "value": 0.4}

        # is not ref-science
        isdiff_flt1 = {"attribute": "isdiffpos", "operator": "!=", "value": "0"}

        # is not ref-science (again!)
        isdiff_flt2 = {"attribute": "isdiffpos", "operator": "!=", "value": "f"}

        # has host galaxy detected (removes orphans)
        distnr_flt = {"attribute": "distnr", "operator": ">", "value": 0}

        self._default_filters = [
            isdetect_flt,
            distnr_flt,
            nbad_flt,
            fwhm_flt,
            elong_flt,
            magdiff_flt,
            isdiff_flt1,
            isdiff_flt2,
        ]

        # conversion from PS1 to SDSS (https://arxiv.org/abs/1203.0297)
        flt = ("g", "r", "i", "z")
        self.A0 = {flt[i]: x for i, x in enumerate((0.013, -0.001, -0.005))}
        self.A1 = {flt[i]: x for i, x in enumerate((0.145, 0.004, 0.011))}
        self.A2 = {flt[i]: x for i, x in enumerate((0.019, 0.007, 0.010))}

        self.logger.info("NuclearFilter initialized")

    # convert from PS1/ZTF to Gaia
    def alert_to_Gaia(self, alert, force_ztf=False):

        # get the PS1 PSF photometry
        ps1_mag = {k: alert.get_values("s" + k + "mag1")[0] for k in ("g", "r", "i")}

        # if we have a problem with the PS1 photometry, try to use ZTF
        if (min([ps1_mag[k] for k in ps1_mag]) < 0) or force_ztf:
            if not force_ztf:
                self.logger.warning(
                    "no decent PS1 photometry ({0:0.2f} {1:0.2f} {2:0.2f}), trying to use ZTF to estimate Gaia flux".format(
                        ps1_mag["g"], ps1_mag["r"], ps1_mag["i"]
                    )
                )

            magnr_g = alert.get_values(
                "magnr", filters={"attribute": "fid", "operator": "==", "value": 1}
            )
            magnr_R = alert.get_values(
                "magnr", filters={"attribute": "fid", "operator": "==", "value": 2}
            )

            if len(magnr_g) and len(magnr_R):

                ztf_g, ztf_r = np.median(magnr_g), np.median(magnr_R)
                gr = ztf_g - ztf_r

                Gaia_from_ztf = (
                    ztf_g - 0.0662 - 0.7854 * gr - 0.2859 * gr**2 + 0.0145 * gr**3
                )  # from Jordi (2010) #https://arxiv.org/abs/1008.0815
                return Gaia_from_ztf
            else:
                self.logger.warning(
                    "not enough ZTF detections (g={0}, R={1}) to compute Gaia mag".format(
                        len(magnr_g), len(magnr_R)
                    )
                )
                return None

        else:
            # convert from PS1 to SDSS (https://arxiv.org/abs/1203.0297)
            gr_ps1 = ps1_mag["g"] - ps1_mag["r"]
            sdss_mag = {
                flt: ps1_mag[flt]
                + self.A0[flt]
                + self.A1[flt] * (gr_ps1)
                + self.A2[flt] * (gr_ps1) ** 2
                for flt in ("g", "i")
            }

            gi = sdss_mag["g"] - sdss_mag["i"]

            Gaia_from_sdss = (
                sdss_mag["g"]
                - 0.098958
                - 0.6758 * gi
                - 0.043274 * gi**2
                + 0.0039908 * gi**3
            )  # from Jordi (2010) #https://arxiv.org/abs/1008.0815
            return Gaia_from_sdss

    def get_gaia(self, alert):
        """
        function to run catsHTM.cone_search and convert to np.array
        """
        # some files are corrupted, we have to catch the exception
        ra, dec = np.median(alert.get_values("ranr")), np.median(
            alert.get_values("decnr")
        )
        try:
            srcs, colnames, colunits = catsHTM.cone_search(
                "GAIADR2",
                np.radians(ra),
                np.radians(dec),
                30,
                catalogs_dir=self.catsHTM_path,
            )
            if len(srcs) == 0:
                self.logger.info('''no Gaia sources within 30"''')
                return None, None

        except OSError as ose:
            self.logger.warning("OSError from catsHTM %s" % str(ose))
            return None, None

        # convert output to numpy array
        gaia_arr = np.zeros(len(srcs[:, 0]), dtype=[(k, "f8") for k in colnames])
        for i, k in enumerate(gaia_arr.dtype.names):
            gaia_arr[k] = srcs[:, i]

        gaia_match_dist = 3600 * np.degrees(
            catsHTM.celestial.sphere_distance_fast(
                gaia_arr["RA"], gaia_arr["Dec"], np.radians(ra), np.radians(dec)
            )
        )
        return gaia_arr, gaia_match_dist

    def gaia_veto(self, alert):

        gaia_arr, gaia_match_dist = self.get_gaia(alert)

        if gaia_arr is None:
            return False

        if len(gaia_arr) > self.maxGaiamatches:
            self.why = (
                """number of Gaia matches within 30" is {0}, which is > {1}""".format(
                    len(gaia_arr), self.maxGaiamatches
                )
            )
            self.logger.info(self.why)
            return True

        idxBright = gaia_match_dist < self.brightObjDist
        if sum(idxBright):
            gaia_brighest_match = np.min(gaia_arr["Mag_G"][idxBright])
        else:
            self.logger.info('''no Gaia matches within 15"''')
            return False

        if gaia_brighest_match < self.brightGaiaMag:
            self.why = """within {0:0.1f}" of  bright star in Gaia: G={1:0.2f}, which is {2:0.1f}""".format(
                self.brightObjDist, gaia_brighest_match, self.brightGaiaMag
            )

        idxNear = gaia_match_dist < self.singleGaiaDist
        if sum(idxNear) > 1:
            self.why = """number of Gaia matches within {0:0.1f}" is {1}, which is one too many""".format(
                self.singleGaiaDist, sum(idxNear)
            )
            self.logger.info(self.why)
            return True
        elif sum(idxNear) == 0:
            self.why = '''no Gaia matches within {0:0.1f}"'''.format(
                self.singleGaiaDist
            )
            self.logger.info(self.why)
            return False

        # with one left, we can look at the nearest match only
        gaiam = gaia_arr[idxNear][0]

        parallax_snr = gaiam["Plx"] / gaiam["ErrPlx"]
        if parallax_snr > self.maxParallaxSNR:
            self.why = """parallax SNR={0:0.1f}, which is > {1:0.1f}""".format(
                parallax_snr, self.maxParallaxSNR
            )
            self.logger.info(self.why)
            return True

        # compute the difference between Gaia and ground-based photometry (PS1/ZTF)
        predicted_G = self.alert_to_Gaia(alert)
        if predicted_G is not None:
            delta_G = gaiam["Mag_G"] - predicted_G
            self.logger.info("Gaia - Ground = {0:0.2f}".format(delta_G))
        else:
            self.why = """Can't compute Gaia magnitude"""
            self.logger.info(self.why)
            return True

        # if the score is in the "safe" range, don't allow veto based on photometry
        sg1 = alert.get_values("sgscore1")[0]
        if sg1 < self.minSafeSgscore:
            self.logger.info(
                "not vetoing this alert because sgscore={0:0.2f}, which is < {1:0.2f}".format(
                    sg1, self.minSafeSgscore
                )
            )
            return False

        if delta_G < self.minDeltaG:
            self.why = "Gaia - Ground = {0:0.2f}, which is < {1:0.2f}".format(
                delta_G, self.minDeltaG
            )
            self.logger.info(self.why)
            return True

    # compute a weighted mean
    def wmean(self, x, w):
        return np.sum(x * w) / np.sum(w)

    def apply(self, alert):
        """
        Mandatory implementation.
        To exclude the alert, return *None*
        To accept it, either
                * return self.on_match_t2_units
                * return a custom combination of T2 unit names

        Make a selection on:
        - the distance between the transient and host in reference image
        - the Real/Bogus sore
        - the distance to a bright star
        """

        # first check we have an extended source (note this can remove flares from faint galaxies that missclassified in PS1)
        # these will have to be dealt with in the orphan/faint filter

        sgscore1 = alert.get_values("sgscore1")
        if len(sgscore1) == 1:

            alert_pps = alert.pps[0]
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

            # collect all PS1 magnitude into a dict
            psmag1 = {k: alert_pps["s" + k + "mag1"] for k in ("g", "r", "i", "z")}
            psmag2 = {k: alert_pps["s" + k + "mag2"] for k in ("g", "r", "i", "z")}
            psmag3 = {k: alert_pps["s" + k + "mag3"] for k in ("g", "r", "i", "z")}

        # check that we dont have very old (pre v1.8) schema
        else:

            self.why = "this schema doesnt have sgscore1, rejecting this alert"
            self.logger.warning(self.why)
            return None

        # also check the sgscore of the close match
        if (sgscore2 > self.minSgscore) and (abs(distpsnr2) < self.closePS1dist):
            self.why = (
                "sgscore2 = {0:0.2f}, which is > {1:0.2f}  (distpsnr2={2:0.2f})".format(
                    sgscore2, self.minSgscore, distpsnr2
                )
            )
            self.logger.info(self.why)
            return None

        # check problems with photometry (likely for saturated stars)
        if (psmag1["r"] < 0) and not (
            (psmag2["r"] > 0) and (abs(distpsnr2) < self.closePS1dist)
        ):
            self.why = "1st PS1 match r-band photometry is faulty: mag={0:0.2f} (dist={1:0.2f})".format(
                psmag1["r"], distpsnr1
            )
            if distpsnr2 < self.closePS1dist:
                self.why += "; 2nd PS1 match r-mag={0:0.2f} (dist={1:0.2f})".format(
                    psmag1["r"], distpsnr1
                )
            self.logger.info(self.why)
            return None

        # check that the nearest PS1 source is not too far
        if abs(distpsnr1) > self.minDistPS1source:
            self.why = (
                "distance to 1st PS1 match is {0:0.2f}, which is > {1:0.2f}".format(
                    distpsnr1, self.minDistPS1source
                )
            )
            self.logger.info(self.why)
            return None

        # check for nearby bright stars
        for k in ("g", "r", "i", "z"):

            if abs(psmag1[k]) < self.brightPS1Mag[k]:
                self.why = "1st PS1 match in {0}-band, mag={1:0.2f}, which is < {2:0.2f} (dist={3:0.2f})".format(
                    k, psmag1[k], self.brightPS1Mag[k], distpsnr1
                )
                self.logger.info(self.why)
                return None

            if psmag2 is not None:  # old schema check
                if (abs(psmag2[k]) < self.brightPS1Mag[k]) and (
                    abs(distpsnr2) < self.brightObjDist
                ):
                    self.why = "2nd PS1 match in {0}-band, mag={1:0.2f}, which is < {2:0.2f} (dist={3:0.2f})".format(
                        k, psmag2[k], self.brightPS1Mag[k], distpsnr2
                    )
                    self.logger.info(self.why)
                    return None

            if psmag3 is not None:
                if (abs(psmag3[k]) < self.brightPS1Mag[k]) and (
                    abs(distpsnr3) < self.brightObjDist
                ):
                    self.why = "3rd  PS1 match in {0}-band, mag={1:0.2f}, which is < {2:0.2f} (dist={3:0.2f})".format(
                        k, psmag3[k], self.brightPS1Mag[k], distpsnr3
                    )
                    self.logger.info(self.why)
                    return None

            # don't use the code below because it will remove sources next to objects
            # that were detected in just one pan-starrs band and thus have srmag=-999
            #
            # if ((psmag2['r']<0) or (psmag['g']<0)) and (abs(distpsnr2)< self.brightObjDist):
            # 	self.why = "2nd PS1 match saturated(?) sgmag={0:0.2f} srmag={1:0.2f} (dist={2:0.2f})".format(psmag2['g'], psmag2['r'], distpsnr2)
            # 	self.logger.info(self.why)
            # 	return None

        # get RealBogus scores for observations, check number of bad pixels
        tuptup = alert.get_ntuples(["rb", "jd", "magnr"], filters=self._default_filters)

        # check that we have anything
        if len(tuptup) == 0:
            self.why = "nothing passed default filter"
            self.logger.info(self.why)
            return None

        # now get the tuples
        rb_arr, jd_arr, magnr_arr = map(np.array, zip(*tuptup))

        # check number of detections in all bands
        if len(jd_arr) < self.minDetections:
            self.why = "only {0} detections pass default filter, but we require at least {1} (in each band)".format(
                len(jd_arr), self.minDetections
            )
            self.logger.info(self.why)
            return None

        # check that source is not too bright in ZTF ref img
        if self.brightRefMag > np.min(magnr_arr) > 0:
            self.why = "min(magnr)={0:0.2f}, which is < {1:0.1f}".format(
                np.min(magnr_arr), self.brightRefMag
            )
            self.logger.info(self.why)
            return None

        # check source density --> todo replace this check on Gaia source density, or allow veto for big galaxies!
        nmtchps = alert.get_values("nmtchps")[0]
        if nmtchps > self.maxPS1matches:
            self.why = "nmtchps={0}, which is > {1}".format(nmtchps, self.maxPS1matches)
            self.logger.info(self.why)
            return None

        # check Galactic Latitude
        ra, dec = alert.get_values("ranr")[0], alert.get_values("decnr")[0]
        gc = astropy.coordinates.SkyCoord(
            ra=ra * astropy.units.deg, dec=dec * astropy.units.deg
        ).galactic
        gal_l, gal_b = gc.galactic.l.deg, gc.galactic.b.deg
        if abs(gal_b) < self.minGalLat:
            self.why = (
                "Too close to the plane; l,b=({0:0.5f} {1:0.5f}), |b|<{2:0.1f}".format(
                    gal_l, gal_b, self.minGalLat
                )
            )
            self.logger.info(self.why)
            return None

        # if we want, only check last observation
        if self.LastOnly:

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
                return None

            rb_arr = [
                rb_arr[np.argmax(jd_arr)]
            ]  # make sure rb check below is only for last detection

        # if no detections pass real bogus, remove
        if max(rb_arr) < self.minRealBogusScore:
            self.why = "max(rb)={0:0.2f}, which is  < {1:0.2f}".format(
                max(rb_arr), self.minRealBogusScore
            )
            self.logger.info(self.why)
            return None

        # do cut on moving sources (with all good detections)
        dt = abs(np.sort(jd_arr)[1:] - np.sort(jd_arr)[0:-1])

        # first check that we dont have bunch of duplicates
        if not sum(dt > 0):
            self.why = "number of detections={0}, but time difference between all is zero".format(
                len(dt)
            )
            self.logger.info(self.why)
            return None

        dt = dt[dt > 0]
        if np.max(dt) < self.minDeltaJD:
            self.why = "potential mover, number of good detections={0}; max(time diff)={1:1.3f} h, which is < {2:0.3f} h".format(
                len(jd_arr), max(dt) * 24, self.minDeltaJD * 24
            )
            self.logger.info(self.why)
            return None

        # Check time between detecetions (this could be made more fancy using upper limits etc.)
        if np.min(dt) > self.maxDeltaJD:
            self.why = "number of good detections={0}; min(time diff)={1:1.1f} h, which is > {2:0.1f} h".format(
                len(jd_arr), min(dt) * 24, self.maxDeltaJD * 24
            )
            self.logger.info(self.why)
            return None

        # get some more arrays
        distnr_arr, magpsf_arr, sigmapsf_arr, rb_arr, fwhm_arr, fid_arr = map(
            np.array,
            zip(
                *alert.get_ntuples(
                    ["distnr", "magpsf", "sigmapsf", "rb", "fwhm", "fid"],
                    filters=self._default_filters,
                )
            ),
        )

        # check sufficient detections in each band
        idx_all = np.repeat(True, len(distnr_arr))
        idx_g = fid_arr == 1
        idx_r = fid_arr == 2

        if (sum(idx_g) < self.minDetections) and (sum(idx_r) < self.minDetections):
            self.why = "number of good detections in (g, r)=({0}, {1}), which is < {2:0}".format(
                sum(idx_g), sum(idx_r), self.minDetections
            )
            self.logger.info(self.why)
            return None

        # apply magnitude cut to brightest detection
        if np.min(magpsf_arr) > self.diffmagLimit:
            self.why = "too faint, min(mag)={0:0.2f}, which is > {1:0.1f}".format(
                np.min(magpsf_arr), self.diffmagLimit
            )
            self.logger.info(self.why)
            return None

        # apply lower limit on flux increase
        if np.min(magpsf_arr - magnr_arr) > self.maxDeltaMag:
            self.why = "flux increase too small, min(magpsf-magnr)={0:0.2f}, which is > {1:0.2f}".format(
                np.min(magpsf_arr - magnr_arr), self.maxDeltaMag
            )
            self.logger.info(self.why)
            return None

        # allow Gaia to veto
        G_veto = self.gaia_veto(alert)
        if G_veto:
            return None

        # if we make it this far, compute the host-flare distance, using only (decent-enough) detections
        ra_arr, dec_arr, ranr_arr, decnr_arr = map(
            np.array,
            zip(
                *alert.get_ntuples(
                    ["ra", "dec", "ranr", "decnr"], filters=self._default_filters
                )
            ),
        )

        offset_sigma = np.clip(
            0.27 + 1.66 * (sigmapsf_arr - 0.1), 0.1, 10
        )  # equation from the NedStark paper
        snr_weight = 1 / offset_sigma**2

        # compute a few different measures of the distance
        # we also compute these for each band seperately
        rb_arr = np.clip(rb_arr, 0.01, 1)  # remove zero scores
        my_weight = 1 / (
            np.clip(fwhm_arr, 1.5, 10) * np.clip(sigmapsf_arr, 0.01, 1)
        )  # combine differen measures for how good the distnr measurement is

        flagagn = 0

        for filterid, filtername in zip([1, 2], ["g", "r"]):
            try:
                time, flux = map(
                    np.array,
                    zip(
                        *alert.get_tuples(
                            "jd",
                            "magpsf",
                            filters=[
                                {
                                    "attribute": "candid",
                                    "operator": "is not",
                                    "value": None,
                                },
                                {
                                    "attribute": "fid",
                                    "operator": "==",
                                    "value": filterid,
                                },
                            ],
                        )
                    ),
                )
                fluxerr = np.asarray(
                    alert.get_values(
                        "sigmapsf",
                        filters=[
                            {
                                "attribute": "candid",
                                "operator": "is not",
                                "value": None,
                            },
                            {"attribute": "fid", "operator": "==", "value": filterid},
                        ],
                    )
                )

                sortind = np.argsort(time)
                time = time[sortind]
                flux = flux[sortind]
                fluxerr = fluxerr[sortind]
                results = qso_fit(
                    time,
                    flux,
                    fluxerr,
                    filter=filtername,
                    sys_err=0.0,
                    return_model=True,
                )

                if results["class"] == "qso":
                    flagagn = 1
            except ValueError:
                pass

        if flagagn == 0:
            self.why = "Ambiguous or non-qso type according to light curve modelling, or not enough points in light curve"
            self.logger.info(self.why)
            return None

        for idx, bnd in zip([idx_g, idx_r], ["g", "r"]):

            if sum(idx) >= self.minDetections:

                mean_dist = 3600 * np.sqrt(
                    np.mean(ra_arr - ranr_arr) ** 2
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

                if mean_dist < self.maxDeltaRad and mean_dist > self.minDeltaRad:
                    self.why = "pass on mean dist={0:0.2f}, band={1}; detections used={2}".format(
                        mean_dist, bnd, sum(idx)
                    )
                    self.logger.info(self.why)
                    return self.on_match_t2_units

                if median_dist < self.maxDeltaRad and median_dist > self.minDeltaRad:
                    self.why = "pass on median dist={0:0.2f}, band={1}; detections used={2}".format(
                        median_dist, bnd, sum(idx)
                    )
                    self.logger.info(self.why)
                    return self.on_match_t2_units

                if (
                    weighted_dist < self.maxDeltaRad
                    and weighted_dist > self.minDeltaRad
                ):
                    self.why = "pass on weighted dist={0:0.2f}, band={1}; detections used={2}".format(
                        weighted_dist, bnd, sum(idx)
                    )
                    self.logger.info(self.why)
                    return self.on_match_t2_units

        # if none of the measures of the host-flare pass the cut, reject this alert
        try:
            self.why = "mean/median/weighted dist in {0}-band = ({1:0.2f}/{2:0.2f}/{3:0.2f}), which is > {4:0.2f}".format(
                bnd, mean_dist, median_dist, weighted_dist, self.maxDeltaRad
            )
            self.logger.info(self.why)
            return None
        except NameError:
            self.why = "offset distance not created, likely not enough data passed; number left in r&g band ={0}".format(
                sum(idx)
            )
            self.logger.info(self.why)
            return None
