#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t0/GaiaVeto.py
# License           : BSD-3-Clause
# Author            : sjoertvv <sjoert@astro.umd.edu>
# Date              : 26.02.2018
# Last Modified Date: 21.10.2022
# Last Modified By  : simeon.reusch@desy.de

from typing import Optional, Sequence, Union, Any

from ampel.ztf.base.CatalogMatchUnit import CatalogMatchUnit
import pandas as pd  # type: ignore
import numpy as np

from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol


class GaiaVetoMixin(CatalogMatchUnit):

    maxGaiamatches: int = 30  #: remove source in dens fields (almost certainly star)
    brightObjDist: float = (
        15.0  #: max distance in arcsec for checking nearby bright stars in Gaia or PS1
    )
    brightGaiaMag: float = 13.5  #: Gaia bright star removal in mulitple bands
    maxParallaxSNR: float = 3.0  #: Gaia parallax cut

    singleGaiaDist: float = (
        1.0  #: Distance to check multiple matches in Gaia (ie, reject double stars)
    )

    minDeltaG: float = 0.2  #: Difference between Gaia mag and PS1_converted_Gaia (large values imply more flux in PS1, ie extended sources) only applied to source with sgscore1 >= minSgscore

    minSafeSgscore: float = (
        0.5  # min star-galaxy score for which we dont allow a Gaia-based veto
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # conversion from PS1 to SDSS (https://arxiv.org/abs/1203.0297)
        flt = ("g", "r", "i", "z")
        self.__A0 = {flt[i]: x for i, x in enumerate((0.013, -0.001, -0.005))}
        self.__A1 = {flt[i]: x for i, x in enumerate((0.145, 0.004, 0.011))}
        self.__A2 = {flt[i]: x for i, x in enumerate((0.019, 0.007, 0.010))}

    def alert_to_Gaia(
        self, alert: AmpelAlertProtocol, force_ztf=False
    ) -> Optional[float]:
        """
        convert from PS1/ZTF to Gaia (used in Gaia veto function)
        """

        # get the PS1 PSF photometry
        ps1_mag = {k: alert.get_values("s" + k + "mag1")[0] for k in ("g", "r", "i")}

        # if we have a problem with the PS1 photometry, try to use ZTF
        if (min([ps1_mag[k] for k in ps1_mag]) < 0) or force_ztf:
            if not force_ztf:
                self.logger.warn(
                    "no decent PS1 photometry ({0:0.2f} {1:0.2f} {2:0.2f}), trying to use ZTF to estimate Gaia flux".format(
                        ps1_mag["g"], ps1_mag["r"], ps1_mag["i"]
                    )
                )

            magnr_g = alert.get_values(
                "magnr", filters=[{"attribute": "fid", "operator": "==", "value": 1}]
            )
            magnr_R = alert.get_values(
                "magnr", filters=[{"attribute": "fid", "operator": "==", "value": 2}]
            )

            if len(magnr_g) and len(magnr_R):

                ztf_g, ztf_r = np.median(magnr_g), np.median(magnr_R)
                gr = ztf_g - ztf_r

                # Gaia_from_ztf = ztf_g -0.0662 -0.7854*gr -0.2859*gr**2 +0.0145*gr**3          # from Jordi (2010) #https://arxiv.org/abs/1008.0815
                Gaia_from_ztf = (
                    ztf_g
                    - 0.038025
                    - 0.76988 * gr
                    - 0.1931 * gr**2
                    + 0.0060376 * gr**3
                )  # from Evan (2018) #https://arxiv.org/abs/1804.09368
                return Gaia_from_ztf
            else:
                self.logger.warn(
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
                + self.__A0[flt]
                + self.__A1[flt] * (gr_ps1)
                + self.__A2[flt] * (gr_ps1) ** 2
                for flt in ("g", "i")
            }

            # convert from SDSS to Gaia
            gi = sdss_mag["g"] - sdss_mag["i"]
            # Gaia_from_sdss = sdss_mag['g']-0.0912   -0.5310*gi  -0.1042*gi**2   +0.0068*gi**3     # (DR0, Jordi+18, https://arxiv.org/abs/1008.0815)
            Gaia_from_sdss = (
                sdss_mag["g"]
                - 0.074189
                - 0.51409 * gi
                - 0.080607 * gi**2
                + 0.0016001 * gi**3
            )  # (DR2, Evans+18, https://arxiv.org/abs/1804.09368, v3)
            return Gaia_from_sdss

    # no useful type hints for astropy.Table
    def get_Gaia(self, alert: AmpelAlertProtocol) -> Optional[pd.DataFrame]:
        """
        function to run catsHTM.cone_search and convert to np.array
        """
        # some files are corrupted, we have to catch the exception
        ra, dec = np.median(alert.get_values("ranr")), np.median(
            alert.get_values("decnr")
        )

        matches = self.cone_search_all(
            ra,
            dec,
            [
                {
                    "name": "GAIADR2",
                    "use": "catsHTM",
                    "rs_arcsec": 30,
                    "keys_to_append": ["Mag_G", "Plx", "ErrPlx"],
                }
            ],
        )[0]

        if not matches:
            return None

        return pd.DataFrame(
            [{"dist_arcsec": m["dist_arcsec"], **m["body"]} for m in matches]
        )

    def gaia_veto(self, alert: AmpelAlertProtocol) -> bool:
        """
        veto IPAC star/galaxy score based on Gaia:
         - parallax
         - too many close matches
         - too many matches
         - compactness: Gaia flux compared to PS1 (or ZTF) flux
        """

        if (gaia_df := self.get_Gaia(alert)) is None:
            return False

        if len(gaia_df) > self.maxGaiamatches:
            self.why = (
                """number of Gaia matches within 30" is {0}, which is > {1}""".format(
                    len(gaia_df), self.maxGaiamatches
                )
            )
            self.logger.info(
                '''too many Gaia matches within 30"''',
                extra={"n_Gaia_matches_within_30ac": len(gaia_df)},
            )
            return True

        idxBright = gaia_df["dist_arcsec"] < self.brightObjDist
        if sum(idxBright):
            gaia_brighest_match = np.min(gaia_df["Mag_G"][idxBright])
        else:
            self.logger.info(
                """no Gaia matches within {0:0.1f}" (not vetoing)""".format(
                    self.brightObjDist
                )
            )
            return False

        if gaia_brighest_match < self.brightGaiaMag:
            self.why = """within {0:0.1f}" of bright star in Gaia: G={1:0.2f}, which is {2:0.1f}""".format(
                self.brightObjDist, gaia_brighest_match, self.brightGaiaMag
            )
            self.logger.info(
                "close to Gaia bright star",
                extra={"gaia_brighest_match": gaia_brighest_match},
            )
            return True

        idxNear = gaia_df["dist_arcsec"] < self.singleGaiaDist
        if sum(idxNear) > 1:
            self.why = """number of Gaia matches within {0:0.1f}" is {1}, which is one too many""".format(
                self.singleGaiaDist, sum(idxNear)
            )
            self.logger.info(
                "Close Gaia match", extra={"n_close_Gaia_matches": int(sum(idxNear))}
            )
            return True
        elif sum(idxNear) == 0:
            self.logger.info("no Gaia match (not vetoing)")
            return False

        # with one left, we can look at the nearest match only
        gaiam = gaia_df[idxNear].iloc[0]

        if (
            gaiam["Plx"] is not None
            and gaiam["ErrPlx"] is not None
            and (parallax_snr := gaiam["Plx"] / gaiam["ErrPlx"]) > self.maxParallaxSNR
        ):
            self.why = """parallax SNR={0:0.1f}, which is > {1:0.1f}""".format(
                parallax_snr, self.maxParallaxSNR
            )
            self.logger.info(
                "gaia parallax, vetoing", extra={"parallax_snr": parallax_snr}
            )
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

        # if not in safe range, compute the difference between Gaia and ground-based photometry (PS1/ZTF)
        predicted_G = self.alert_to_Gaia(alert)
        if predicted_G is not None:
            delta_G = gaiam["Mag_G"] - predicted_G
            self.why = "Gaia - Ground = {0:0.2f}".format(delta_G)
            self.logger.info(None, extra={"delta_G": delta_G})
        else:
            self.why = """Can't compute Gaia magnitude from PS1 or ZTF data, vetoing"""
            self.logger.info(self.why)
            return True

        if delta_G < self.minDeltaG:
            self.why = """Gaia - Ground = {0:0.2f}, which is < {1:0.2f}""".format(
                delta_G, self.minDeltaG
            )
            self.logger.info("too compact in Gaia, vetoing", extra={"delta_G": delta_G})
            return True

        return False
