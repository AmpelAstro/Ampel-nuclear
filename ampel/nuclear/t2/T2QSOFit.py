#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : /ampel/nuclear/t2/T2QSOFit.py
# License           : BSD-3-Clause
# Author            : matteo.giomi@desy.de
# Date              : 12.11.2018
# Last Modified Date: 10.26.2020
# Last Modified By  : svv (sjoert@nyu.edu)


from typing import Dict, List, Any
from astropy.table import Table  # type: ignore

from ampel.log import VERBOSE
from ampel.view.LightCurve import LightCurve
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.contrib.ztfbh.charlotte_.qso_fit import qso_fit  # type: ignore


class T2QSOFit(AbsLightCurveT2Unit):
    """
    Fit the light curve of this transient with a damped random model. The goal
    is to use this to help discriminating QSO like variability.
    These functions are copied from the contributions in Amepl-contrib-ZTFBh.
    """

    filter_names: Dict[int, str] = {1: "g", 2: "r", 3: "i"}
    filter_ids: List[int] = [1, 2, 3]
    min_ndet: int = 3

    def run(self, light_curve: LightCurve):
        """
        Parameters
        -----------
        light_curve: `ampel.view.LightCurve` instance.
                See the LightCurve docstring for more info.

        run_config: `dict` or None
                configuration parameter for this job. If none is given, the
                default behaviour would be to compute the metrics for the light
                curve in all the three bands (if the corresponding light curves have
                some detection), to use zero-order (step-like) interpoaltion
                between the LC points, and to exclude points with negative detections
                (having isdiffpos in ['f', 0]).

                These defaults can be changed by the following keys of the
                run_config dictionary:

                        filter_ids: `list` or `tuple`
                                list of ints in the range 1 to 3 that specify the filter
                                ids for which this job has to run. 1=g, 2=r, and 3=i

                        min_ndet: `int`
                                minimum number of points in the light curve (upper limits are
                                not included) needed to perform the fit. Default to self.default_min_ndet

        :returns: the returned dictionary containing, for each of the ZTF band in the LC, the following info:

        chi2/nu - classical variability measure
        chi2_qso/nu - for goodness of fit given fixed parameters
        chi2_qso/nu_extra - for parameter fitting, add to chi2/nu
        chi^2/nu_NULL - expected chi2/nu for non-qso variable

        signif_qso - significance chi^2/nu<chi^2/nu_NULL (rule out false alarm)
        signif_not_qso - significance chi^2/nu>1 (rule out qso)
        signif_vary - significance that source is variable
        class - resulting source type (ambiguous, not_qso, qso)

        model - time series prediction for each datum given all others (iff return_model==True,
        disabled for now)
        dmodel - model uncertainty, including uncertainty in data

        If, for a given band, the number of points in the light curve is below threshold,
        then the corresponding entry in the dict would be None.

        E.g.:
        {
                'g':
                        {'lvar': -4.1892290115356445,
                        'ltau': 2.9646145057678224,
                        'chi2/nu': 7.734741788611833,
                        'nu': 33.0,
                        'chi2_qso/nu': 7.718420939516062,
                        'chi2_qso/nu_NULL': 13.629683453346471,
                        'signif_qso': 1.9511898541054034,
                        'signif_not_qso': 5.618242539463743,
                        'signif_vary': 12.592220943514539,
                        'class': 'not_qso'
                        },
                'r':
                        {
                        'lvar': -4.401520004272461,
                        'ltau': 3.1661400032043456,
                        'chi2/nu': 77.98072995085302,
                        'nu': 36.0,
                        'chi2_qso/nu': 58.64461068420829,
                        'chi2_qso/nu_NULL': 74.84431170466726,
                        'signif_qso': 1.1980841481386015,
                        'signif_not_qso': 10.066851343992596,
                        'signif_vary': 51.23134323248932,
                        'class': 'not_qso'
                        },
                'i': None
        }
        """

        # run on the single bands individually
        out_dict: Dict[str, Any] = {}
        # debug = self.logger.level > VERBOSE
        for fid in self.filter_ids:

            filters = [{"attribute": "fid", "operator": "==", "value": fid}]

            # if debug:
            self.logger.debug(
                f"Fitting QSO variability to light curve for filter id {fid}"
                f" ({self.filter_names[fid]}-band)"
            )
            self.logger.debug("Applying filter: %s" % repr(filters))

            # get detections and errors time series
            pps = light_curve.get_ntuples(("jd", "magpsf", "sigmapsf"), filters=filters)

            if not pps or len(pps) < self.min_ndet:

                self.logger.debug(
                    "Lightcurve has only %d points, you asked for %d at least"
                    % (len(pps) if pps else 0, self.min_ndet)
                )
                out_dict[self.filter_names[fid]] = None
                continue

            # sort according to time
            photo = Table(rows=pps, names=("jd", "mag", "dmag"))
            photo.sort("jd")

            # run damped random walk fit & save result
            results = qso_fit(
                photo["jd"],
                photo["mag"],
                photo["dmag"],
                filter=self.filter_names[fid],
                sys_err=0.0,
                return_model=False,
            )
            out_dict[self.filter_names[fid]] = results

        self.logger.debug("output", extra=out_dict)

        # return the info as dictionary
        return out_dict
