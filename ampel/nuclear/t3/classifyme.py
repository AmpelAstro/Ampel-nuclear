#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t3/classifyme.py
# License           : BSD-3-Clause
# Author            : mitchell@nyu.edu & sjoert@nyu.edu
# Date              : 08.06.2020
# Last Modified Date: 11.12.2020
# Last Modified By  : sjoert@nyu.edu

"""
container with a function to get classification information 
this is needed for mulitple T3s of ZTFbh
"""

import json
import io


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
    elif hasattr(obj, "isoformat"):
        return obj.isoformat()
    else:
        return obj


def get_class(tran_view, dropbox, logger, write=False):
    """
    Collect classifications from T2 transient view info,
    and our own neoWISE measurment.

    Input is:
     - trans_view: the t2 output which contains the catalogs matches and Fritz/TNS/Growth info
     - dropbox: instance to read/write (from dropboxIO, inititated by each T2)
     - logger: to keep chatting

    The following order applies to accepting classification:
     Fritz
     Growth
     TNS
     static catalog matches from T2
     neoWISE variability measurments form our T3

    This function also returns some extra information about the neoWISE measurements.

    >> class_str, extra_info = get_class(tran_view dropbox, logger, write=True)

    """

    growth_report = tran_view.extra.get("GROWTHMarshalReport")
    tns_report = reports[0] if (reports := tran_view.extra.get("TNSReports")) else None
    fritz_report = tran_view.extra.get("FritzReport")

    alert_base = (
        "/mampel/alerts/" + "20" + tran_view.stock["name"][0][3:5]
    )  # lets assume ZTF ends before 2099 :)

    # read  WISE light curve and color information
    wise_path = alert_base + "/{0}/{0}_neoWISE.json".format(tran_view.stock["name"][0])
    try:
        infile = dropbox.read_file(wise_path)
        if isinstance(infile, str):
            wise_json = json.loads(infile)
        else:
            wise_json = infile.json()
        logger.info("loading wise_class for {}".format(tran_view.stock["name"][0]))
        wise_class = wise_json["wise_class"]
        wise_info = wise_json["info_list"][
            0
        ]  # this now this just contains the neoWISE light curve info
    except FileNotFoundError:
        wise_class = None
        wise_info = ""

    # read the t2 records, take care of rouge dicts that mess with json
    # ampel_matches = tran_view.get_t2_result("T2CatalogMatch") or {}
    ampel_matches = tran_view.get_latest_t2_body(unit="T2CatalogMatch")

    xmatches = {}
    for key in ampel_matches:
        try:
            xmatches[key] = json.loads(json.dumps(ampel_matches[key]))
        # only report distance if rest of content cant be jsonified
        except TypeError:
            xmatches[key] = {"dist2transient": ampel_matches[key]["dist2transient"]}

    # write classifications
    if write:
        class_match = {
            "GROWTH": growth_report,
            "TNS": tns_report,
            "Fritz": fritz_report,
            "xmatch": xmatches,
        }

        buf = io.BytesIO()
        buf.write(bytes(json.dumps(jsonify(class_match), indent=3), "utf-8"))
        buf.seek(0)
        dropbox.put(
            alert_base + "/" + tran_view.stock["name"][0] + "/class_match.json",
            buf.read(),
        )

    # ---
    # Classify

    # if we have a Fritz classification, it is selected first ...
    if fritz_report:
        if "classifications" in fritz_report:

            # Frits can have mulitple classifications,
            # for now, just pick the first one
            if classification := next(
                iter(fritz_report["classifications"] or []), None
            ):

                # shorten some labels
                fritz_class = classification["classification"]
                fritz_class = fritz_class.replace("Tidal Disruption Event", "TDE")

                # deal with lower probability classifications
                if (prob := classification["probability"]) is not None and prob < 0.9:
                    fritz_class += "?"

                    return fritz_class, wise_info

    # ... otherwise Marshal ...
    if growth_report:
        if growth_report["classification"] is not None:
            return growth_report["classification"], wise_info

    # ... then TNS ...
    if tns_report:
        if tns_report["object_type"]["name"] is not None:
            return tns_report["object_type"]["name"], wise_info

    # ... or finally, do classification based on catalog matches
    classification = "None"
    if ampel_matches["milliquas"]:
        classification = "AGN"

    if ampel_matches["wise_color"]:
        classification = "AGN"

    if ampel_matches["varstars"]:
        if classification == "None":
            classification = "varstar"
        elif classification == "AGN":
            classification = "QSO"

    # ... or finally finally, use my neoWISE class
    if wise_class is not None:
        if classification == "None":
            classification = wise_class

    return classification, wise_info
