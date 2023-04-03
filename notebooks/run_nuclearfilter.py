import argparse
import collections
import json
import os
import pickle
import time
from pathlib import Path

import numpy as np
import requests
from ampel.log.AmpelLogger import AmpelLogger
from ampel.nuclear.t0.NuclearFilter import NuclearFilter
from ampel.ztf.alert.ZiAlertSupplier import ZiAlertSupplier
from ampel.ztf.t0.load.ZTFArchiveAlertLoader import ZTFArchiveAlertLoader
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper
from astropy.time import Time
from tqdm import tqdm

LOGGER = AmpelLogger.get_logger()

API_TOKEN = os.environ["AMPEL_API_ARCHIVE_TOKEN_PASSWORD"]

LOGGER.level = 1 << 12


class APIError(Exception):
    pass


def api_stream(token: str, date_start: str, date_end: str, token_only=False) -> str:
    """
    Initiate a stream query
    """
    endpoint = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3//streams/from_query"

    t_min_jd = Time(date_start).jd
    t_max_jd = Time(date_end).jd

    # t_min_jd = 2459170.663773
    # t_max_jd = 2459170.6637732

    query = {
        "jd": {
            "$gt": t_min_jd,
            "$lt": t_max_jd,
        },
        "candidate": {
            "rb": {"$gt": 0.3},
            "magpsf": {"$lte": 20},
            "ndethist": {"$gt": 7},
            "distpsnr1": {"$lte": 1},
            "isdiffpos": {"$in": ["t", "1"]},
            "nmtchps": {"$lte": 100},
        },
    }

    header = {"Authorization": "bearer " + token}

    response = requests.post(endpoint, json=query, headers=header)

    if not response.ok:
        print(f"Accessing stream not successful. Response: {response.json()}")
        raise APIError(f"{response.json()['detail'][0]['msg']}")

    else:
        resume_token = response.json()["resume_token"]

    print("Stream initiated.")
    if token_only:
        print(resume_token)
        quit()

    return resume_token


def run_filter(
    resume_token: str, date_start: str, date_end: str, sleepytime: int = 20
) -> list:
    """
    Get all alerts from stream (with the token) and
    run them through the nuclear filter
    """
    print("Attempting API request")
    stream_config = {
        "archive": "https://ampel.zeuthen.desy.de/api/ztf/archive/v3",
        "stream": resume_token,
    }

    basedir = Path("data")
    # if not basedir.is_dir():
    basedir.mkdir(parents=True, exist_ok=True)

    infile = basedir / f"{date_start}_to_{date_end}.pickle"

    # infile = Path("BS")

    if infile.is_file():
        print("Reading local file")
        with open(infile, "rb") as file:
            all_alerts_shaped = pickle.load(file)

    else:
        print(
            f"Now we get so sleep for {sleepytime} seconds to let the Archive ramp up."
        )
        time.sleep(sleepytime)

        alertloader = ZTFArchiveAlertLoader(**stream_config)

        all_alerts = []
        all_ztfids = []

        for i, alert in enumerate(alertloader.get_alerts()):
            print(i)
            all_ztfids.append(alert["objectId"])
            all_alerts.append(alert)

        if all_ztfids:
            save_list(
                names=all_ztfids,
                title="first_stage",
                date_start=date_start,
                date_end=date_end,
            )

        all_alerts_shaped = [
            ZiAlertSupplier.shape_alert_dict(a, ["FilterTest"]) for a in all_alerts
        ]

        outfile = basedir / f"{date_start}_to_{date_end}.pickle"
        with open(outfile, "wb") as file:
            pickle.dump(all_alerts_shaped, file)

    # # Only defaults

    filter_config = {}
    t0filter = NuclearFilter(
        **filter_config,
        resource={
            "ampel-ztf/catalogmatch": "https://ampel.zeuthen.desy.de/api/catalogmatch/"
        },
        logger=LOGGER,
    )
    t0filter.post_init()
    return_vals = []

    ztf_ids = []

    for alert in tqdm(all_alerts_shaped):
        returncode = t0filter.process(alert)
        if returncode == True:
            ztfid = ZTFIdMapper.to_ext_id(alert.stock)
            ztf_ids.append(ztfid)

    ztf_ids = list(set(ztf_ids))
    print(len(ztf_ids))

    print(
        f"Processed {len(all_alerts_shaped)} alerts, found {len(ztf_ids)} candidates."
    )

    return ztf_ids


def initiate_token(
    initiate: bool, date_start: str, date_end: str, token_only: bool
) -> str:
    """
    Get a stream query token from the Ampel API and save token to disk
    If not passed --initiate, load token from disk
    """
    if initiate:
        resume_token = api_stream(
            API_TOKEN, date_start=date_start, date_end=date_end, token_only=token_only
        )
        with open("resume_token.json", "w") as outfile:
            json.dump({"resume_token": resume_token}, outfile)
        print(f"Generated resume token and wrote to file: {resume_token}")
        return resume_token
    else:
        infile = Path("resume_token.json")
        if infile.is_file():
            with open("resume_token.json", "r") as infile:
                resume_token_dict = json.load(infile)
                resume_token = resume_token_dict["resume_token"]
                print(f"Using saved resume token: {resume_token}")
                return resume_token
        else:
            raise ValueError("No token found, run with --initiate first")
        return None


def save_list(names, title, date_start, date_end):
    """
    Save list of ztf_ids to textfile
    """
    to_save = np.asarray(list(set(names)))
    outfile = Path("data") / f"{date_start}_to_{date_end}_{title}.txt"
    np.savetxt(
        outfile,
        to_save,
        delimiter=",",
        fmt="%s",
    )


if __name__ == "__main__":
    """ """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--initiate",
        action="store_true",
        help="Initiate a new stream (and resume token)",
    )
    parser.add_argument(
        "--token",
        action="store_true",
        help="Only print and save a new token, do nothing else",
    )
    cli_args = parser.parse_args()
    initiate = cli_args.initiate
    token_only = cli_args.token

    date_start = "2022-03-01"
    date_end = "2022-03-31"

    resume_token = initiate_token(
        initiate, date_start=date_start, date_end=date_end, token_only=token_only
    )

    passed = run_filter(
        resume_token=resume_token,
        date_start=date_start,
        date_end=date_end,
        sleepytime=1,
    )

    if len(passed) > 0:
        save_list(
            names=passed, title="passed", date_start=date_start, date_end=date_end
        )
