#!/usr/bin/env python
import argparse
import datetime
import os
import time
from typing import Any, Optional

import matplotlib  # type: ignore
import requests
from ampel.cli.config import get_user_data_config_path
from ampel.dev.DevAmpelContext import DevAmpelContext
from ampel.log.AmpelLogger import AmpelLogger
from ampel.nuclear.t0.NuclearFilter import NuclearFilter
from ampel.secret.AmpelVault import AmpelVault
from ampel.secret.DictSecretProvider import DictSecretProvider
from astropy.time import Time  # type: ignore

matplotlib.use("Agg")

logger = AmpelLogger.get_logger()

purge_db = False

temp_dir_base = os.path.join(os.getcwd(), "temp")


def run(
    initiate: bool = False,
    date: Optional[str] = None,
    daysago: Optional[int] = None,
    push_to_dropbox: bool = False,
):
    # Decent filter parameters
    filter_config = {
        "minDetections": 3,
        "maxDistPS1source": 0.5,
        "closePS1dist": 0.5,
        "maxDeltaRad": 0.5,
        "diffmagLimit": 20,
        "maxDeltaMag": 2.5,
        "resource": {
            "ampel-ztf/catalogmatch": "https://ampel.zeuthen.desy.de/api/catalogmatch/"
        },
    }

    t0filter = NuclearFilter(**filter_config, logger=logger)
    t0filter.post_init()

    if date:
        date_format = "%Y-%m-%d"
        req_date = str(datetime.datetime.strptime(date, date_format))
        startdate_jd = Time(req_date, format="iso", scale="utc").jd
        delta_t = 1
        enddate_jd = startdate_jd + delta_t

    elif daysago:
        enddate_jd = Time.now().jd
        startdate_jd = enddate_jd - float(daysago)

    else:
        startdate_jd = 2459899.04167
        enddate_jd = 2459899.045167

    query = {
        "jd": {
            "$gt": startdate_jd,
            "$lt": enddate_jd,
        },
        "candidate": {
            "distpsnr1": {"$lt": filter_config["maxDistPS1source"]},
            "rb": {"$gt": 0.3},
            "magpsf": {"$lte": 20},
            "sgscore1": {"$lte": 0.9},
            "ndethist": {
                "$gte": filter_config["minDetections"],
            },
            "isdiffpos": {"$in": ["t", "1"]},
            "nmtchps": {"$lte": 100},
        },
    }

    if initiate:
        endpoint = (
            "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/streams/from_query?"
        )
        header = {"Authorization": "bearer " + os.environ["AMPEL_ARCHIVE_TOKEN"]}

        response = requests.post(endpoint, headers=header, json=query)

        stream_token = response.json()["resume_token"]

        print("Sleeping for 10 seconds to give the archive DB some time")
        print(f"Your stream_token is: {stream_token}")
        with open("stream_token.txt", "w") as f:
            f.write(str(stream_token))

        time.sleep(10)

    else:
        stream_token_infile = open("stream_token.txt", "r")
        stream_token = stream_token_infile.read()
        print(f"Your stream_token is: {stream_token}")

    # Create a secret vault
    secrets = {
        "ztf/archive/token": os.environ["AMPEL_ARCHIVE_TOKEN"],
        "dropbox/token": os.environ["DROPBOX_TOKEN"],
        "fritz/token": os.environ["FRITZ_TOKEN"],
    }

    vault = AmpelVault([DictSecretProvider(secrets)])

    channel = "AMPEL_NUCLEAR_TEST"

    config = get_user_data_config_path()

    if not os.path.isfile(config):
        raise ImportError("No AMPEL config found. Run 'ampel config install' first.")

    ctx = DevAmpelContext.load(
        config=config,
        db_prefix="ampel-nuclear-test",
        purge_db=purge_db,
        one_db=True,
        vault=vault,
    )
    ctx.add_channel(name=channel, access=["ZTF", "ZTF_PUB", "ZTF_PRIV"])

    # Will use NED for spectroscopic redshifts.
    cat_conf = {
        "catalogs": {
            "SDSS_spec": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": ["z", "bptclass", "subclass"],
                "all": False,
            },
            "NEDz": {
                "use": "catsHTM",
                "rs_arcsec": 2.0,
                "keys_to_append": ["ObjType", "Velocity", "z"],
            },
            "GLADEv23": {
                "use": "extcats",
                "rs_arcsec": 2,
                "keys_to_append": ["z", "dist", "dist_err", "flag1", "flag2", "flag3"],
            },
            "LSPhotoZZou": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": [
                    "photoz",
                    "ra",
                    "dec",
                    "e_photoz",
                    "specz",
                    "_6",
                    "logMassBest",
                    "logMassInf",
                    "logMassSup",
                ],
                "pre_filter": None,
                "post_filter": None,
                "all": False,
            },
            "wiseScosPhotoz": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": [
                    "zPhoto_Corr",
                    "ra",
                    "dec",
                    "wiseID",
                    "w1mCorr",
                    "w2mCorr",
                ],
                "pre_filter": None,
                "post_filter": None,
            },
            "twoMPZ": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": ["zPhoto", "ra", "dec", "zSpec"],
                "pre_filter": None,
                "post_filter": None,
            },
            "PS1_photoz": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": [
                    "raMean",
                    "decMean",
                    "z_phot",
                    "z_photErr",
                    "z_phot0",
                    "_2",
                ],
                "pre_filter": None,
                "post_filter": None,
            },
            "PS1": {
                "use": "catsHTM",
                "rs_arcsec": 1,
            },
            "brescia": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": [],
            },
            "milliquas": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": ["broad_type", "ref_name"],
            },
            "portsmouth": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": ["sigma_stars", "sigma_stars_err", "bpt"],
            },
            "ptfvar": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": [],
            },
            "varstars": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": [],
            },
            "wise_color": {
                "use": "extcats",
                "rs_arcsec": 2.0,
                "keys_to_append": ["W1mW2"],
            },
        }
    }

    flexfit_conf = {
        "oldest_upper_limits": 14,
        "max_post_peak": 200,
    }

    directives = [
        {
            "channel": channel,
            "filter": {
                "unit": "NuclearFilter",
                "config": filter_config,
                "on_stock_match": "bypass",
            },
            "ingest": {
                "mux": {
                    "unit": "ZiArchiveMuxer",
                    "config": {"history_days": 999, "future_days": 0},
                    "combine": [
                        {
                            "unit": "ZiT1Combiner",
                            "state_t2": [
                                {
                                    "unit": "T2FlexFit",
                                    "config": flexfit_conf,
                                },
                                {
                                    "unit": "T2SimpleMetrics",
                                },
                                {
                                    "unit": "T2LightCurveSummary",
                                },
                            ],
                        }
                    ],
                    "insert": {
                        "point_t2": [
                            {
                                "unit": "T2CatalogMatch",
                                "config": cat_conf,
                                "ingest": {
                                    "filter": "PPSFilter",
                                    "sort": "jd",
                                    "select": "first",
                                },
                            },
                        ],
                    },
                }
            },
        }
    ]

    loader_config = {
        "archive": "https://ampel.zeuthen.desy.de/api/ztf/archive/v3",
        "stream": stream_token,
    }

    ac = ctx.new_context_unit(
        unit="AlertConsumer",
        process_name="AP_test",
        iter_max=1000000000,
        log_profile=os.environ.get("log_profile", "debug"),
        shaper="ZiDataPointShaper",
        compiler_opts="ZiCompilerOptions",
        supplier={
            "unit": "ZiAlertSupplier",
            "config": {
                "deserialize": None,
                "loader": {"unit": "ZTFArchiveAlertLoader", "config": loader_config},
            },
        },
        directives=directives,
    )

    n = ac.run()

    print(f"Processed {n} alerts locally (based on query output).")

    t2w = ctx.new_context_unit(
        unit="T2Worker",
        process_name="T2Worker_test",
        log_profile=os.environ.get("log_profile", "default"),
    )

    t2w.run()

    db_config: dict[str, Any] = {
        "dropbox_token": {"label": "dropbox/token"},
    }
    if not push_to_dropbox:
        db_config.update(
            {"dryRun": True, "dryRunDir": temp_dir_base, "date": date if date else None}
        )

    t3_metrics_config = db_config.copy()
    t3_ranking_config = db_config.copy()
    t3_plotneowise_config = db_config.copy()

    t3_metrics_config.update({"verbose": True})
    t3_plotneowise_config.update(
        {
            "apply_qcuts": True,
            "plot_allWISE": False,
            "verbose": True,
        }
    )

    execute_config_dict = {
        "unit": "T3ReviewUnitExecutor",
        "config": {
            "supply": {
                "unit": "T3DefaultBufferSupplier",
                "config": {
                    "select": {
                        "unit": "T3StockSelector",
                        "config": {"channel": channel},
                    },
                    "load": {
                        "unit": "T3SimpleDataLoader",
                        "config": {
                            "directives": [
                                "TRANSIENT",
                                "T2RECORD",
                                "DATAPOINT",
                                "COMPOUND",
                            ],
                            "channel": channel,
                        },
                    },
                    "complement": [
                        {"unit": "TNSReports", "config": {}},
                        {
                            "unit": "GROWTHMarshalReport",
                            "config": {},
                        },
                        {
                            "unit": "FritzReport",
                            "config": {"token": {"label": "fritz/token"}},
                        },
                    ],
                },
            },
            "stage": {
                "unit": "T3SequentialStager",
                "config": {
                    "execute": [
                        {
                            "unit": "T3MetricsPlots",
                            "config": t3_metrics_config,
                        },
                        {
                            "unit": "T3Ranking",
                            "config": t3_ranking_config,
                        },
                        {
                            "unit": "T3PlotNeoWISE",
                            "config": t3_plotneowise_config,
                        },
                    ],
                },
            },
        },
    }

    t3p = ctx.new_context_unit(
        unit="T3Processor",
        process_name="T3Processor_test",
        execute=[execute_config_dict],
    )

    t3p.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TDE filter")
    parser.add_argument(
        "-i",
        "--initiate",
        action="store_true",
        help="Initiate a new stream",
    )

    parser.add_argument(
        "-d",
        "--date",
        type=str,
        default=None,
        help="Enter a date in the form YYYY-MM-DD",
    )

    parser.add_argument(
        "--daysago",
        type=int,
        default=None,
        help="Starting from night, get the last n days",
    )

    parser.add_argument(
        "-p", "--push", action="store_true", help="NO dry run, push to dropbox instead"
    )

    args = parser.parse_args()

    run(
        initiate=args.initiate,
        date=args.date,
        daysago=args.daysago,
        push_to_dropbox=args.push,
    )
