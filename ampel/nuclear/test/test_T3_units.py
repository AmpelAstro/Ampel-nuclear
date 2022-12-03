#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, pickle, json, datetime

import pytest

from ampel.log.AmpelLogger import AmpelLogger
from ampel.secret.AmpelVault import AmpelVault
from ampel.secret.DictSecretProvider import DictSecretProvider

logger = AmpelLogger.get_logger()

now = datetime.datetime.now()


@pytest.fixture
def single_TransientView_fixture():
    with open(
        os.path.join(os.path.dirname(__file__), "ZTF20aaaaflr_transientview.pickle"),
        "rb",
    ) as f:
        tview = pickle.load(f)
    return tview


@pytest.fixture
def multi_TransientView_fixture():
    with open(
        os.path.join(os.path.dirname(__file__), "transientview_list.pickle"), "rb"
    ) as f:
        tview_list = pickle.load(f)
    return tview_list


@pytest.fixture
def dropbox_token():
    from ampel.secret.NamedSecret import NamedSecret

    return NamedSecret(label="dropbox/token", value=os.environ["DROPBOX_TOKEN"])


@pytest.fixture
def testdir():
    midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
    seconds_to_midnight = (midnight - now).seconds

    out_dict = {
        "testdir": now.strftime("/mampel/test/%Y-%m-%d_%H%M"),
        "force_date": now.strftime("%Y-%m-%d"),
        "force_year": now.strftime("%Y"),
        "force_md": now.strftime("%m-%d"),
        "datetime": now,
    }

    if seconds_to_midnight / 60 < 10:
        print(
            "Be wary: we are within {0:0.1f} minutes of midnight; testing runs for fixed date: {1}".format(
                seconds_to_midnight60, out_dict["force_date"]
            )
        )

    return out_dict


@pytest.fixture
def true_output():
    # In case we want to add other comparison files in the future, unused currently
    truth = {}
    with open("ampel/nuclear/test/ZTF20aaaaflr_wisedump.json") as f:
        truth["wise"] = json.loads(f.read())

    return truth


def test_t3_ranking(multi_TransientView_fixture, dropbox_token, testdir):

    from ampel.nuclear.t3.T3Ranking import T3Ranking

    # don't run this right before midnight :)
    # path = "{0}/ranking/{1}/{2}_everything.txt".format(
    #     testdir["testdir"], testdir["force_year"], testdir["force_md"]
    # )
    path = f"{testdir['testdir']}/ranking/{testdir['force_year']}/{testdir['force_md']}_everything.txt"

    print("test_t3_ranking: creating ranking files in:\n{0}".format(path.split("_")[0]))

    sum_dir = (
        f"{testdir['testdir']}/sum_plots/{testdir['force_year']}/{testdir['force_md']}"
    )

    print("test_t3_ranking: creating summary plots in:\n{0}".format(sum_dir))

    Rank_T3 = T3Ranking(
        logger=logger,
        dryRun=False,
        base_location=testdir["testdir"],
        dropbox_token=dropbox_token,
    )
    Rank_T3.post_init()

    # test running with no input
    Rank_T3.process([])

    # need to init again
    Rank_T3.post_init()

    # test running all transients
    Rank_T3.process(multi_TransientView_fixture)

    # check that ranking and summary plot files are created
    print(Rank_T3.stats["files"])
    assert Rank_T3.stats["files"] > 6  # check all of the files were created

    # check for ranking
    nlines = Rank_T3.read_file(path).text.count("ZTF")
    print(Rank_T3.read_file(path).text)
    print(nlines)
    assert nlines > 45  # check number of lines in ranking, at least 46 lines long

    # check for summary
    assert Rank_T3.exists(sum_dir + "/rise_color.pdf")
    assert Rank_T3.exists(sum_dir + "/rise_fade.pdf")
    assert Rank_T3.exists(sum_dir + "/color_change.pdf")


def test_t3_metricsplots(single_TransientView_fixture, dropbox_token, testdir):

    from ampel.nuclear.t3.T3MetricsPlots import T3MetricsPlots

    # Test MetricsPlots
    Metrics_T3 = T3MetricsPlots(
        logger=logger,
        verbose=True,
        dryRun=False,
        base_location=testdir["testdir"],
        dropbox_token=dropbox_token,
    )
    Metrics_T3.post_init()

    path = f"{testdir['testdir']}/alerts/2020/ZTF20aaaaflr/ZTF20aaaaflr_flex.json"

    print("test_t3_metricsplots: creating single flexfit output at\n{0}".format(path))

    Metrics_T3.process([single_TransientView_fixture])

    metrics_dump = Metrics_T3.read_file(path).json()
    assert metrics_dump["name"] == "ZTF20aaaaflr"
    assert Metrics_T3.stats["files"] == 5  # ensure we save every file properly


def test_t3_plot_neowise(single_TransientView_fixture, dropbox_token, testdir):

    from ampel.nuclear.t3.T3PlotNeoWISE import T3PlotNeoWISE

    # Test PlotNeoWISE
    WISE_T3 = T3PlotNeoWISE(
        logger=logger,
        verbose=True,
        base_location=testdir["testdir"],
        apply_qcuts=True,
        plot_allWISE=False,
        dryRun=False,
        dropbox_token=dropbox_token,
    )
    WISE_T3.post_init()

    path = f"{testdir['testdir']}/alerts/2020/ZTF20aaaaflr/ZTF20aaaaflr_neoWISE.json"

    print("test_t3_plot_neowise: creating new neoWISE output in:\n{0}".format(path))

    WISE_T3.process([single_TransientView_fixture])

    wise_dump = WISE_T3.read_file(path).json()
    assert WISE_T3.stats["files"] >= 3
    assert list(wise_dump.keys()) == ["info_list", "wise_class", "out_dict"]
