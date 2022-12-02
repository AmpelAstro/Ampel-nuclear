#!/usr/bin/env python

import os, sys, pickle, json

import pytest
from ampel.log.AmpelLogger import AmpelLogger

logger = AmpelLogger.get_logger()

# import pickle

# with open(
#     "/Users/simeon/ampel-nuclear/ampel/nuclear/test/test_tv.pickle", "wb"
# ) as f:
#     pickle.dump(transients[0], f)


@pytest.fixture
def lightCurve_fixture():
    with open(
        os.path.join(os.path.dirname(__file__), "ZTF20aaadtjq_transientview.pickle"),
        "rb",
    ) as f:
        tv = pickle.load(f)

    print(tv)
    lightCurve = tv.get_lightcurves()[0]

    return lightCurve


@pytest.fixture
def true_output():
    truth = {}
    truth["flex"] = json.load(
        open(
            os.path.join(os.path.dirname(__file__), "ZTF20aaadtjq_flexresult.json"),
            "rb",
        )
    )
    truth["simple"] = json.load(
        open(
            os.path.join(os.path.dirname(__file__), "ZTF20aaadtjq_simpleresult.json"),
            "rb",
        )
    )
    return truth


def test_t2_flexfit(lightCurve_fixture, true_output):
    from ampel.nuclear.t2.T2FlexFit import T2FlexFit

    # test T2FlexFit
    flexT2 = T2FlexFit(logger=logger, oldest_upper_limits=14, max_post_peak=200)
    flex_result = flexT2.process(lightCurve_fixture)
    flex_truth = true_output["flex"]

    assert list(flex_result.keys()) == [
        "fit_params",
        "lightcurve_data",
        "plot_info",
    ]  # check to make sure output dict is full
    assert flex_result["fit_params"]["name"] == flex_truth["fit_params"]["name"]


def test_t2_simplemetrics(lightCurve_fixture, true_output):
    from ampel.nuclear.t2.T2SimpleMetrics import T2SimpleMetrics

    #     # test T2SimpleMetrics
    simpleT2 = T2SimpleMetrics(logger=logger)
    simple_result = simpleT2.process(lightCurve_fixture)
    simple_truth = true_output["simple"]
    assert list(simple_result.keys()) == [
        "metrics",
        "plot_info",
    ]  # Ensure proper output

    assert simple_result["metrics"]["name"] == simple_truth["metrics"]["name"]
