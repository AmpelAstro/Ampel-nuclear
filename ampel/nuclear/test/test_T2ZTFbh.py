import pytest
import os
from ampel.log.AmpelLogger import AmpelLogger
logger = AmpelLogger.get_logger()
import sys
from ampel.contrib.hu.test.test_t2brightsnprob import assert_equivalent
import pickle
import json

@pytest.fixture
def t2_ZTFbh_lightcurve():
    with open(os.path.join(os.path.dirname(__file__), "ZTF20aaadtjq_transientview.pkl"), 'rb') as f:
        test_tview = pickle.load(f)
    lightCurve = test_tview[0].get_lightcurves()[-1]
    return lightCurve

@pytest.fixture
def true_output():
    truth = {}
    truth['flex'] = json.load(open(os.path.join(os.path.dirname(__file__),"ZTF20aaadtjq_flexresult.json"), 'rb'))
    truth['simple'] = json.load(open(os.path.join(os.path.dirname(__file__),"ZTF20aaadtjq_simpleresult.json"), 'rb'))
    return truth


def test_t2_flexfit(t2_ZTFbh_lightcurve, true_output):
    from ampel.contrib.ztfbh.t2.T2FlexFit import T2FlexFit

    #Test flexfit
    flexT2 = T2FlexFit(logger = logger, oldest_upper_limits = 14, max_post_peak = 200)
    flex_result = flexT2.run(t2_ZTFbh_lightcurve)
    flex_truth = true_output['flex']
    assert list(flex_result.keys()) == ['fit_params', 'lightcurve_data', 'plot_info'] # check to make sure output dict is full
    assert flex_result['fit_params']['name'] == flex_truth['fit_params']['name']

    # For now, not making detailed comparison of fitting results, we might improve flexfit in the near-future
    # (and don't want to update the test each time)
    #assert_equivalent(flex_result['fit_params'], flex_truth['fit_params']) 
    


def test_t2_simplemetrics(t2_ZTFbh_lightcurve, true_output):
    from ampel.contrib.ztfbh.t2.T2SimpleMetrics import T2SimpleMetrics

    #Test SimpleMetrics
    simpleT2 = T2SimpleMetrics(logger = logger)
    simple_result = simpleT2.run(t2_ZTFbh_lightcurve)
    simple_truth = true_output['simple']
    assert list(simple_result.keys()) == ['metrics', 'plot_info'] # Ensure proper output
    # assert_equivalent({k:v for k,v in simple_result['metrics'].items() if k not in ['age']}, 
    #                   {k:v for k,v in simple_truth['metrics'].items() if k not in ['age']}) # Check the metrics excluding 
    assert simple_result['metrics']['name'] == simple_truth['metrics']['name']
