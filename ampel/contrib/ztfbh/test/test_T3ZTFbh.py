import pytest
import os
from ampel.log.AmpelLogger import AmpelLogger
import sys
logger = AmpelLogger.get_logger()
#from ampel.contrib.hu.test.test_t2brightsnprob import assert_equivalent # not needed any more
import pickle
import json
import datetime

now = datetime.datetime.now()

@pytest.fixture
def single_t3_ZTFbh_TransientView():
    with open(os.path.join(os.path.dirname(__file__), 'ZTF20aaaaflr_transientview.pkl'), 'rb') as f:
        tview = pickle.load(f)
    return tview

@pytest.fixture
def multi_t3_ZTFbh_TransientViews():
    with open(os.path.join(os.path.dirname(__file__), '50_test_tviews.pkl'), 'rb') as f:
        tviews = pickle.load(f)
    return tviews

@pytest.fixture
def secret_provider():
    from ampel.dev.DictSecretProvider import DictSecretProvider
    if not (secrets_file := os.getenv('AMPEL_SECRETS')):
        raise RuntimeError(
            "The enviroment variable AMPEL_SECRETS should point to a (possibly "
            "sops-encrypted) YAML file containing the DropBox access token"
        )
    return DictSecretProvider.load(secrets_file)


@pytest.fixture
def testdir():
    midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
    seconds_to_midnight = (midnight - now).seconds

    out_dict = {
            'testdir':now.strftime('/mampel/test/%Y-%m-%d_%H%M'),
            'force_date':now.strftime('%Y-%m-%d'),
            'force_year':now.strftime('%Y'),
            'force_md':now.strftime('%m-%d'),
            'datetime':now
                }

    if seconds_to_midnight/60<10:
        print ('Be wary: we are within {0:0.1f} minutes of midnight; testing runs for fixed date: {1}'.format(
                seconds_to_midnight60, out_dict['force_date']))

    return out_dict

@pytest.fixture
def true_output():
    # In case we want to add other comparison files in the future, unused currently
    truth = {}
    with open('ampel/contrib/ztfbh/test/ZTF20aaaaflr_wisedump.json') as f:
        truth['wise'] = json.loads(f.read())

    return truth


def test_t3_ranking(multi_t3_ZTFbh_TransientViews, secret_provider, testdir):
    
    from ampel.contrib.ztfbh.t3.T3Ranking import T3Ranking

    # don't run this right before midnight :) 
    path = '{0}/ranking/{1}/{2}_everything.txt'.format(testdir['testdir'], testdir['force_year'], testdir['force_md'])
    print ('test_t3_ranking: creating ranking files in:\n{0}'.format(path.split('_')[0]))

    sum_dir = '{0}/sum_plots/{1}/{2}'.format(testdir['testdir'], testdir['force_year'], testdir['force_md'])
    print ('test_t3_ranking: creating summary plots in:\n{0}'.format(sum_dir))

    Rank_T3 = T3Ranking(logger = logger, dryRun = False, 
                        base_location = testdir['testdir'], 
                        dropbox_token = secret_provider.get("dropbox/ztfbh", str))

    # test running with no input
    Rank_T3.add(multi_t3_ZTFbh_TransientViews[0:0])

    # test running 5 transients
    # Rank_T3.add(multi_t3_ZTFbh_TransientViews[0:5])

    # #test running in chunks
    for i in range(0,len(multi_t3_ZTFbh_TransientViews),10):
        Rank_T3.add(multi_t3_ZTFbh_TransientViews[i:i+10])
    Rank_T3.add(multi_t3_ZTFbh_TransientViews[i+10:])
    
    Rank_T3.done()

    # check that ranking and summary plot files are created
    print (Rank_T3.stats['files'])
    assert Rank_T3.stats['files'] > 6 # check all of the files were created
    


    # check for ranking    
    nlines = Rank_T3.read_file(path)[1].text.count('ZTF')
    assert nlines > 45 # check number of lines in ranking, at least 46 lines long

    # check for summary
    assert Rank_T3.exists(sum_dir+'/rise_color.pdf')
    assert Rank_T3.exists(sum_dir+'/rise_fade.pdf')
    assert Rank_T3.exists(sum_dir+'/color_change.pdf')


def test_t3_metricsplots(single_t3_ZTFbh_TransientView, secret_provider, testdir):
    
    from ampel.contrib.ztfbh.t3.T3MetricsPlots import T3MetricsPlots

    #Test MetricsPlots
    Metrics_T3 = T3MetricsPlots(logger = logger, verbose = True, dryRun = False, 
                                            base_location = testdir['testdir'],
                                            dropbox_token = secret_provider.get("dropbox/ztfbh", str))
    
    path = '{0}/alerts/2020/ZTF20aaaaflr/ZTF20aaaaflr_flex.json'.format(testdir['testdir'])
    print ('test_t3_metricsplots: creating single flexfit output at\n{0}'.format(path))
    
    Metrics_T3.add(single_t3_ZTFbh_TransientView)
    
    metrics_dump = Metrics_T3.read_file(path)[1].json()
    assert metrics_dump['name'] == "ZTF20aaaaflr"
    assert Metrics_T3.stats['files'] == 5 #ensure we save every file properly


def test_t3_plot_neowise(single_t3_ZTFbh_TransientView, secret_provider, testdir):
    
    from ampel.contrib.ztfbh.t3.T3PlotNeoWISE import T3PlotNeoWISE

    #Test PlotNeoWISE
    WISE_T3 = T3PlotNeoWISE(logger = logger, verbose = True, 
                            base_location = testdir['testdir'],
                            apply_qcuts = True, plot_allWISE = False, dryRun = False, 
                            dropbox_token = secret_provider.get("dropbox/ztfbh", str))
    
    path = '{0}/alerts/2020/ZTF20aaaaflr/ZTF20aaaaflr_neoWISE.json'.format(testdir['testdir'])
    print ('test_t3_plot_neowise: creating new neoWISE output in:\n{0}'.format(path))
    WISE_T3.add(single_t3_ZTFbh_TransientView)

    wise_dump = WISE_T3.read_file(path)[1].json()
    assert WISE_T3.stats['files'] >=3
    assert list(wise_dump.keys()) == ['info_list', 'wise_class', 'out_dict']
    #assert_equivalent(wise_dump['out_dict'], wise_truth) #don't do this, WISE is dynamic

