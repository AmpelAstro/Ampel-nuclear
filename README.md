# Ampel-nuclear
Central repository to host AMPEL code to search for and analyze nuclear transients. At the moment, this is exclusively code from the ZTFbh science working group.

## Installation
### Prerequisites
You need to export environment variables for the [AMPEL ZTF archive](https://ampelproject.github.io/astronomy/ztf/index) (tokens are available [here](https://ampel.zeuthen.desy.de/live/dashboard/tokens)), for [Fritz](https://fritz.science/), and for the dropbox API (ask Sjoert). 

Furthermore, you need a running instance of [MongoDB](https://www.mongodb.com/docs/manual/installation/). On macOS, make sure you have the command line tools installed (in doubt, run `xcode-select –install`).

### Setup
Create a fresh Python 3.10 conda env
```
conda create -n tde_filter_upgrade python=3.10
conda activate tde_filter_upgrade
```
Install is done via poetry:
```
pip install poetry 
git clone https://github.com/AmpelProject/ampel-nuclear
cd Ampel-nuclear
poetry install
```
Now we have to build the ampel config. Issue
```
ampel config build -out ampel_conf.yaml
```
Note: this will throw some import errors, but you can ignore these because those packages are not needed locally. 

Now you need to export the following tokens
```
export AMPEL_ARCHIVE_TOKEN='' 
export DROPBOX_TOKEN=''
export FRITZ_TOKEN=''
```

## Test
To run the test, start your local MongoDB. And then issue

```
./run_tde_scan.py -i
```

If you cannot execute the file, issue `chmod +x run_tde_scan.py`.

Note: To push the result of a run to the dropbox, add `-p`.

The `-i` initiates (and saves) a new archive API stream token. To request one day, use `-d YYYY-MM-DD` for a certain day. The script will request alerts for the 24 hours after this date.

Your can also use `--daysago n` to scan the last `n` days. 

Note: When requesting a full day with `-d` or the last `n` days with `--daysago n` from the archive, the first run will probably fail, as the archive database has not fully ramped up yet (`URL is locked`). In this case, just rerun `./run_tde_scan.py` (without any parameters to prevent requesting a new stream token and overwriting the current one) until the archive starts serving alerts (you will see them getting ingested).

To check the output, go to the `temp` directory that gets created when script is run without `-p` (push to dropbox), or check the dropbox.

To see all available commands of the test script, run `./run_tde_scan.py -h`.

### Examples
```
./run_tde_scan.py -i --daysago 4 -p

```

This will perform a search for the last 4 days and push the results to the dropbox.
```
./run_tde_scan.py -i -d 2022-10-06

```
This will perform a search for October 6, 2022 and save the result in a local directory.
