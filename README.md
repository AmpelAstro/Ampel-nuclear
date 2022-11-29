# Ampel-ZTFbh
Central repository to host private AMPEL code from the ZTFbh SWG.

## Instructions
Create a fresh Python 3.10 conda env
```
conda env create -n tde_filter_upgrade python=3.10
conda activate tde_filter_upgrade
```
Install is done via poetry:
```
pip install poetry 
mkdir tde_filter_upgrade
git clone https://github.com/AmpelProject/ampel-HU-astro
cd Ampel-HU-astro
poetry install
cd ..
git clone https://github.com/robertdstein/Ampel-contrib-ZTFbh
cd Ampel-contrib-ZTFbh
pip install -e .
```
Now we have to build the ampel config. Issue
```
ampel config build -out ampel_conf.yaml
```
Now you need to export the following tokens
```
export AMPEL_ARCHIVE_TOKEN, DROPBOX_TOKEN and FRITZ_TOKEN
```
To run the test, issue
`./run_tde_test.py -i`
The `-i` initiates a new stream token. To change the date, use `-d YYYY-MM-DD` for a certain day. The script will request alerts for the 24 hours after this date.
Note: When requesting a full day with `-d` from the archive, the first run will probably fail, as it's not ready yet (`URL is locked`). Just rerun `./run_tde_test.py -d YYYY-MM-DD` in that case until it works.

To check output, go to the `temp` directory that gets created.