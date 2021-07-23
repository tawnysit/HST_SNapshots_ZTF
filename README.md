## Scripts to search for supernova candidates in ZTF for the HST Rolling Snapshots program

#### Getting Started
---
You will need to add your Fritz credentials to `set_params.json`:
```
"user":
    {"user_last_name": "USER_LAST_NAME_HERE",
     "FritzToken": "FRITZ_AUTHENTICATION_TOKEN_HERE"
    },
```
The last name is there but not strictly necessary to run the script, but you will need your command-line authentication token.

#### Adjusting Search Parameters
---
By default, the script searches all candidates saved to the CLU and RCF programs in the last 15 days. These parameters can be adjusted by editing `saved_days` and `programs` in `set_params.json`:
```
"query_params":
    {"saved_days": NUMBER_OF_DAYS_BEFORE_TODAY_TO_QUERY,
     "programs": {"CLU":43, "RCF":41, PROGRAM_NAME:PROGRAM_ID}
    },
```
The current criteria for selecting a source as a HST Rolling Snapshots candidate is that the source must be within 150 Mpc, brighter than 19.5 in apparent magnitude, and younger than 2 weeks. These parameters can also be adjusted by editing the `selection_params` in `set_params.json`:
```
"selection_params":
    {"dist_Mpc": 150.0,
     "age_days": 14,
     "peak_mag": 19.5
    }
```
Note that you will have to go into `templates/candtable.html` to adjust the description of selection criteria if these parameters are adjusted. The same file can also be edited to display different parameters in the table on the local webpage.

#### Running the Scripts
---
After editing `set_params.json`, just navigate to the main directory and run 
```
python hst_candidates.py
```
from the command line. If there is already an `.ascii` data file of candidates in the directory `candidate_data`, it will prompt you to ask if you would like to use this file `[y]` or search again `[n]`. If you choose to search, or there is no previous data, the script automatically runs `query_fritz_for_hst.py`. Then, a local webpage should open via Flask, displaying a table with candidates, direct links to source pages on Fritz, and some important information.

If you would not like to open the webpage, you can also just run
```
python query_fritz_for_hst.py
```
from the command line. This will simply output an `.ascii` file with candidate information in the directory `candidate_data`.

#### Additional Notes
---
If you would like to search for non-SN objects, eg. TDE, you will have to edit the filter in line 383 of the `main()` function of `query_fritz_for_hst.py` directly. Currently, the script searches for classifications that include the substrings `Ia`, `Ib`, `Ic`, `II`, or `Ca-rich` which should encompass all available SN classifications on Fritz.
