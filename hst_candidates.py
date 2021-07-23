import query_fritz_for_hst as query_fritz
from astropy.io import ascii
import numpy as np
from natsort import natsorted
import os
from flask import Flask, render_template
import webbrowser
from threading import Timer
from datetime import date

app = Flask(__name__) # initiate Flask so we can open a table of candidates in the browser and click through directly to Fritz

today = date.today()
files = natsorted(os.listdir('candidate_data'))
if len(files)>0:
    file = files[-1]
    print (f"Most recent candidate search: {file}")
else:
    file = None
    print("No recent searches.")
    

if file:
    load_new = input('Continue using this file? [y/n] ')
    if load_new=='y':
        print(f"Loading data from {file}")
        data = ascii.read(f'candidate_data/{file}')
    elif load_new=='n':
        print('Running Fritz querying script...')
        query_fritz.main()
        file = f'HST_Candidates_{today.year}_{today.month}_{today.day}.ascii'
        print(f"Loading data from {file}")
    data = ascii.read(f'candidate_data/{file}')
    print("")
    
else:
    print("Running Fritz querying script...")
    query_fritz.main()
    file = f'HST_Candidates_{today.year}_{today.month}_{today.day}.ascii'
    print(f"Loading data from {file}")
    data = ascii.read(f'candidate_data/{file}')
    print("")

def open_browser():
      webbrowser.open_new('http://127.0.0.1:5000/')

@app.route("/")
def load_webpage():
    theta = data[np.argsort(data['age'])]
    data_date = f'{int(file.split("_")[2])}-{int(file.split("_")[3])}-{int(file.split("_")[4].split(".")[0])}'
    return render_template('candtable.html', posts=theta, date=data_date)

if __name__=='__main__':
    Timer(1, open_browser).start()
    #app.run(port=5000, debug=True)
    app.run(port=5000, debug=True, use_reloader=False)