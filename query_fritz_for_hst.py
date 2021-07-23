import numpy as np
import os, glob
import requests, json
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy import constants as const
from astropy.table import Table, Column, MaskedColumn
import threading
import time
from tqdm import tqdm
import warnings
from datetime import date, datetime, timedelta
from astropy.time import Time
from astropy.io import ascii

today = date.today()
warnings.filterwarnings('ignore')

# Define global parameters
global GETTOKEN, BASEURL, QUERY_DAYS, QUERY_PROGRAMS, SELECT_PARAMS
with open('set_params.json') as usr:
    param_data = json.load(usr)

GETTOKEN = param_data['user']['FritzToken']
QUERY_DAYS = param_data['query_params']['saved_days']
QUERY_PROGRAMS = param_data['query_params']['programs']
SELECT_PARAMS = param_data['selection_params']
BASEURL = 'https://fritz.science/'

# Empty dict.
db = dict()
db['source'] = []

def api(method, endpoint, data=None):
    ''' Info : Basic API query, takes input the method (eg. GET, POST, etc.), the endpoint (i.e. API url)
               and additional data for filtering
        Returns : response in json format
        CAUTION! : If the query doesn't go through, try putting the 'data' input in 'data' or 'params'
                    argument in requests.request call
    '''
    headers = {'Authorization': f'token {GETTOKEN}'}

    response = requests.request(method, endpoint, json=data, headers=headers)

    return response.json()

def dist_mod_mag(app_mag, distance):
    """
    Calculate the absolute magnitude via the distance modulus.

    Input:
    -----
    app_mag (float): apparent magnitude of the source
    distance (float): distance of the source in Mpc

    Output:
    ------
    abs_mag (float): absolute magnitude
    """
    return (app_mag - 5*np.log10(distance)-25)

def redshift_to_distance(z):
    """
    Convert a given redshift to distance in Mpc, assuming H0=67.7 km/s/Mpc, and T_cmb=2.725 K
    Current cosmological constants used are similar to those on fritz.science

    Input:
    -----
    z (float): redshift of the source

    Output:
    ------
    distance (float): distance in Mpc
    """
    cosmo = FlatLambdaCDM(H0=67.66 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.30966, Neff=3.05, m_nu=[0., 0., 0.06]*u.eV, Ob0=0.049)
    return cosmo.luminosity_distance(z).value

def get_all_sources(program_id):
    """
    Fetch all the obj_id and save dates of all sources for a given program_id.
    program_id (int): program_id of your program (i.e 43 for Caltech Census of the Local Univese)
    """
    url = BASEURL + f'api/sources?group_ids={program_id}&saveSummary=true'
    #print (f"Query API URL: {url}")
    response = api('GET',url)
    return (response)

def get_source_api(ztfname, comments=False):
    ''' Info : Query a single source, takes input ZTF name
        comments (bool): If True it will also include the source comments
        Returns : all basic data of that source
    '''
    if comments:
        url = BASEURL + f'api/sources/{ztfname}?includeComments=true'
    else:
        url = BASEURL+f'api/sources/{ztfname}'

    response = api('GET', url)
    return response['data']

def get_source_photometry(ztfname, extrapolate_photometry=False):
    """ Fetch the photometry issued by ZTF for a given ZTF source.
    Params:
    -------
    ztfname (str): obj_id of the given ZTF source
    extrapolate_photometry (bool): Extrapolate the photometry from the .json file and convert it to astropy.Table format

    Output:
    ------
    astropy.Table containing mjd, filter, mag, mag_err (apparent magnitude), if extrapolate_photometry: True
    .json containing all magnitude information (including upper limits), if extrapolate_photometry: False
    """
    url = BASEURL + "api/sources/" + ztfname + "/photometry" # base-url to fetch data
    response = api("GET", url)

    if extrapolate_photometry:
        phot = response['data'] # .json file with all the photometric information
        mjd = [alert['mjd'] for alert in phot if alert['mag']!=False]
        filt = [alert['filter'] for alert in phot if alert['mag']!=False]
        mag = [alert['mag'] for alert in phot if alert['mag']!=False]
        mag_err = [alert['magerr'] for alert in phot if alert['mag']!=False]

        if len(mag):
            return (Table([mjd, filt, mag, mag_err], names=("mjd", "filter", 'mag', 'mag_err')))
        else:
            return (False)
    else:
        return (response['data'])

def fetch_z(source):
    """Fetch redshift"""

    # Fetch source data from Fritz
    # source = get_source_api(ZTF_id)

    if source['redshift']:
        return source['redshift'] # True redshift

    # There is no redshift entry on Fritz -- if in CLU, check through annotations
    else:
        # Check if there are any annotations
        if source['annotations']:
            for s in source['annotations']:
                if s['origin']=='clu:clu': # found CLU xmatch
                        z_list, z_sep_list = [], []
                        try:
                            for val in zip(s['data']['CLU_z'], s['data']['CLU_distance_arcsec']):
                                # append those redshift and seperations
                                z_list.append(val[0])
                                z_sep_list.append(val[1])

                            # Convert those redshifts and seperations into arrays
                            z_list, z_sep_list = np.array(z_list), np.array(z_sep_list)

                            # sort by the smallest seperation
                            min_to_max_z = np.argsort(z_sep_list)
                            z_list_sorted, z_sep_sorted = z_list[min_to_max_z], z_sep_list[min_to_max_z]

                            # the first index is likely the cloest match based on the redshift, but check within 30 kpc first
                            d0_kpc = z_sep_sorted[0]  * (abs(z_list_sorted[0]) / 0.05)
                            if d0_kpc<=30.:
                                return z_list_sorted[0]
                            else:
                                return None
                            
                        except:
                            return None
                else:
                    return None

def get_peak_app_mag(ztfname):
    photometry = get_source_photometry(ztfname, extrapolate_photometry=True)
    
    if photometry:
        try:
            nan_index = photometry['mag'].data!=None # select index where no None's

            filt = photometry['filter'].data[nan_index]
            Mapp = photometry['mag'].data[nan_index]
            
            peak_app = np.argmin(Mapp)

            return (Mapp[peak_app], filt[peak_app])
        
        except:
            return (None, None)
    else:
        return (None, None)

def get_phot_endpoints(source_phot, first=False, last=True):
    """Retuns the latest photometric detection alongside the band it was detected at
    source_phot is the result from get_source_photometry with extrapolate_photometry=False
    """
    mt, mf, mjde = [], [], []
    for val in source_phot:
        k0 = False
        for l in val['groups']:
            if l['name']=='Sitewide Group':
                k0 = True

        if k0==False:
            if val['mag']:
                mt.append(val['mag'])
                mf.append(val['filter'])
                mjde.append(val['mjd'])

    mt = np.array(mt) # magnitude test
    mf = np.array(mf) # magnitude test filter
    mjde = np.array(mjde) # magnitude jde

    delta_time = Time(f"{today.year}-{today.month}-{today.day}T17:00").mjd - mjde # time since today in MJD note time in in UTC
    try:
        if last:
            dt = np.argmin(delta_time)
        elif first:
            dt = np.argmax(delta_time)
            
        return (mt[dt], mf[dt], delta_time[dt])
    
    except:
        return (0, None, 0)
    
def get_age(source_phot):
    """
    Returns age, defined as the midpoint between the first detection and last non-detection, of a given ZTF source
    If the first detection is the first entry on the light curve, then age is just time since first detection
    source_phot is result of get_source_photometry with extrapolate_photometry=False
    """
    all_mjds = np.array([val['mjd'] for val in source_phot])
    all_mags = np.array([val['mag'] for val in source_phot])
    
    if len(all_mjds)>0:
        sort_by_mjd = np.argsort(all_mjds)
        
        if any(all_mags):
            first_det = (next(i for i, mag in enumerate(all_mags[sort_by_mjd]) if mag is not None))
            
            if first_det==0:
                age = Time(f"{today.year}-{today.month}-{today.day}T17:00").mjd - all_mjds[sort_by_mjd][first_det]
                return age
            else:
                age = Time(f"{today.year}-{today.month}-{today.day}T17:00").mjd - np.mean((all_mjds[sort_by_mjd][first_det-1], all_mjds[sort_by_mjd][first_det]))
                return age
        else:
            return None
    else:
        return None

def saved_in_clu_or_rcf(api_data):
    groups = []
    
    for val in api_data['groups']:
        groups.append(val['nickname'])
    
    in_clu = 'clu' in groups
    in_rcf = 'rcf' in groups
    
    return (in_clu, in_rcf)
    
def query_source(ZTF_name):
    # Fetch source infomation
    source_api_data = get_source_api(ZTF_name, comments=False) 
    source_photometry = get_source_photometry(ZTF_name, extrapolate_photometry=False)

    # Fetch redshift and classification
    param_z = fetch_z(source_api_data)
    clf = source_api_data['classifications']
    
    # source has to be classified and have a redshift so we can calculate a distance
    if param_z and clf:
        dist_Mpc = redshift_to_distance(param_z)
        
        if dist_Mpc<SELECT_PARAMS['dist_Mpc']: # HST candidates within 150 Mpc
            param_age = get_age(source_photometry)
            
            if param_age<SELECT_PARAMS['age_days']: # HST candidates younger than 2 weeks
                peak_apparent_magnitude, peak_filt = get_peak_app_mag(ZTF_name)
                
                if peak_apparent_magnitude<SELECT_PARAMS['peak_mag']: # HST candidates brighter than 19.5 apparent mag
                    peak_absolute_mag = dist_mod_mag(peak_apparent_magnitude, dist_Mpc)

                    param_classification = clf[-1]['classification']
                    param_prob = clf[-1]['probability']

                    param_in_clu, param_in_rcf = saved_in_clu_or_rcf(source_api_data)

                    param_obj_id = ZTF_name
                    param_ra = source_api_data['ra']
                    param_dec = source_api_data['dec']
                    last_mag, last_filt, last_det = get_phot_endpoints(source_photometry, first=False, last=True) # last detection
                    first_mag, first_filt, first_det = get_phot_endpoints(source_photometry, first=True, last=False) # first detection

                    db['source'].append({"ZTF_id": param_obj_id, "ra":param_ra, "dec":param_dec, 
                                         "z":param_z, "d_Mpc":dist_Mpc, 'age':param_age,
                                         "classification": param_classification, "classification_prob": param_prob,
                                         'peak_abs_mag': peak_absolute_mag, "peak_app_mag": peak_apparent_magnitude, 'peak_mag_filt':peak_filt,
                                         'last_mag': last_mag, 'last_filt':last_filt, 'last_det':last_det, 
                                         'first_mag':first_mag, 'first_filt':first_filt, 'first_det':first_det,
                                         'in_clu':param_in_clu, 'in_rcf':param_in_rcf})
                else:
                    pass
            else:
                pass
        else:
            pass
    else:
        pass
    
def main():
    program_names, program_ids = zip(*QUERY_PROGRAMS.items())
    
    print(f'Searching sources in saved in the last {QUERY_DAYS} days from {program_names}...')
    date_sort = Time(datetime.now() - timedelta(days=QUERY_DAYS)) # only query sources saved in the last query_from weeks
    
    # get names 
    print('Loading ZTF IDs...')
    for i, program_id in enumerate(program_ids):
        sources = get_all_sources(program_id) # fetch all sources for a given program ID (CLU: 43, RCF: 41)

        # All sources (names & dates) of all saved candidates in the program specified by program_id
        names = np.array([s['obj_id'] for s in sources['data']['sources']])
        dates = np.array([s['saved_at'] for s in sources['data']['sources']])
        w0_t = np.where(dates>=date_sort)
        
        if i==0: # if first program we're querying, need to initialize all_names array
            all_names = names[w0_t]
            print(f'{len(all_names)} sources in {program_names[i]}')
        else:
            names = names[w0_t]
            all_names = np.concatenate((all_names, names))
            print(f'{len(names)} sources in {program_names[i]}')
    
    all_names = np.unique(all_names) # removes duplicate names
    print(f'{len(all_names)} unique sources')
    print('')
    print('Querying Fritz...')
    
    list_threads = [] # list of threads for pooling
    buff = 0 # buffer time (NOTE: API sometimes crashes if request is to frequent)
    for i in tqdm(range(len(all_names))):
        buff +=1
        if buff>=5:
            time.sleep(5) 
            buff = 0 # restart buffer

        t = threading.Thread(target=query_source, args=[all_names[i]])
        t.start()
        list_threads.append(t)

    for thread in list_threads:
        thread.join()
        
    param_names = ("ZTF_id", "ra", "dec", "z", "d_Mpc", 'age',
                   "classification", "classification_prob",
                   'peak_abs_mag', "peak_app_mag", 'peak_mag_filt',
                   'last_mag', 'last_filt', 'last_det', 
                   'first_mag', 'first_filt', 'first_det',
                   'in_clu', 'in_rcf')
        
    cols = []
    for i, param_name in enumerate(param_names):
        param_vals = np.asarray([s[param_name] for s in db['source']])
        none_vals = np.isin(param_vals, [None, '', 'None']) # test if any of the None/no value flags exist, return True at index if it does
            
        # if there's a bad value, make a mask and save as astropy MaskedColumn
        if np.any(none_vals):
            col = MaskedColumn(param_vals, name=param_name, mask=none_vals)
        # if all the values are fine, save as regular astropy Column
        else:
            col = Column(param_vals, name=param_name)
                
        cols.append(col)
            
    data_table = Table(cols)
    data_table.sort('age')
    
    print('')
    print(f'Checking for supernova...')
    
    SNs = np.array([i for i, clf in enumerate(data_table['classification']) if ('Ia' in clf) or ('Ib' in clf) or ('Ic' in clf) or ('II' in clf) or ('Ca-rich' in clf)])
    cands = data_table[SNs]
    
    print(f'Found {len(cands)} supernova in {len(data_table)} candidates for HST rolling snapshots program as of {today.year}-{today.month}-{today.day}: ')
        
    for i in range(len(cands)):
        print(f"{cands[i]['ZTF_id']}: {cands[i]['classification']}")
        
    print('')
    print('Remember to double-check source pages!')
    
    data_table.write(f"candidate_data/HST_Candidates_{today.year}_{today.month}_{today.day}.ascii", format='ascii')
    
    print('')
    
if __name__=='__main__':
    main()