#create or edit an wiki page
#don't upload it to github
import sys
sys.path.append('.')
import requests
from datetime import datetime
from stix.core import mongo_db as db
from stix.spice import solo
mdb = db.MongoDB()
S = requests.Session()
HOST= "https://pub023.cs.technik.fhnw.ch"
URL=HOST+"/wiki/api.php"


def format_page(flare):
    peak_utc=flare['peak_utc']
    emph=solo.get_solo_ephemeris(peak_utc,peak_utc)
    tdiff=emph.get('light_time_diff','NA')+' s'
    sun_solo_r=emph.get('sun_solo_r','NA')+' AU'
    uid=str(flare['peak_counts']).replace('.','')
    #unique_id, anti-bots
    now=datetime.now().strftime('%c')
    

    title=f"STIX_Flare:_{flare['flare_id']}"
    content=(f'''[[Category:STIX_Flares]]\n'''
            f'''<!-- The session below was created by the STIX data center pipeline bot at {now}. \n '''
            f'''  Please contact hualin.xiao@fhnw.ch if you wish to contribute code to the bot. --> '''
            f'''==Solar Orbiter observations==\n=== STIX ===\n'''
            f"* Flare ID: {flare['flare_id']}\n"
            f"* Light travel time difference: {tdiff}\n"
    f"* Flare Peak UTC: {flare['peak_utc']}\n"
    f'* STIX QL light curves\n'
    f'<img src="{HOST}/request/image/flare?id={flare_id}&type=stixlc&uid={uid}&p=0"></img>\n'
    f'<img src="{HOST}/request/image/flare?id={flare_id}&type=loc&uid={uid}&p=0"></img>\n'
    f'===EUI===\n'
    f'<img src="{HOST}/request/image/flare?id={flare_id}&type=eui&uid={uid}&p=0"></img>\n'
    f'===EPD===\n'
    f'<img src="{HOST}/request/image/flare?id={flare_id}&type=epd&uid={uid}&p=0"></img>\n'
    f'==GOES X-ray flux==\n'
    f'<img src="{HOST}/request/image/flare?id={flare_id}&type=goes&uid={uid}&p=0"></img>\n'
    f'==AIA==\n'
    f'<img src="{HOST}/request/image/flare?id={flare_id}&type=aia&uid={uid}&p=0"></img>\n'
    f'<!-- end of first section --->\n')

    )
    return title, content

def create_page(title, content):
    # Step 1: GET request to fetch login token
    PARAMS_0 = {
        "action": "query",
        "meta": "tokens",
        "type": "login",
        "format": "json"
    }

    R = S.get(url=URL, params=PARAMS_0)
    DATA = R.json()
    print(DATA)

    LOGIN_TOKEN = DATA['query']['tokens']['logintoken']

    # Step 2: POST request to log in. Use of main account for login is not
    # supported. Obtain credentials via Special:BotPasswords
    # (https://www.mediawiki.org/wiki/Special:BotPasswords) for lgname & lgpassword
    PARAMS_1 = {
        "action": "login",
        "lgname": "stix_data_center",
        "lgpassword": "gostix2021",
        "lgtoken": LOGIN_TOKEN,
        "format": "json"
    }

    R = S.post(URL, data=PARAMS_1)

    # Step 3: GET request to fetch CSRF token
    PARAMS_2 = {
        "action": "query",
        "meta": "tokens",
        "format": "json"
    }

    R = S.get(url=URL, params=PARAMS_2)
    DATA = R.json()
    print(DATA)

    CSRF_TOKEN = DATA['query']['tokens']['csrftoken']
    print(CSRF_TOKEN)

    # Step 4: POST request to edit a page
    PARAMS_3 = {
        "action": "edit",
        "bot":"true",
        "title": title,
        "token": CSRF_TOKEN,
        "format": "json",
        "recreate":"true",
        "appendtext":content 

    }

    R = S.post(URL, data=PARAMS_3)
    DATA = R.json()
    print(DATA)
    PARAMS_3['captchaid']=DATA['edit']['captcha']['id']
    PARAMS_3['captchaword']='stix'
    R = S.post(URL, data=PARAMS_3)
    DATA = R.json()
    print(DATA)
def create_flare_wiki(_id):
    fdb=mdb.get_collection('flares_tbc')
    doc=fdb.find_one({'_id':_id, 'hidden':False})
    print(doc)
    if doc:
        title,content=format_page(doc)
        print(content)
        create_page(title, content)
if __name__=='__main__':
    if len(sys.argv)==2:
        flare_id=int(sys.argv[1])
        create_flare_wiki(flare_id)
