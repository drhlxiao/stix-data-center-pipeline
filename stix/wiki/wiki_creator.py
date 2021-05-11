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
    uid=str(flare['peak_counts']).replace('.','')
    #unique_id, anti-bots
    now=datetime.now().strftime('%c')
    try:
        emph=solo.get_solo_ephemeris(peak_utc,peak_utc)
        tdiff=f"{emph['light_time_diff'][0]} s"
        angle=f"{emp['earth_sun_solo_angle'][0]} deg"
    except:
        tdiff='N/A'
        angle='N/A'
    

    title=f"STIX_Flare:_{flare['flare_id']}"
    content=(f'''[[Category:STIX_Flares]]\n'''
            f'''<!-- The session below was created by the STIX data center pipeline bot at {now}. \n '''
            f'''  Please contact hualin.xiao@fhnw.ch if you wish to contribute code to the bot. -->\n'''
            f'''==Solar Orbiter observations==\n=== STIX ===\n'''
            f"* Flare ID: {flare['flare_id']}\n"
            f"* Light travel time difference: {tdiff}\n"
            f"* Earth-Sun-SC angle: {angle}\n"
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
    f'<!-- end of section --->\n')
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
def create_wiki_for_flare(_id):
    fdb=mdb.get_collection('flares_tbc')
    doc=fdb.find_one({'_id':_id, 'hidden':False})
    if not doc:
        print(f'Flare {_id} not file in db!')
        return
    title,content=format_page(doc)
    print(content)
    create_page(title, content)
def create_wiki_for_file(file_id):
    fdb=mdb.get_collection('flares_tbc')
    flares=fdb.find({'run_id':file_id, 'hidden':False})
    if not flares:
        print(f'Flare {file_id} not file in db!')
        return
    for doc in flares:
        print(f'create wiki for {doc["_id"]}')
        title,content=format_page(doc)
        print(content)
        create_page(title, content)

if __name__=='__main__':
    if len(sys.argv)==2:
        fid=int(sys.argv[1])
        create_wiki_for_file(fid)
    elif len(sys.argv)==3:
        start_flare_id=int(sys.argv[1])
        end_flare_id=int(sys.argv[2])
        for i in range(start_flare_id, end_flare_id+1):
            create_wiki_for_file(fid)
    else:
        print('usage:\nmain <file id> \n main <start file_id> <end file_id>')
        print('file id not provided')
