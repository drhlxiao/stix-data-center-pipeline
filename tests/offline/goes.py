import urllib.request, json 
goes_url='https://services.swpc.noaa.gov/json/goes/primary/xrays-3-day.json'
with urllib.request.urlopen(goes_url) as url:
    data = json.loads(url.read().decode())
    print(data)
