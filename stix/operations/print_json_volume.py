"""
  print data volume information in json files 
  Feb 22, 2022
  Hualin Xiao
"""
import glob
import json
print('Filename, QL Estimated volume, Conservative ')
for fname in  glob.glob('*json'):
    f=open(fname)
    jf=json.load(f)
    vol_a, vol_b=jf['predicted_total_volume'],jf["total_volume_upper_limit"]
    print(f'{fname}, {vol_a/(1024):.0f}, {vol_b/1024:.0f}')
