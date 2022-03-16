"""
  Monitor HK paramters
"""
import math
from stix.core import mongo_db as db

mdb = db.MongoDB()
class HKMonitor(object):
    _limits_raw={
            'NIXD0502':{
                'range': [20, math.inf], #archive memory usage
                'message': 'Archive memory usage below the limit'
                },
            'NIXD0068':{
                'range': [1,1], #att open
                'message': 'Attenuator status changed!'
            }
            }
    
    @staticmethod
    def scan(file_id):
        warnings=[]
        limit_raw_keys=HKMonitor._limits_raw.keys()
        cursor=mdb.get_houskeeping_packets_of_file(file_id)
        for packet in cursor:
            for param in packet['parameters']:
                if param[0] in limit_raw_keys:
                    if param[1] > _limits_raw[param[0]]['range'][1] or param[1]< _limits_raw[param[0]]['range'][0]:
                        warnings.append(f"[WARNING] At {packet['header']['utc']}, {_limits_raw[param[0]]['message']}")
        return warnings


if __name__=='__main__':
    msg=HKMonitor.scan(1029)
    print(msg)
    #test






                    

        


