"""
  Monitor HK paramters
"""
import math
from stix.core import mongo_db as db
from stix.spice import time_utils as sdt

mdb = db.MongoDB()
MODES= {
            0: "RESET",
            1: "BOOT",
            2: "SAFE",
            3: "MAINT",
            4: "CONFIG",
            5: "NOMINAL"
        }

class HKMonitor(object):
    _limits_raw = {
        'NIXD0502': {
            'range': [20, math.inf],  # archive memory usage
            'message': 'Archive memory usage below the limit',
            'event_type':"BEYOND_LIMITS"
        },
        'NIXD0068': {
            'range': [1, 1],  # att open
            'message': 'Attenuator status changed!',
            'event_type':'ATT_STATUS_CHANGE'
        }
    }

    @staticmethod
    def scan(file_id):
        warnings = []
        limit_raw_keys = HKMonitor._limits_raw.keys()
        last_mode = None
        last_att_stat=None
        cursor = mdb.get_houskeeping_packets_of_file(file_id)
        current_att_status=None
        num=0
        for packet in cursor:
            header = packet['header']
            parameters = packet['parameters']
            current_mode = MODES[int(parameters[3][1])]
            current_time = header['unix_time']
            utc= header.get('UTC',None)
            num+=1
            if not utc:
                utc=sdt.unix2datetime(current_time)

            if current_mode != last_mode and last_mode is not None:
                doc={'time': utc,  'unix_time': current_time, 
                     'last_mode': last_mode, 'mode': current_mode, 'event_type':'MODE_CHANGE', 'file_id':file_id}
                mdb.update_hk_events(doc)
                warnings.append(
                            f"[WARNING] Mode Changed from {last_mode} to {current_mode} at {utc}"
                        )
                print(doc)
            last_mode = current_mode
            for param in packet['parameters']:
                if param[0] in limit_raw_keys:
                    if param[0]=='NIXD0068':
                        current_att_status=param[1]

                    if param[1] > HKMonitor._limits_raw[param[0]]['range'][1] or param[
                            1] < HKMonitor._limits_raw[param[0]]['range'][0]:

                        warnings.append(
                            f"[WARNING] At {utc}, {HKMonitor._limits_raw[param[0]]['message']}"
                        )

            if current_att_status!=last_att_stat and last_att_stat is not None:
                doc={'utc': utc,  'unix_time': current_time, 
                                 'ATT_STAT': 'OPEN' if current_att_status==1 else 'CLOSE', 'event_type':'ATT_STAT_CHANGE', 'file_id':file_id}
                mdb.update_hk_events(doc)
            last_att_stat=current_att_status

        print('Number of packets:', num)
        return warnings






        
        



if __name__ == '__main__':
    for i in range(200,7800):
        print("File i:",i)
        HKMonitor.scan(i)
    # test
