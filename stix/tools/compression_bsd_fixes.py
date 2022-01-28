import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
packet_db= db['packets']
cursor=packet_db.find({'header.SPID':54115})
#cursor=packet_db.find({'_id':771287})

print('number of packets:',cursor.count())
from stix.core import decompressor
from stix.core import datatypes as std

s=0
k=4
m=4
trig_parameters=[
            'NIX00242',
            'NIX00243',
            'NIX00244',
            'NIX00245',
            'NIX00246',
            'NIX00247',
            'NIX00248',
            'NIX00249',
            'NIX00250',
            'NIX00251',
            'NIX00252',
            'NIX00253',
            'NIX00254',
            'NIX00255',
            'NIX00256',
            'NIX00257'
            ]




num=0

for pkt in cursor:
    packet=std.Packet(pkt)
    s=packet.get_one('NIXD0007')[1]
    k=packet.get_one('NIXD0008')[1]
    m=packet.get_one('NIXD0009')[1]

    packet.get_one('NIXD0011')[1]=k
    packet.get_one('NIXD0012')[1]=m
    #print(s,k,m)
    #branch=pkt[13][3]
    for param in  pkt['parameters'][13][3]:
        if param[0] in trig_parameters:
            raw=param[1]
            #origin=int(param[2])
            param[2]=decompressor.decompress(raw, s, k, m)
            #print(origin,'=>', param[2])

    num += 1
    packet_db.replace_one({'_id':pkt['_id']}, pkt)
    print('fixing:', num)
            


    

    
    
