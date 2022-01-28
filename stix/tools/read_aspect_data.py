from pymongo import MongoClient
import pprint
client = MongoClient()
db = client['stix']
packets_col = db['packets']
runs_col = db['runs']
#last_run_id=runs_col.find()#.sort({'_id':-1})#.limit(1)
#pprint.pprint(last_run_id)
#packets=runs_col.find().sort({'run_id', last_run_id})
#pprint.pprint(packets)

#print(runs_col.find({},{"_id":-1}))
run_id = 42
SPID = 54102
parameters = ['NIX00078', 'NIX00079', 'NIX00080', 'NIX00081', 'NIXD0041']
data_time = []
data = dict()
for e in parameters:
    data[e] = []

print('start')
docs = packets_col.find({'header.SPID': SPID, 'run_id': run_id})
print('processing')
for packet in docs:
    data_time.append(packet['header']['time'])
    for e in packet['parameter']:
        name = e['name']
        if name in parameters:
            data[name].extend(e['raw'])

f = open('aspect.csv', 'w')
f.write('SCET(s),Photodiode A0,Photodiode A1,Photodiode B0,Photodiode B1\n')
print('num:{}\n'.format(len(data_time)))
for i in range(0, len(data_time)):
    f.write('{},{},{},{},{}\n'.format(
        data_time[i] - data_time[0], data[parameters[0]][i],
        data[parameters[1]][i], data[parameters[2]][i],
        data[parameters[3]][i]))
print('done')
