from pprint import pprint
f = open('stix/data/MIB/css.dat')
lines = f.readlines()
output = {}
last_seq = ''

for l in lines:
    line = l.split('\t')
    if line[0] not in output:
        output[line[0]] = []
    timestamp = line[8].strip().split('.')
    #print(line[8], timestamp)

    delta_time = ''
    if len(timestamp) == 3:
        delta_time = 3600 * int(timestamp[0]) + 60 * int(timestamp[1]) + int(
            timestamp[2])

    output[line[0]].append(
        (line[4], line[1], line[2], line[7], line[8], delta_time))
pprint(output)
