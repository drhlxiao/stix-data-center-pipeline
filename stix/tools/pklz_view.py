import pprint
import gzip
import pickle


def view(file_in):
    print('pickle file loader')
    f = None
    if file_in.endswith('.pklz'):
        f = gzip.open(file_in, 'rb')
    elif file_in.endswith('.pkl'):
        f = open(file_in, 'rb')
    else:
        print('unknown file-type')
        return
    raw = pickle.load(f)
    print(raw)
    try:
        data = raw['packet']
        for p in data:
            print('*' * 20)
            pprint.pprint(p['header'])
            print('-' * 20)
            pprint.pprint(p['parameter'])
    except:
        pprint.pprint(raw)


import sys

if len(sys.argv[1]) < 2:
    print('pklz view <FILE>')
else:
    view(sys.argv[1])
