import sys

import pickle
from core import stix_writer


def push(fname):
    with open(fname, 'rb') as f:
        data = pickle.load(f)
        print('loading data')
        packets = data['packets']
        run = data['run']
        st_writer = stix_writer.StixMongoWriter()
        in_filename = run['input']
        file_size = 1024
        print('pushing data')
        for i in range(500):
            comment = 'push ' + str(i)
            print(comment)
            st_writer.register_run(in_filename, file_size, comment)
            st_writer.write_all(packets)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(' filename is needed ')
    else:
        fname = sys.argv[1]
        push(fname)
