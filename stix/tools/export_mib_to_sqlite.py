"""  A script to export mib to a sqlite3 database. 
"""

import sqlite3
import glob
import os
MIB_FOLDER = './mib'
OUTPUT = 'mib.sqlite'


def start(folder=MIB_FOLDER, output=OUTPUT):
    pattern = folder + '/*.dat'

    file_list = (glob.glob(pattern))
    conn = sqlite3.connect('mib.sqlite')
    cur = conn.cursor()
    create_table = open('create_mib.sql', 'r').read()
    print('creating database')
    try:
        cur.executescript(create_table)
    except Exception as e:
        print(e)
        raise

    for fname in file_list:
        name = os.path.splitext(os.path.basename(fname))[0]

        f = open(fname, 'r')
        try:
            cursor = cur.execute('select * from ? limit 1', (name, ))
        except:
            print('Error: %s is not imported' % name)
            continue
        names = list(map(lambda x: x[0], cursor.description))
        num = len(names)
        for line in f:
            qmark = '?'
            cols = [e.strip() for e in line.split("\t")]
            if num > len(cols):
                cols.extend(['NULL'] * (num - len(cols)))

            for i in range(1, len(cols)):
                qmark += ',?'

            sql = ('''insert into {tb} values ({q})'''.format(
                tb=name, q=qmark))
            if num != len(cols):
                print(names, cols)
            else:
                cur.execute(sql, cols)

    conn.commit()


if __name__ == "__main__":
    start()
