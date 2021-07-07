import os
import sqlite3
import threading
from dateutil import parser as dtparser
filename="grb.sqlite"
con=None
cur=None
try:
    con = sqlite3.connect(filename)
    cur = con.cursor()
except sqlite3.Error as e:
    print("Error %d: %s" % (e.args[0],e.args[1]))


def excute_sql(sql):
    print('executing sql...')
    cur.execute(sql)
    print('commit...')
    con.commit()


def contains(a, b):
    if b in a:
        return True
    return False
def append_grb_info(det, trigid, entryid):
    more=", %s %s)"%(det,trigid)
    sql='update burst_events set comments=concat(comments,"%s") where id=%d'%(more, entryid)
    excute_sql(sql)



def grb_exists(det, trigid, burst_time):
    name=f"{det} {trigid}"
    sql="select * from burst_events where name='{name}'"
    res=fetchall(sql)
    return res
    

def excute_fetchall(sql):
    cur.execute(sql)
    res=cur.fetchall()
    return res
def fetchall(sql):
    return excute_fetchall(sql)

def main():
    sql='select start_utc, name  from  burst_events'
    rows=fetchall(sql)
    for row in rows:
        ut=dtparser.parse(row[0]).timestamp()
        print(ut)
if __name__=='__main__':
    main()
