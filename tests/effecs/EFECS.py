#!/usr/bin/env python
# coding: utf-8

# In[57]:


f='EFECS_M04_V02.xml'
import sys
import xmltodict
from dateutil import parser as dtparser
from datetime import datetime, timedelta
fd=open(f)
doc = xmltodict.parse(fd.read())
#import 
for t in doc['eventfile']['events']['PASS']:
    start=t['@time']
    #sdt=dtparser.parse(start)
    sdt=datetime.strptime(start,'%Y-%jT%H:%M:%SZ')
    duration=t['@duration']
    end=sdt+timedelta(seconds=int(duration))
    
    print(sdt, end)
#print(x)


# In[17]:


print(req)


# In[ ]:




