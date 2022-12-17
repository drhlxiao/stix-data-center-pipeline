#!/bin/bash
#daily routines 
#rsync -av hxiao@pub099.cs.technik.fhnw.ch:/var/www/data/fits/L1/* /data/pub099/fits/L1/
rsync -av hxiao@pub099.cs.technik.fhnw.ch:/var/www/data/fits/L2/* /data/pub099/fits/L2/ --exclude="*hk*"
date>>sync.log
/usr/bin/python3 /opt/stix/parser/stix/data_imports/data_archive.py
/usr/bin/python3 /opt/stix/parser/stix/fits/fits_creator.py -br 14 13
/usr/bin/python3 /opt/stix/parser/stix/tools/find_nofits_bsdc.py

