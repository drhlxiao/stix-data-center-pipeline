#!/bin/bash
#daily routines 
cd /opt/stix/pystix
echo "####################################################">>sync.log
date>>sync.log
echo `pwd`
echo "changing directory to /opt/stix/pystix"  >>sync.log
rsync -av hxiao@pub099.cs.technik.fhnw.ch:/var/www/data/fits/L1/* /data/pub099/fits/L1/
rsync -av hxiao@pub099.cs.technik.fhnw.ch:/var/www/data/fits/L2/* /data/pub099/fits/L2/ --exclude="*hk*"


