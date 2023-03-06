#!/bin/bash
#daily routines 
cd /opt/stix/pystix
echo "####################################################">>sync.log
date>>sync.log
echo `pwd`
echo "changing directory to /opt/stix/pystix"  >>sync.log
rsync -av hxiao@pub099.cs.technik.fhnw.ch:/var/www/data/fits/L1/* /data/pub099/fits/L1/
rsync -av hxiao@pub099.cs.technik.fhnw.ch:/var/www/data/fits/L2/* /data/pub099/fits/L2/ --exclude="*hk*"
cd /opt/stix/pystix

echo "run daily pipeline sub tasks..."
nohup /usr/bin/python3 ./bin/stixpl_daily.py &
#download 
echo "Importing aspect to data..."
nohup /usr/bin/python3 ./stix/data_imports/data_archive.py &
sleep 10
echo "starting pipeline to recreate fits files for old data " >> sync.log
nohup /usr/bin/python3 ./stix/fits/fits_creator.py -br 14 13 &
sleep 3600
kill -9 `pgrep -f  "find_no_fits_bsd.py"` >>sync.log
nohup /usr/bin/python3 stix/find_no_fits_bsd.py &
echo "Done">>sync.log

sleep 3600
nohup /usr/bin/python3 ./stix/flare_pipeline/imaging_task_manager.py &
