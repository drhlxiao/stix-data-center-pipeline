# run idl image reconstruction pipeline
#!/bin/bash
#daily routines 
cd /opt/stix/pystix

nohup python ./bin/run_imaging_server &
nohup python ./stix/flare_pipeline/plot_idl.py &

cd stix/idl/
nohup ./sswrun.csh top.pro &
