# run idl image reconstruction pipeline
#!/bin/bash
#daily routines 
cd /opt/stix/pystix
./bin/run_imaging_server

cd stix/idl/
./sswrun.csh top.pro
