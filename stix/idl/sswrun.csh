#!/bin/tcsh
# Set SSW home
setenv SSW /opt/ssw
# Set SSW instruments
setenv SSW_INSTR "stix spex xray hessi"
# Setup needed environment variables
source $SSW/gen/setup/setup.ssw
# Setup IDL environment
setenv IDL_DIR /opt/idl88/idl88
setenv IDL_PATH /opt/stix/parser/stix/idl
# Startup SSW IDL
#so that idl can find the scripts
echo "Executing IDL:$1 ..."
sswidl  $1
echo "IDL existing..."
