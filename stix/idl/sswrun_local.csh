#!/bin/tcsh
# Set SSW home
setenv SSW /data2/ssw
# Set SSW instruments
setenv SSW_INSTR "stix"
# Setup needed environment variables
source $SSW/gen/setup/setup.ssw
# Setup IDL environment
setenv IDL_DIR /data2/idl/idl88
# Startup SSW IDL
#sswidl /opt/stix/parser/stix/idl/top.pro
echo "Executing IDL:$1 ..."
sswidl  $1
echo "IDL existing..."
