#!/bin/tcsh
# Set SSW home
setenv SSW /opt/ssw
# Set SSW instruments
setenv SSW_INSTR "stix"
# Setup needed environment variables
source $SSW/gen/setup/setup.ssw
# Setup IDL environment
setenv IDL_DIR /opt/idl88/idl88
# Startup SSW IDL
sswidl /opt/stix/parser/stix/idl/top.pro
