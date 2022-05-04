;include paths for any packages we are loading

;add any other paths we need to the IDL path
!PATH=!PATH
.run /opt/stix/parser/stix/idl/stx_image_reconstruction.pro
.run /opt/stix/parser/stix/idl/imaging_daemon.pro
.run /opt/stix/parser/stix/idl/stixmap2fits.pro
run_daemon(0)
exit
