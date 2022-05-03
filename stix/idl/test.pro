;include paths for any packages we are loading

;add any other paths we need to the IDL path
!PATH=!PATH
.run stx_image_reconstruction.pro
.run imaging_daemon.pro
run_daemon(0)
exit
