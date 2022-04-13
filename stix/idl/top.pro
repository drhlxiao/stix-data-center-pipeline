;include paths for any packages we are loading

;add any other paths we need to the IDL path
!PATH=!PATH
;run user scripts
.run /home/xiaohl/FHNW/STIX/SolarFlareAnalysis/idl/hissw/stix_image_reconstruction.pro
stix_image_construct, ; parameters go here
exit
