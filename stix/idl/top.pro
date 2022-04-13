;include paths for any packages we are loading

;add any other paths we need to the IDL path
!PATH=!PATH
.run /opt/stix/parser/stix/idl/stix_image_reconstruction.pro
stx_image_reconstruct, "/data/fits/solo_L1A_stix-sci-xray-l1-2202150022_20220215T070411-20220215T075211_051272_V01.fits","/data/fits/solo_L1A_stix-sci-xray-l1-2202155513_20220215T171055-20220215T195235_052542_V01.fits","2022-02-15T18:11:49.918","2022-02-15T18:12:49.918",1,7,"sci_9111_uid_2202155513_0_0_fwfit.fits","sci_9111_uid_2202155513_0_0_bp.fits",-17.21031407198973,-3.1232662770871786,663.0281591454938,-3.004685269372137
exit
