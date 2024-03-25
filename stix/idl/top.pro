;include paths for any packages we are loading

;add any other paths we need to the IDL path
!PATH=!PATH



.run /opt/stix/parser/stix/idl/stx_bproj_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_em_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_estimate_flare_location_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_make_map_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_map2fits_v5.pro
.run /opt/stix/parser/stix/idl/stx_mem_ge_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_vis_clean_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_vis_fwdfit_pso_imaging_pipeline.pro
.run /opt/stix/parser/stix/idl/stx_image_reconstruction.pro


.run /opt/stix/parser/stix/idl/spex__save_autofit_fits.pro
.run /opt/stix/parser/stix/idl/stx_auto_fit_ssw.pro
.run /opt/stix/parser/stix/idl/stx_ospex_pipeline_wrapper.pro
.run /opt/stix/parser/stix/idl/imaging_spectroscopy_daemon.pro
run_daemon(0)
exit


