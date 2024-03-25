PRO stx_image_reconstruction, path_bkg_file, path_sci_file, $
                           start_utc, end_utc, $
                           energy_range_lower_limit_keV, energy_range_upper_limit_keV, $
                           energy_range_full_disk_bp_map_lower_limit_keV, energy_range_full_disk_bp_map_upper_limit_keV, $
    
                           vis_fwdfit_configuration, $
                          
                           full_disk_bp_map_solo_hpc_filename,$
                           bp_map_solo_hpc_filename, $
                           vis_fwdfit_map_solo_hpc_filename, $
                           em_map_solo_hpc_filename, clean_map_solo_hpc_filename, $
                           mem_ge_map_solo_hpc_filename,$
                          
                           bp_map_stx_filename, $
                           vis_fwdfit_map_stx_filename, $
                           em_map_stx_filename, clean_map_stx_filename, $
                           mem_ge_map_stx_filename, $
                          
                           vis_filename, aux_filename, L0, B0, RSUN, roll_angle, $
                           x_offset_arcsec, y_offset_arcsec, gui

  
  clean_uniform_weighting=0
  full_disk_bp_map_size=[512,512]
  full_disk_bp_map_mapcenter=[0.,0.]
  map_size=[128,128]
  pixel_size=[2.,2.]
	full_disk_bp_map_subc_index=stx_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c'])
	subc_index=stx_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c',$
		                        '6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c'])
	clean_niter  = 200    
	clean_gain   = 0.1    
	clean_beam_width = 20.

	;;***** Parameters
	energy_range = [energy_range_full_disk_bp_map_lower_limit_keV,energy_range_full_disk_bp_map_upper_limit_keV]
	time_range   = anytim([start_utc,end_utc])
	
	;;***** Auxiliary data
	aux_data = {stx_pointing:[x_offset_arcsec,y_offset_arcsec], ROLL_ANGLE: roll_angle, L0: L0, 	B0: B0, RSUN: RSUN}
	
	;;***** Estimate flare location
	stx_estimate_flare_location_imaging_pipeline, path_sci_file, time_range, aux_data, flare_loc=flare_loc, $
	                             imsize = full_disk_bp_map_size, mapcenter = full_disk_bp_map_mapcenter, $
	                             subc_index=full_disk_bp_map_subc_index, path_bkg_file=path_bkg_file, $
	                             /silent, energy_range=energy_range, bp_map_full_disk=full_disk_bp_map
	
	;;***** Compute visibilities
	mapcenter = stx_hpc2stx_coord(flare_loc, aux_data)
	xy_flare  = mapcenter
	
	energy_range = [energy_range_lower_limit_keV,energy_range_upper_limit_keV]
	vis=stx_construct_calibrated_visibility(path_sci_file, time_range, energy_range, mapcenter, subc_index=subc_index, $
	                                        path_bkg_file=path_bkg_file, xy_flare=xy_flare, /silent)
	;;***** Save visibility and auxiliary fits file                                
	save, vis, filename = vis_filename 
	save, aux_data, filename = aux_filename 
	
	;;******* Compute the Back Projection map (around the flare location)
	bp_map = stx_bproj_imaging_pipeline(vis,map_size,pixel_size,aux_data)

	;;******* Compute the FWDFIT reconstruction (around the flare location)
	vis_fwdfit_pso_map = stx_vis_fwdfit_pso_imaging_pipeline(vis_fwdfit_configuration, vis, aux_data, imsize=map_size, pixel=pixel_size, /silent)

	;;******* Compute the CLEAN reconstruction (around the flare location)
	clean_map=stx_vis_clean_imaging_pipeline(vis,aux_data, niter=clean_niter,image_dim=map_size[0],PIXEL=pixel_size[0],uni=clean_uniform_weighting,gain=clean_gain,beam_width=clean_beam_width)
  
	;;******* Compute MEM_GE map
	total_flux = max(abs(vis.obsvis))
	mem_ge_map=stx_mem_ge_imaging_pipeline(vis,map_size,pixel_size,aux_data,total_flux=total_flux, /silent)
	
	;;******* Compute the EM reconstruction (around the flare location)
	pixel_data_summed = stx_construct_pixel_data_summed(path_sci_file, time_range, energy_range, path_bkg_file=path_bkg_file, $
	                                                    xy_flare=xy_flare, /silent, subc_index=subc_index)

	em_map = stx_em_imaging_pipeline(pixel_data_summed, aux_data, imsize=map_size, pixel=pixel_size, mapcenter=mapcenter, /silent)
	

	xy_shift=[x_offset_arcsec, y_offset_arcsec]

;;---- SAVE MAPS (SOLO HPC REFERENCE FRAME)	
	stx_map2fits_v5, full_disk_bp_map, full_disk_bp_map_solo_hpc_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift 
	stx_map2fits_v5, bp_map.solo_hpc_ref_map, bp_map_solo_hpc_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift 
	stx_map2fits_v5, vis_fwdfit_pso_map.solo_hpc_ref_map, vis_fwdfit_map_solo_hpc_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift 
	stx_map2fits_v5, clean_map.solo_hpc_ref_map, clean_map_solo_hpc_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift 
	stx_map2fits_v5, mem_ge_map.solo_hpc_ref_map, mem_ge_map_solo_hpc_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift 
	stx_map2fits_v5, em_map.solo_hpc_ref_map, em_map_solo_hpc_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift 

	;;---- SAVE MAPS (STIX REFERENCE FRAME)
	stx_map2fits_v5, bp_map.stx_ref_map, bp_map_stx_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift
	stx_map2fits_v5, vis_fwdfit_pso_map.stx_ref_map, vis_fwdfit_map_stx_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift
	stx_map2fits_v5, clean_map.stx_ref_map, clean_map_stx_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift
	stx_map2fits_v5, mem_ge_map.stx_ref_map, mem_ge_map_stx_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift
	stx_map2fits_v5, em_map.stx_ref_map, em_map_stx_filename, path_sci_file, path_bkg_file=path_bkg_file, xy_shift=xy_shift

	if (gui eq 1) then  begin
		
		chs2=1.5
		
		window,0
		cleanplot
		loadct, 3
		plot_map,full_disk_bp_map,charsize=chs2,title='Full disk Back Projection', /limb,grid_spacing=15
		
		
		window,1,xsize=1300,ysize=700
		cleanplot
		!p.multi=[0,5,2]
		loadct, 5
		
		;; SOLO HPC reference frame
		
		plot_map,bp_map.solo_hpc_ref_map,charsize=chs2,title='Back Projection',/limb,grid_spacing=5
		plot_map,em_map.solo_hpc_ref_map,charsize=chs2,title='EM',/limb,grid_spacing=5
		plot_map,vis_fwdfit_pso_map.solo_hpc_ref_map,charsize=chs2,title='VIS_FWDFIT_PSO',/limb,grid_spacing=5
		plot_map,mem_ge_map.solo_hpc_ref_map,charsize=chs2,title='MEM_GE',/limb,grid_spacing=5
		clean_maps_solo_hpc_ref = clean_map.solo_hpc_ref_map
		plot_map,clean_maps_solo_hpc_ref[0],charsize=chs2,title='CLEAN',/limb,grid_spacing=5
		
		
		plot_map,bp_map.stx_ref_map,charsize=chs2,title='Back Projection';,/limb,grid_spacing=5
		plot_map,em_map.stx_ref_map,charsize=chs2,title='EM';,/limb,grid_spacing=5
		plot_map,vis_fwdfit_pso_map.stx_ref_map,charsize=chs2,title='VIS_FWDFIT_PSO';,/limb,grid_spacing=5
		plot_map,mem_ge_map.stx_ref_map,charsize=chs2,title='MEM_GE';,/limb,grid_spacing=5
		clean_maps_solo_hpc_ref = clean_map.stx_ref_map
		plot_map,clean_maps_solo_hpc_ref[0],charsize=chs2,title='CLEAN';,/limb,grid_spacing=5
	endif
END

