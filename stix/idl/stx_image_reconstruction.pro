PRO stx_image_reconstruct, path_bkg_file, path_sci_file, $
	flare_start_UTC, flare_end_UTC, $
	energy_range_lower_limit_keV, energy_range_upper_limit_keV, $
	energy_range_full_disk_bp_map_lower_limit_keV, energy_range_full_disk_bp_map_upper_limit_keV, $
	full_disk_bp_map_filename,    $
	bp_map_filename, $
	vis_fwdfit_map_filename, vis_fwdfit_source_type, $
	em_map_filename, $
	clean_map_filename,   $
	L0, B0, RSUN, roll_angle, dsun, $
	x_offset_arcsec, y_offset_arcsec                     


	clean_uniform_weighting=0
	full_disk_bp_map_size=[512,512]
	full_disk_bp_map_mapcenter=[0.,0.]
	map_size=[256,256]
	pixel_size=[1.,1.]
	full_disk_bp_map_subc_index=stix_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c'])
	subc_index=stix_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c',$
		'6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c'])
	clean_niter  = 200    
	clean_gain   = 0.1    
	clean_beam_width = 20.


	;;***** Parameters
	energy_range = [energy_range_full_disk_bp_map_lower_limit_keV,energy_range_full_disk_bp_map_upper_limit_keV]
	time_range   = anytim([flare_start_UTC,flare_end_UTC])

	;;******* Compute the visibility values
	vis  = stix2vis_sep2021(path_sci_file, time_range, energy_range, full_disk_bp_map_mapcenter, subc_index=full_disk_bp_map_subc_index, $
		path_bkg_file=path_bkg_file, /silent)

	;;******* Use Giordano's (u,v) points: no need to perform projection correction (see Giordano et al., 2015)
	subc_str = stx_construct_subcollimator()
	uv       = stx_uv_points_giordano()
	u        = -uv.u * subc_str.phase
	v        = -uv.v * subc_str.phase
	vis.u    = u[full_disk_bp_map_subc_index]
	vis.v    = v[full_disk_bp_map_subc_index]

	duration=anytim2tai( flare_end_UTC) - anytim2tai(flare_start_UTC)
	;;******* Compute the Back Projection map
	pixel_size_full_disk_bp_map = RSUN * 2.6 / full_disk_bp_map_size
	full_disk_bp_map = stx_bproj(vis,full_disk_bp_map_size,pixel_size_full_disk_bp_map)
	full_disk_bp_map.xc += x_offset_arcsec
	full_disk_bp_map.yc += y_offset_arcsec

	full_disk_bp_map.time= flare_start_UTC
	full_disk_bp_map.dur=duration
	add_prop,full_disk_bp_map, DSUN= DSUN

	;;******* Compute the coordinates of the maximum value of the Back Projection map, i.e. of the location of the flare
	max_bp       = max(full_disk_bp_map.data, ind_max)
	ind_max      = array_indices(full_disk_bp_map.data, ind_max)
	max_bp_coord = [(ind_max[0]-full_disk_bp_map_size[0]/2)*pixel_size_full_disk_bp_map[0]+full_disk_bp_map_mapcenter[0], $
		(ind_max[1]-full_disk_bp_map_size[1]/2)*pixel_size_full_disk_bp_map[1]+full_disk_bp_map_mapcenter[1]]

	;;***** Re-compute the visibilities (needed for setting the map center correctly and for using subcollimators 3 to 10)
	energy_range = [energy_range_lower_limit_keV,energy_range_upper_limit_keV]
	vis  = stix2vis_sep2021(path_sci_file, time_range, energy_range, max_bp_coord, $
		xy_flare=max_bp_coord, path_bkg_file=path_bkg_file, subc_index=subc_index, /silent)

	;;******* Compute the Back Projection map (around the flare location)
	bp_map = stx_bproj(vis,map_size,pixel_size)
	bp_map.L0 = L0
	bp_map.B0 = B0
	bp_map.RSUN = RSUN
	bp_map.roll_angle = roll_angle
	bp_map.xc += x_offset_arcsec
	bp_map.yc += y_offset_arcsec



	bp_map.time= flare_start_UTC
	bp_map.dur=duration
	add_prop,bp_map, DSUN= DSUN

	;;******* Compute the FWDFIT reconstruction (around the flare location)
	vis_fwdfit_pso_map = stx_vis_fwdfit_pso(vis_fwdfit_source_type, vis, imsize=map_size, pixel=pixel_size, /silent)
	vis_fwdfit_pso_map.L0 = L0
	vis_fwdfit_pso_map.B0 = B0
	vis_fwdfit_pso_map.RSUN = RSUN
	vis_fwdfit_pso_map.roll_angle = roll_angle
	vis_fwdfit_pso_map.xc += x_offset_arcsec
	vis_fwdfit_pso_map.yc += y_offset_arcsec


	vis_fwdfit_pso_map.time= flare_start_UTC
	vis_fwdfit_pso_map.dur=duration
	add_prop,vis_fwdfit_pso_map, DSUN= DSUN

	;;******* Compute the CLEAN reconstruction (around the flare location)

	clean_map=stx_vis_clean(vis,niter=clean_niter,image_dim=map_size[0],PIXEL=pixel_size[0],uni=clean_uniform_weighting,gain=clean_gain,beam_width=clean_beam_width)

	clean_map.L0 = L0
	clean_map.B0 = B0
	clean_map.RSUN = RSUN

	clean_map.time= flare_start_UTC
	clean_map.dur=duration
	add_prop,clean_map, DSUN= DSUN

	clean_map.roll_angle = roll_angle
	clean_map.xc += x_offset_arcsec
	clean_map.yc += y_offset_arcsec


	;;******* Compute the EM reconstruction (around the flare location)
	data = stix_compute_vis_amp_phase(path_sci_file,time_range,energy_range,xy_flare=max_bp_coord,bkg_file=path_bkg_file, /silent)

	em_map = stx_em(data.RATE_PIXEL,energy_range,time_range,IMSIZE=map_size,PIXEL=pixel_size,MAPCENTER=max_bp_coord, $
		xy_flare=max_bp_coord, subc_index=subc_index, /silent)
	em_map.L0 = L0
	em_map.B0 = B0
	em_map.RSUN = RSUN


	em_map.roll_angle = roll_angle
	em_map.xc += x_offset_arcsec
	em_map.yc += y_offset_arcsec

	em_map.time= flare_start_UTC
	em_map.dur=duration
	add_prop,em_map, DSUN= DSUN

	stixmap2fits, full_disk_bp_map, full_disk_bp_map_filename
	stixmap2fits, bp_map, bp_map_filename
	stixmap2fits, vis_fwdfit_pso_map, vis_fwdfit_map_filename
	stixmap2fits, clean_map, clean_map_filename
	stixmap2fits, em_map, em_map_filename
END

