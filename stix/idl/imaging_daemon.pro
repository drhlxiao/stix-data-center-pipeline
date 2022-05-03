FUNCTION run_daemon, i
host='http://localhost:8022'
url_request=host+'/request/imaging/task/last'
url_post=host+"/request/imaging/task/update"
obj = IDLnetUrl()
i=0

WHILE (1 ne 0) DO BEGIN
	print, 'Loop '+i+'...'
	i++
	print, url_request
	json= wget(url_request,/string_array)
	data=JSON_PARSE(json,/TOSTRUCT)
	wait, 1
	i++
	if (data.pending eq 0) then  begin
		wait, 10
		continue
	endif

	sig_fname=data.filename
	bkg_fname=data.background.filename

	time_range=data.utc_range
	;filename_prefix=data.get('filename_prefix')

	path_sci_file=data.filename
	path_bkg_file=data.background.filename

	start_utc=data.utc_range[0]
	end_utc=data.utc_range[1]

	bp_elow=6
	bp_ehigh=10

	elow=data.energy_range[0]
	ehigh=data.energy_range[1]

	out_folder=data.idl_config.folder+'/'+data.idl_config.prefix
	aux=data.aux
	L0=aux.L0
	B0=aux.B0
	RSUN=aux.rsun
	roll_angle=aux.roll
	x_offset_arcsec=aux.sun_center[0]
	y_offset_arcsec=aux.sun_center[1]    



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

	vis_fwdfit_source_type= data.idl_config.fwdfit_shape ; multi or ellipse


	bp_fname=out_folder + 'bp_map.fits' 
	full_disk_bp_fname=out_folder + 'full_disk_bp_map.fits' 
	vis_fwdfit_fname=out_folder + 'vis_fwdfit_map.fits' 
	em_fname=out_folder + 'em_map.fits'
	clean_fname=out_folder + 'clean_map.fits'

	print, bp_fname+','+full_disk_bp_fname

	clean_uniform_weighting=0
	stx_image_reconstruct, path_bkg_file, path_sci_file, $
		print, path_bkg_file, path_sci_file, $
		start_utc, end_utc, $
		elow, ehigh, $
		bp_elow, bp_ehigh, $
		full_disk_bp_fname, full_disk_bp_map_size, full_disk_bp_map_subc_index, full_disk_bp_map_mapcenter, $
		map_size, pixel_size, subc_index, $
		bp_fname, $
		vis_fwdfit_fname, vis_fwdfit_source_type, $
		em_fname, $
		clean_fname, clean_niter, clean_gain, clean_beam_width, clean_uniform_weighting, $
		L0, B0, RSUN, roll_angle, $
		x_offset_arcsec, y_offset_arcsec     
	resp="_id="+string(data._id)+"&image_bp="+bp_fname+"&vis_fwdfit="+vis_fwdfit_fname+"&image_em="+em_fname+"&image_clean="+clean_fname+"&image_full_disk="+full_disk_bp_fname
	
	ret=obj->Put(resp, /buffer, /post, url=url_post)
	print, "done"
ENDWHILE 
return, 1
END




