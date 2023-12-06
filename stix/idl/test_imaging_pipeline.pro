
;;********** Download files

out_dir = "/Users/admin/Documents/scientific_projects/STIX/data_center/data"

website_url = 'https://datacenter.stix.i4ds.net/download/fits/bsd/'

; UID of the science fits file to be dowloaded from the website
uid_sci_file = "1178428688"
; Download the science fits file (if not already stored in out_dir)
sock_copy, website_url + uid_sci_file, out_name, status = status, out_dir = out_dir, $
  local_file=path_sci_file, clobber=0

; UID of the background fits file to be dowloaded from the website
uid_bkg_file = "1178082832"
; Download the background fits file (if not already stored in out_dir)
sock_copy, website_url + uid_bkg_file, out_name, status = status, out_dir = out_dir, $
  local_file=path_bkg_file, clobber=0

; URL of the server containing the L2 auxiliary fits files
website_url = 'http://dataarchive.stix.i4ds.net/fits/L2/'
; Filename of the auxiliary L2 fits file to be downloaded
file_name    = '2020/06/07/AUX/solo_L2_stix-aux-ephemeris_20200607_V01.fits'
; Download the L2 auxiliary fits file (if not already stored in out_dir)
sock_copy, website_url + file_name, out_name, status = status, out_dir = out_dir, $
  local_file=path_aux_fits_file, clobber=0


;;********** Define parameters

start_utc='2020-06-07T21:39:00.000000'
end_utc  ='2020-06-07T21:42:49.000000'

; Parameters used for full-disk image reconstruction
bp_elow=6 ; back-projection energy range lower limit in units of keV
bp_ehigh=10 ; back-projection energy range upper limit in units of keV

; Energy range for EM, BP and forward-fit algorithms 
elow=6.0
ehigh=10.0

vis_fwdfit_configuration='ellipse'
;Change the source shape type if necessary, source shape type options:  'circle', 'ellipse', 'loop' or combinations 
;(e.g., ['circle','circle'], ['circle','ellipse'], etc.)


bp_fname="bp_map.fits" 
full_disk_bp_fname="full_disk_bp_map.fits" 
vis_fwdfit_fname= "vis_fwdfit_map.fits" 
em_fname= "em_map.fits"
clean_fname="clean_map.fits"
;Output filenames

	
stx_image_reconstruct, path_bkg_file, path_sci_file, path_aux_fits_file, $
  start_utc, end_utc, $
  elow, ehigh, $
  bp_elow, bp_ehigh, $
  full_disk_bp_fname,    $
  bp_fname, $
  vis_fwdfit_fname, vis_fwdfit_configuration, $
  em_fname, $
  clean_fname, 1  

end
