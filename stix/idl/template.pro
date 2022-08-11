;A idl script to reproduce STIX preview images 
;Created at 2022-05-06T16:06:32.856475 by STIX data center online image reconstruction software

;To run this script, sswidl and stix idl software must be installed on your computer.  
; Download FITS files from STIX data center
wget("https://datacenter.stix.i4ds.net/download/fits/filename?uid=2202155513&filename=solo_L1A_stix-sci-xray-l1-2202150022_20220215T070411-20220215T075211_051272_V01.fits",filename="solo_L1A_stix-sci-xray-l1-2202150022_20220215T070411-20220215T075211_051272_V01.fits") ; background file,
wget("https://datacenter.stix.i4ds.net/download/fits/filename?uid=2202155513&filename=solo_L1A_stix-sci-xray-l1-2202155513_20220215T171055-20220215T195235_052542_V01.fits", filename="solo_L1A_stix-sci-xray-l1-2202155513_20220215T171055-20220215T195235_052542_V01.fits" ; signal file,
; Uncomment the following two line if you don't have stix_image_reconstruction.pro and stixmap2fits.pro on your local disk,
;wget("https://datacenter.stix.i4ds.net/pub/misc/stix_imaging/stx_image_reconstruction.pro", filename="stx_image_reconstruction.pro"),
;wget("https://datacenter.stix.i4ds.net/pub/misc/stix_imaging/stixmap2fits.pro", filename="stixmap2fits.pro"),

sig_fname="/data/fits/solo_L1A_stix-sci-xray-l1-2202155513_20220215T171055-20220215T195235_052542_V01.fits"
bkg_fname="/data/fits/solo_L1A_stix-sci-xray-l1-2202150022_20220215T070411-20220215T075211_051272_V01.fits"

time_range=["2022-02-15 17:58:42.7144","2022-02-15 18:23:47.5432"]
; time range for image reconstruction

path_sci_file=""
path_bkg_file=""
start_utc="2022-02-15 17:58:42.7144"
end_utc="2022-02-15 18:23:47.5432"

bp_elow=6 ; back-projection energy range lower limit
bp_ehigh=10
; energy range in  units of keV, used to make a back-project full image
; the result will be used to locate the source(s)

elow=13.0
ehigh=28.0
; used to make other images


L0=-17.217203594422983

B0=-3.1216889460075663
RSUN=1325.640342183433
; apparent radius  of the sun in  units of arcsec
dsun=108246894279.66695
;distance between the sun and s/c in units of meters

roll_angle=-3.006307334797977
; spacecraft roll angle in units of degrees. Extracted from SPICE kernel data 

x_offset_arcsec=-0.07648256123004528
y_offset_arcsec=0.15199631645222014
;Note that aspect system correction was not applied when calculating the pointing offsets.


vis_fwdfit_source_type='circle'
;Change the source shape if necessary. The source shape can be also "ellipse" or "multi"  (multi-circle).


bp_fname="bp_map.fits" 
full_disk_bp_fname="full_disk_bp_map.fits" 
vis_fwdfit_fname= "vis_fwdfit_map.fits" 
em_fname= "em_map.fits"
clean_fname="clean_map.fits"
;Output filenames

stx_image_reconstruct, path_bkg_file, path_sci_file, $
	start_utc, end_utc, $
	elow, ehigh, $
	bp_elow, bp_ehigh, $
	full_disk_bp_fname,  $
	bp_fname, $
	vis_fwdfit_fname, vis_fwdfit_source_type, $
	em_fname, $
	clean_fname,  $
	L0, B0, RSUN, roll_angle, dsun, $
	x_offset_arcsec, y_offset_arcsec     

end
