;The files used for this demonstration are hosted on a STIX server and downloaded when the demo is first run
site = 'http://dataarchive.stix.i4ds.net/data/demo/ospex/'
website_url = 'https://datacenter.stix.i4ds.net/download/fits/bsd/'

;The OSPEX folder in under stx_demo_data will usually start off empty on initial installation of the STIX software
out_dir = '/home/xiaohl/FHNW/STIX/SolarFlareAnalysis/ospex_pipeline'

;if the OSPEX demo database folder is not present then create it

;As an example a spectrogram (Level 4) file for a flare on 8th February 2022 is used
;https://datacenter.stix.i4ds.net/view/list/bsd/id/7867
l4_filename = '/home/xiaohl/FHNW/STIX/SolarFlareAnalysis/ospex_pipeline/solo_L1A_stix-sci-xray-l1-2205111438_20220511T180554-20220511T184715_062294_V01.fits'

;Download the spectrogram  fits file to the stix/dbase/demo/ospex/ directory
sock_copy, site + l4_filename, status = status, out_dir = out_dir

;An observation of a non-flaring quiet time close to the flare observation can be used as a background estimate
;https://datacenter.stix.i4ds.net/view/list/bsd/id/8204
bk_filename  = 'solo_L1A_stix-sci-xray-l1-2202090020_20220209T002720-20220209T021400_036307_V01.fits'
sock_copy, site + bk_filename, status = status, out_dir = out_dir


;Now they have been dowloaded set the paths of the science data files
fits_path_data_l4 = l4_filename
;loc_file(l4_filename, path = out_dir )
fits_path_bk   = loc_file(bk_filename, path = out_dir )


;The fit interval selected for this demonstration covers one minute over the first non-thermal peak
;spex_fit_time_interval = ['11-May-2022 18:30:00.000', '11-May-2022 18:40:00.000']
spex_fit_time_interval = ['2022-05-11T18:30:00.000', '2022-05-11T18:40:00.000']

flare_start_utc = spex_fit_time_interval[0]
flare_end_utc = spex_fit_time_interval[1]

results_filename = 'fit_results.fits'
stx_auto_fit_ssw,fits_path_data = fits_path_data_l4, fits_path_bk =  fits_path_bk, flare_start_utc = flare_start_utc , $
  flare_end_utc = flare_end_utc, results_filename=results_filename

;end


