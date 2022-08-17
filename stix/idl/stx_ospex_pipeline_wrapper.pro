pro stx_ospex_pipeline_wrapper,fits_path_data = fits_path_data, fits_path_bk =  fits_path_bk, flare_start_utc = flare_start_utc , $
  flare_end_utc = flare_end_utc,  results_filename=results_filename
	catch, error
    if error ne 0 then begin
        catch, /cancel
        print, 'OSPEX pipeline error: ' + !error_state.msg
        return, []
    endif


	stx_auto_fit_ssw,fits_path_data = path_sci_file, fits_path_bk =  path_bkg_file, flare_start_utc = start_utc, $
	  flare_end_utc = end_utc, time_shift=0, results_filename=spectral_fitting_results_filename

end
