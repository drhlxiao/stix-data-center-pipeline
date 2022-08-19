FUNCTION stx_ospex_pipeline_wrapper, fits_path_data , fits_path_bk, flare_start_utc , flare_end_utc, require_nonthermal,   results_filename
	catch, error
    if error ne 0 then begin
		PRINT, 'OSPEX Error index: ',error 
        print, 'OSPEX pipeline error: ' + !error_state.msg
        catch, /cancel
		return, ''
    endif

	stx_auto_fit_ssw,fits_path_data = fits_path_data, fits_path_bk =  fits_path_bk, flare_start_utc = flare_start_utc, $
	  flare_end_utc = flare_end_utc, time_shift=0, require_nonthermal=require_nonthermal, results_filename=results_filename
	return, results_filename
END
