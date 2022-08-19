pro stx_auto_fit_ssw,fits_path_data = fits_path_data, fits_path_bk =  fits_path_bk, flare_start_utc = flare_start_utc , $
  flare_end_utc = flare_end_utc, distance = distance, time_shift = time_shift, require_nonthermal = require_nonthermal, results_filename = results_filename 

  spex_fit_time_interval = [flare_start_utc, flare_end_utc]


  default, require_nonthermal, 0

  stx_get_header_corrections, fits_path_data, distance = header_distance, time_shift = header_time_shift
  default, distance, header_distance
  default, time_shift, header_time_shift
  

  ;For Reading the STIX specific spectrum and response matrix files
  spex_file_reader= 'stx_read_sp'

  ;As the STIX instrument is still being fully calibrated an uniform systematic uncertainty of 5% is applied here
  spex_uncert= 0.0500000

  emax_thermal = 15
  use_thermal_only = 0
  min_thick_norm = 1e-10
  max_nt_index = 9
  max_lec = 50

  nan = !values.f_nan


  !null = mrdfits(fits_path_data, 0, primary_header)
  orig_filename = sxpar(primary_header, 'FILENAME')
  rate_str = stx_read_fits(fits_path_data, 1, rate_header)
  request_id = strtrim(rate_str.request_id, 2)

  if strpos(orig_filename, 'cpd') gt -1 or strpos(orig_filename, 'xray-l1') gt -1 then begin
    stx_convert_pixel_data, fits_path_data = fits_path_data, $
      fits_path_bk =  fits_path_bk, $
      distance = distance, $
      time_shift = time_shift,$
      ospex_obj = ospex_obj, $
      plot = 0
  endif else if strpos(orig_filename, 'spec') gt -1 or strpos(orig_filename, 'spectrogram') gt -1 then begin
    stx_convert_spectrogram, fits_path_data = fits_path_data, $
      fits_path_bk =  fits_path_bk, $
      distance = distance, $
      time_shift = time_shift,$
      ospex_obj = ospex_obj, $
      plot = 0
  endif else begin
    message, 'ERROR: the FILENAME field in the primary header should contain either cpd, xray-l1 or spec'
  endelse
   
  set_logenv, 'OSPEX_NOINTERACTIVE', '1'

  ;set the values as defined in section 2 above for this object
  ospex_obj-> set, spex_file_reader = spex_file_reader
  ospex_obj-> set, spex_fit_time_interval = spex_fit_time_interval
  ospex_obj-> set, spex_uncert = spex_uncert
  ospex_obj-> set, spex_fit_manual=0, spex_fit_reverse=0
  ospex_obj-> set, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0, spex_fit_progbar=0
  ospex_obj-> set, fit_function= 'vth'
  ospex_obj-> set, fit_comp_params = [0.01, 1.0, 1.000]
  ospex_obj-> set, spex_fit_auto_erange = 1
  ospex_obj-> set, fit_comp_free_mask = [1B, 1B, 0B]
  ospex_obj-> set, fit_comp_spectrum = ['full']
  ospex_obj-> set, fit_comp_model = ['chianti']
  ospex_obj-> set, mcurvefit_itmax= 100
  ospex_obj-> dofit, /all
  energy_range = ospex_obj-> get(/spex_erange)
  ; if (energy_range[1] gt 28) and  (energy_range[1] le 36) then energy_range[1] = 28
  fit_param_thermal_full = ospex_obj->get(/fit_comp_params)
 
  
  ospex_obj-> set, fit_function= 'vth+thick2'
  ospex_obj-> set, fit_comp_params = [0.01, 1.0, 1.000,  1e-10, 5., 1e+6, 6.0, 15., 3.2e+4]
  ospex_obj-> set, fit_comp_free_mask = [1B, 1B, 0B, 0B, 0B, 0B, 0B, 0B, 0B]
  ospex_obj-> set, fit_comp_spectrum = [ 'full', '']
  ospex_obj-> set, fit_comp_model = [ 'chianti', '']
  ospex_obj-> set, spex_fit_auto_erange = 0
  ospex_obj-> set, spex_erange = [4,emax_thermal]
  ospex_obj->set, mcurvefit_itmax=1000
  ospex_obj->set, mcurvefit_tol=1e-4
  ospex_obj-> dofit

  nchan = n_elements(ospex_obj->get(/spex_summ_energy))/2


  fit_param_thermal_low = ospex_obj->get(/fit_comp_params)
  summ_vth = ospex_obj -> get(/spex_summ)
  vth_func_components = ospex_obj -> calc_func_components(spex_units ='flux')

  mil = 1
  start_time_string =  time2fid( atime((summ_vth.spex_summ_time_interval)[0]), /time, mil = mil)
  start_time_string = start_time_string.substring(0,-3)
  start_time_string = start_time_string.replace('_','T')

  end_time_string =  time2fid( atime((summ_vth.spex_summ_time_interval)[1]), /time, mil = mil)
  end_time_string = end_time_string.substring(0,-3)
  end_time_string = end_time_string.replace('_','T')

 ; results_filename = 'fit_results_'+start_time_string+'-'+end_time_string+'.fits'


  ospex_obj-> set, spex_erange = [emax_thermal,energy_range[1]]
  ospex_obj-> set, fit_comp_free_mask = [0B, 0B, 0B, 1B, 1B, 0B, 0B, 1B, 0B]
  ospex_obj-> set, fit_comp_params = [fit_param_thermal_low[0:2], 0.1, 5., 1e+6, 6.0, 15., 3.2e+4]
  ospex_obj-> dofit
  ospex_obj-> set, spex_erange = [4,energy_range[1]]
  ospex_obj-> set, fit_comp_free_mask = [1B, 1B, 0B, 1B, 1B, 0B, 0B, 1B, 0B]
  ospex_obj-> dofit

  vth_thick_func_components = ospex_obj -> calc_func_components(spex_units ='flux')

  summ_vth_thick = ospex_obj -> get(/spex_summ)
  params_vth_thick = summ_vth_thick.spex_summ_params

  chisq_vth = summ_vth.spex_summ_chisq
  chisq_vth_thick = summ_vth_thick.spex_summ_chisq


  chi1 = summ_vth.spex_summ_chisq
  dof1 =  total(summ_vth.spex_summ_emask)- total(summ_vth.spex_summ_free_mask)

  chi2 = summ_vth_thick.spex_summ_chisq
  dof2 =  total(summ_vth_thick.spex_summ_emask) - total(summ_vth_thick.spex_summ_free_mask)
  f = (chi1/dof1) / (chi2/dof2)

  prob = mpftest(f, dof1, dof2, /cleve)

  if emax_thermal ge energy_range[1] $
    or chisq_vth le chisq_vth_thick $
    or params_vth_thick[3] lt min_thick_norm $
    or params_vth_thick[4] gt max_nt_index $
    or params_vth_thick[7] gt max_lec $
    or require_nonthermal $
    then use_thermal_only = 1

  if use_thermal_only then begin
    fit_function= 'vth'
    free_mask= [1B, 1B, 0B, 0B, 0B, 0B, 0B, 0B, 0B]
    params = fltarr(9,/no)+nan
    params[0:2] = fit_param_thermal_low
    resid = summ_vth.spex_summ_resid
    sigmas = summ_vth.spex_summ_sigmas
    convfac  =  summ_vth.spex_summ_conv
    phmodel  =  summ_vth.spex_summ_ph_model
    startpar =  summ_vth.spex_summ_starting_params
  endif 
  
  chisq       = use_thermal_only ? chisq_vth                        : chisq_vth_thick
  model_total = use_thermal_only ? (vth_func_components.yvals)[*,0] : (vth_thick_func_components.yvals)[*,0]
  model_vth   = use_thermal_only ? (vth_func_components.yvals)[*,1] : (vth_thick_func_components.yvals)[*,1]
  model_thick = use_thermal_only ? fltarr(nchan)                    : (vth_thick_func_components.yvals)[*,2]

  goodness_of_fit  = chisq le 2. ? 'Acceptable Fit' : 'Poor Fit'

  fit_out = { $
    model_total : model_total, $
    model_vth   : model_vth, $
    model_thick : model_thick, $
    goodness_of_fit : goodness_of_fit $
  }

  resolve_routine,'spex::fitswrite', /either, /compile_full_file
  ospex_obj->save_autofit_fits, fit_out, outfile = results_filename

  if use_thermal_only then begin
    rate = mrdfits(results_filename,1)
    fxhmodify, results_filename, 'FITFUNC', 'vth', 'Fit function used            ',EXTENSION = 1
    rate.chisq     = chisq
    rate.convfac   = convfac
    rate.phmodel   = phmodel
    rate.residual  = resid
    rate.freemask  = free_mask
    rate.params    = params
    rate.startpar  = startpar
    rate.sigmas    = sigmas
    modfits, results_filename, rate, exten_no = 1
  endif


  specpath = ospex_obj->get(/spex_specfile)
  drmpath = ospex_obj->get(/spex_drmfile)

  if file_exist(specpath) then file_delete, specpath
  if file_exist(drmpath) then file_delete, drmpath
  obj_destroy, ospex_obj
  
  set_logenv, 'OSPEX_NOINTERACTIVE', '0'
  stop
end
