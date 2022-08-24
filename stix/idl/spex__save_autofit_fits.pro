pro spex::fit_results_auto_ext, fit_out, fptr, rate_struct, ext, sdate, version
  ; Get header related information
  fptr->setheader, hdr, /init
  self->add_primary_keywords, fptr, rate_struct, sdate, version
  fptr->Set, create = 0
  fptr->Set, modify = 1
  fptr->Set, extension= ext
  ; Add keywords to the  header
  fptr->Addpar, 'EXTNAME', 'Fit Components', 'Extension name'
  fptr->Addpar, 'ORIGIN', 'Unknown', 'Origin of FITS file'
  fptr->Addpar, 'INSTRUME', 'STIX', 'Instrument name'
  fptr->Addpar, 'TIMVERSN', 'OGIP/93-003', 'OGIP memo number where the convention used'
  ;fptr->Addpar, 'VERSION', '0.2.0', 'File format version number'
  currtime = strmid(anytim(!stime, /ccsds), 0, 19)
  fptr->Addpar, 'DATE', currtime, 'File creation date (YYYY-MM-DDThh:mm:ss UTC)'
  fit_out = str_sub2top(fit_out, err_msg=err_msg, err_code=err_code)
  fptr->write, fit_out

end


pro spex::save_autofit_fits,fit_out, outfile=outfile

  version = '1.0'  ; version number of file format

  summ = self->get(/spex_summ)

  if (n_elements(summ.spex_summ_time_interval) gt 1) then begin
    if not is_string(outfile) then begin
      ;       atim = strlowcase( trim( str_replace(anytim(!stime, /vms, /date), '-', '_') ))
      outfile = 'ospex_fit_results_' + time2file( (summ.spex_summ_time_interval)[0], /sec ) + '.fits'
    ENDIF


    IF outfile ne '' THEN BEGIN
      dummy = self->get(/spex_summ_energy)
      nchan = n_elements(dummy)/2   ; spex_summ_energy is 2-d array

      ; Create a structure with keyword values - used for RATE extension
      timeInterval = self->get(/spex_summ_time_interval)
      ntimes = n_elements(timeInterval)/2   ; spex_summ_time_interval is 2-d array

      sdate = anytim(/ccsds, minmax(anytim( summ.SPEX_SUMM_TIME_INTERVAL, /sec )))

      rate_struct = self->defineRateKeys(timeInterval, version)
      rate_struct.detchans = nchan
      rate_struct.fitFunction = self->get(/spex_summ_fit_function)
      rate_struct.area = self->get(/spex_summ_area)


      ; Create a fitswrite object
      fptr = self->createFits(rate_struct, sdate, version, outfile, header)

      ; Create a header with fit results, add keywords and write to file
      self->fit_results_write, fptr, header, rate_struct, 1, sdate, timeInterval, version

      ; Create ebounds extension, add keywords and write to file
      self->createEboundsExt, fptr, rate_struct, 2, sdate, version, 0

      ; Create OSPEX control parameter extension
      self->fit_results_auto_ext, fit_out, fptr, rate_struct, 3, sdate, version

      obj_destroy, fptr

      message, 'Fit results saved in file ' + outfile, /info



    ENDIF ELSE message,'No output file name selected.  Not saving results.'

  ENDIF ELSE BEGIN
    message, 'No Fit Results to save.  Aborting.'
    outfile = ' '
  ENDELSE

end



