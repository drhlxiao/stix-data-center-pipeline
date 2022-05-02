jsonfile='/tmp/sdc.json'
host='http://localhost:5000'
url_request=host+'/request/imaging/task/last'
url_post=host+"/request/imaging/task/update"
iruns=0
obj = OBJ_NEW('IDLnetUrl')
;while 1 ne 0 do begin
	wait, 2
	json= wget(URL=url_request,/string)
	data=JSON_PARSE(json,/TOSTRUCT)
	if data.pending eq 0 then continue 

	sig_fname=data.filename
	bkg_fname=data.background.filename
	
	time_range=data.utc_range
	;filename_prefix=data.get('filename_prefix')
	aux=data.aux
	L0=aux.L0
	B0=aux.B0
	rsun=data.rsun
	roll=data.roll
	print, sig_fname, bkg_fname, time_range,  aux
	iruns++
	print, 'Executing system command'

	resp="files=/data/test.fits, /data3/test3"
	print, resp
	obj->Put(resp, /buffer, /post, url=url_post)

;end
