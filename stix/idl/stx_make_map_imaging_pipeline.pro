;+
;
; NAME:
;
;   stx_make_map
;
; PURPOSE:
;
;   Create a map from a STIX reconstruction (SOLO HPC coordinate). Observed visibilities used to reconstruct the map
;   and visibilities predicted from the map are saved in the map structure
;
; CALLING SEQUENCE:
;
;   stx_map = stx_make_map(im_map, aux_data, pixel, method, vis)
;
; INPUTS:
;
;   im_map: image reconstruted with a STIX imaging method (STIX coordinate frame)
;
;   aux_data: auxiliary data structure containing information on STIX pointing and spacecraft ephemeris
;
;   pixel: bi-dimensional array containing the pixel size in arcsec
;
;   method: string indicating which method has been used for reconstructing the image passed in im_map
;
;   vis: visibility structure used for reconstructing an image
;
; OUTPUTS:
;
;   Map structure. The following fields are added:
;
;     - ENERGY_RANGE: array containing the lower and upper edge of the energy interval used for reconstructing 'im_map' (keV)
;
;     - OBS_VIS: visibility structure used for reconstructing 'im_map'
;
;     - PRED_VIS: structure containing the complex values of the visibilities predicted from 'im_map'. Used for displaying data
;                 fitting plots and for computing the chi2 value associated with 'im_map'
;
;     - AUX_DATA: auxiliary data srtucture used for reconstructing 'im_map'
;
;     - RSUN: apparent radius of the SUN (arcsec)
;
;     - L0: Heliographic longitude (degrees)
;
;     - B0: Heliographic latitude (degrees)
;
;     - COORD_FRAME: string indicating the coordinate system of the output map (SOLO HPC)
;
;     - UNITS: string indicating the units of reconstructed map (photons cm^-2 asec^-2 s^-1 keV^-1)
;
; HISTORY: September 2022, Massa P., created
;
; CONTACT:
;   paolo.massa@wku.edu
;-

function stx_make_map_imaging_pipeline, im_map, aux_data, pixel, method, vis

  imsize = size(im_map, /dim)
  
  ;;----------- MAP 1: conceived in the STIX reference frame
  map1 = make_map(im_map)

  id_map = 'STX ' + method + ': '
  map1.ID = id_map

  map1.dx = pixel[0]
  map1.dy = pixel[1]

  time_range     = vis[0].TIME_RANGE
  this_time_range=stx_time2any(time_range,/vms)

  map1.time = anytim((anytim(this_time_range[1])+anytim(this_time_range[0]))/2.,/vms)
  map1.dur  = anytim(this_time_range[1])-anytim(this_time_range[0])

  ;; Mapcenter coordinates
  mapcenter = vis[0].xyoffset

  map1.xc = mapcenter[0]
  map1.yc = mapcenter[1]

  ;; Add properties
  energy_range   = vis[0].ENERGY_RANGE

  add_prop, map1, energy_range   = energy_range          ; energy range in keV
  add_prop, map1, time_range     = time_range            ; Time range considered for image reconstruction
  add_prop, map1, aux_data       = aux_data 
  add_prop, map1, rsun = aux_data.RSUN
  add_prop, map1, b0   = aux_data.B0
  add_prop, map1, l0   = aux_data.L0
  add_prop, map1, coord_frame = 'STIX reference frame'
  add_prop, map1, units = 'counts cm^-2 asec^-2 s^-1 keV^-1'
  
  ;;----------- MAP 2: conceived in the STIX reference frame

  map2 = make_map(rotate(im_map,1))

  id_map = 'STX ' + method + ': '
  map2.ID = id_map

  map2.dx = pixel[0]
  map2.dy = pixel[1]

  time_range     = vis[0].TIME_RANGE
  this_time_range=stx_time2any(time_range,/vms)

  map2.time = anytim((anytim(this_time_range[1])+anytim(this_time_range[0]))/2.,/vms)
  map2.dur  = anytim(this_time_range[1])-anytim(this_time_range[0])

  ;; Mapcenter coordinates
  mapcenter = stx_rtn2stx_coord(vis[0].xyoffset, aux_data, /inverse)

  map2.xc = mapcenter[0]
  map2.yc = mapcenter[1]
  map2=rot_map(map2,-aux_data.ROLL_ANGLE,rcenter=[0.,0.])
  map2.roll_angle = 0.

  ;; Add properties
  energy_range   = vis[0].ENERGY_RANGE

  add_prop, map2, energy_range   = energy_range          ; energy range in keV
  add_prop, map2, time_range     = time_range            ; Time range considered for image reconstruction
  add_prop, map2, aux_data       = aux_data 
  add_prop, map2, rsun = aux_data.RSUN
  add_prop, map2, b0   = aux_data.B0
  add_prop, map2, l0   = aux_data.L0
  add_prop, map2, coord_frame = 'SOLO HPC'
  add_prop, map2, units = 'counts cm^-2 asec^-2 s^-1 keV^-1'

  maps = {stx_ref_map: map1, $
          solo_hpc_ref_map: map2}
  
  return, maps

end
