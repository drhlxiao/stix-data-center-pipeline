# using the chi2 minimization algorithm to find the best fit solar flare location

import math
import random
import numpy as np
from stix.flare_pipeline import skylut
from stix.spice import solo

#open area ratio obtained by  the simulator
#detectorIndex: ratio

INF=1e20
SAMPLE_POINTS=200*200
def get_chi2(area, use_small_pixels, use_det_fluence, y, yerr, cfl_cnts_sum, cfl_cnts_error_sum, det_fluence, det_fluence_error):
    area_sum = np.sum(area)
    if area_sum==0:
        return INF
    normalized_areas=area/area_sum
    det_count_chi2=0

    if use_det_fluence: #the 13th element in the  pattern
        cfl_fluence=cfl_cnts_sum/area_sum 
        cfl_fluence_error=cfl_cnts_error_sum/area_sum
        sigma=cfl_fluence_error**2+det_fluence_error**2
        #fluence counts/cm2
        det_count_chi2 = 0 if  sigma==0 else (cfl_fluence-det_fluence)**2 / sigma
    obs=y if use_small_pixels else y[0:8]
    exp=normalized_areas if use_small_pixels else normalized_areas[0:8] #expected pattern
    error=yerr if use_small_pixels else yerr[0:8]
    if np.sum(error)==0:
        error=np.sqrt(obs)
    num_elements=8 if use_small_pixels else 12
    out=np.divide((obs-exp)**2, error**2,  where=error!=0)
    return np.sum(out)+ det_count_chi2

def normalize(counts, counts_error):
    ce=np.array(counts)
    ce_err=np.array(counts_error)
    csum=np.sum(ce)
    if csum>0:
        return ce/csum, ce_err/csum
    return np.zeros(ce.size)+INF, np.zeros(ce.size)


def fit_location(counts, count_errors, mean_fluence, mean_fluence_error, flare_utc,  use_small_pixels=True, use_det_fluence=True):
    #counts: CFL pixel counts
    #count_errors: CFL pixel count errors

    #mean_fluence: detector mean fluence = pixel total counts / open area ratio
    #mean_fluence_error: detector mean fluence error



    max_pixel_index=12 if use_small_pixels else 8
    normalized_counts, normalized_count_errors=normalize(counts, count_errors)

    solo_location=solo.get_solo_ephemeris(flare_utc, flare_utc, num_steps=1)
    sun_angular_diameter=solo_location.get('sun_angular_diameter',0)
    sun_angular_radius=60
    if sun_angular_diameter:
        sun_angular_radius=sun_angular_diameter[0]*0.5*60
    #print(f'{sun_angular_diameter=}')
    #print(f'{solo_location=}')
    #in units of arsec
    x_angles = np.linspace(skylut.data['x'][0]*3600,skylut.data['x'][1]*3600,skylut.data['x'][2])
    y_angles = np.linspace(skylut.data['y'][0]*3600,skylut.data['y'][1]*3600,skylut.data['y'][2])
    dx=(skylut.data['x'][1]*3600.-skylut.data['x'][0]*3600)/skylut.data['x'][2]
    dy=(skylut.data['y'][1]*3600.-skylut.data['y'][0]*3600)/skylut.data['y'][2]
    x0=skylut.data['x'][0]*3600.
    y0=skylut.data['y'][0]*3600.

    #convert to arc min
    #contour_levels=np.array([2.3,6.17,11.8])#1,2,3 sigma
    x_steps = len(x_angles)
    y_steps = len(y_angles)
    chi2 = np.zeros((x_steps, y_steps))+INF
    margin=1.5
    min_chi2 = np.inf
    min_i=-1
    min_j=-1
            # here is caliste coordinates 
    best_fit_area=np.zeros(12)
    cfl_counts_sum=sum(counts)
    cfl_counts_error_sum=math.sqrt(sum([x**2 for x in count_errors]))
    min_x=0
    min_y=0
    #for row in skylut.data['skylut']:
    #grid search
    total_steps=x_steps*y_steps

    #print('start searching')
    #for _i in range(SAMPLE_POINTS):
    #    if _i/1000==0:
    #        print(_i)
    #    irow=random.randrange(total_steps)
    for row in skylut.data['skylut']:
        #row=skylut.data['skylut'][irow]
        x_ang=row[0]
        y_ang=row[1]
        if x_ang>margin*sun_angular_radius or x_ang<-margin*sun_angular_radius: 
            continue
        if y_ang>margin*sun_angular_radius or y_ang<-margin*sun_angular_radius: 
            continue
        area_pattern = row[2:]
        val = get_chi2(area_pattern, use_small_pixels, use_det_fluence,
                    normalized_counts, normalized_count_errors, cfl_counts_sum,
                    cfl_counts_error_sum, mean_fluence, mean_fluence_error)
        
        i=int((x_ang-x0)/dx)
        j=int((y_ang-y0)/dy)
        chi2[i][j] = val
        if val<=min_chi2:
            min_i,min_j, min_chi2,min_x,min_y=i,j, val, x_ang,y_ang
            best_fit_area=np.copy(area_pattern)


    delta_chi2 = chi2-min_chi2
    selected_indexes=np.where(delta_chi2<50)
    selected_values=delta_chi2[selected_indexes]
    indexes_list=[x.tolist() for x in selected_indexes]
    #min_pos = [x_angles[min_i], y_angles[min_j]]
    min_pos = [min_x, min_y] 
    sum_best_fit=sum(best_fit_area)
    normalized_area=[x/sum_best_fit for x in best_fit_area]
    cfl_fluence=sum(counts[0:max_pixel_index])/sum(best_fit_area[0:max_pixel_index]) 
    cfl_fluence_error=math.sqrt(sum([x*x for x in count_errors[0:max_pixel_index]]))/sum(best_fit_area[0:max_pixel_index]) 
    

    return delta_chi2, {'x': skylut.data['x'].tolist(),
            'y': skylut.data['y'].tolist(),
            'min_pos': min_pos,
            'min_chi2': min_chi2,
            'pixel_area':  best_fit_area, #before normalized
            'pixel_area_norm':  normalized_area, #normalized
            'det_mean_fluence':mean_fluence, #experimental detector mean fluence (counts/mm2)
            'det_mean_fluence_error':mean_fluence_error,  
            'cfl_fluence':cfl_fluence,  #cfl fluence (counts/mm2)
            'cfl_fluence_error':cfl_fluence_error,
            'norm_counts': normalized_counts.tolist(),
            'norm_counts_errors': normalized_count_errors.tolist(),
            'sun_angular_diameter':2*sun_angular_radius,
            'angle_unit':'arcsec',
            'delta_chi2': [indexes_list, selected_values.tolist()],
            }


def get_cfl_pattern(x_arcsec, y_arcsec, counts, counts_error):
    #calculate pattern for a source 
    def find_nearest(array, value):
        abs_val=np.abs(array-value)
        return np.where(abs_val==np.min(abs_val))[0][0]
    x_angles = np.linspace(skylut.data['x'][0]*3600,skylut.data['x'][1]*3600,skylut.data['x'][2]) #in units of arcsec
    y_angles = np.linspace(skylut.data['y'][0]*3600,skylut.data['y'][1]*3600,skylut.data['y'][2])
    x_steps = len(x_angles)
    y_steps = len(y_angles)
    i=find_nearest(x_angles, x_arcsec)
    j=find_nearest(y_angles, y_arcsec)
    #i_skylut = j*x_steps+i
    i_skylut = i*x_steps+j
    #need to be checked
    #why it is not i*x_step+j
    pattern=skylut.data['skylut'][i_skylut][2:]
    obs=np.array(counts)
    error=np.array(counts_error)
    exp=np.array(pattern[0:len(counts)])
    chi2=np.sum((obs-exp)**2/(error**2))
    result={'pattern': pattern, 'chi2':round(chi2,2)}
    return result
