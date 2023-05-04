"""
STIX science data products
"""
from datetime import timedelta, datetime
from itertools import chain

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table.operations import unique, vstack
from astropy.table.table import QTable
from astropy.time.core import Time
from stix.core import config

from stix.fits.calibration.integer_compression import decompress
from stix.spice.time_utils import scet_to_datetime
from stix.fits.products.common import _get_pixel_mask, _get_detector_mask, _get_compression_scheme
from stix.fits.products.common import get_energies_from_mask, get_energy_channels,get_energies_from_edges

from stix.core.logger import get_logger
logger = get_logger(__name__)




class Control(QTable):
    def __repr__(self):
        return f'<{self.__class__.__name__} \n {super().__repr__()}>'

    def _get_time(self):
        # Replicate packet time for each sample
        base_times = Time(list(chain(
            *[[scet_to_datetime(f'{self["scet_coarse"][i]}:{self["scet_fine"][i]}')]
              * n for i, n in enumerate(self['num_samples'])])))
        # For each sample generate sample number and multiply by duration and apply unit
        start_delta = np.hstack(
            [(np.arange(ns) * it) for ns, it in self[['num_samples', 'integration_time']]])
        # hstack op loses unit
        start_delta = start_delta.value * self['integration_time'].unit

        duration = np.hstack([np.ones(num_sample) * int_time for num_sample, int_time in
                              self[['num_samples', 'integration_time']]])
        duration = duration.value * self['integration_time'].unit

        # TODO Write out and simplify
        end_delta = start_delta + duration

        # Add the delta time to base times and convert to relative from start time
        times = base_times + start_delta + (end_delta - start_delta) / 2
        # times -= times[0]
        return times, duration

    @classmethod
    def from_packets(cls, packets):

        control = cls()

        # Control
        control['tc_packet_id_ref'] = np.array(packets.get('NIX00001'), np.int32)
        control['tc_packet_seq_control'] = np.array(packets.get('NIX00002'), np.int32)
        control['request_id'] = np.array(packets.get('NIX00037'), np.uint32)
        control['compression_scheme_counts_skm'] = _get_compression_scheme(packets, 'NIXD0007',
                                                                           'NIXD0008', 'NIXD0009')
        control['compression_scheme_triggers_skm'] = _get_compression_scheme(packets, 'NIXD0010',
                                                                             'NIXD0011', 'NIXD0012')
        control['time_stamp'] = np.array(packets.get('NIX00402'))

        control['num_structures'] = np.array(packets.get('NIX00403'), np.int32)

        return control

def qadd(table, key, value):
    """append key to qtable
    """
    try:
        table[key]=value
    except ValueError:
        for k in table.keys():
            logger.info(f'{k}, {table[k].shape}')
        logger.error(f'{value}')
        raise ValueError(f"Failed to add {key}, shape {value.shape}, \ntable dtype: {table.dtype}")



class Data(QTable):
    def __repr__(self):
        return f'<{self.__class__.__name__} \n {super().__repr__()}>'

    @classmethod
    def from_packets(cls, packets):
        pass


class Product(object):
    def __init__(self, control, data):
        """
        Generic product composed of control and data

        Parameters
        ----------
        control : parser.products.quicklook.Control
            Table containing control information
        data : parser.products.quicklook.Data
            Table containing data
        """
        self.type = 'sci'
        self.control = control
        self.data = data

        self.obs_beg = self.data['time'][0] - self.data['timedel'][0] / 2
        self.obs_end = self.data['time'][-1] + self.data['timedel'][-1] / 2
        self.obs_avg = self.obs_beg + (self.obs_end - self.obs_beg) / 2

    def __add__(self, other):
        other.control['index'] = other.control['index'] + self.control['index'].max() + 1
        control = vstack((self.control, other.control))
        cnames = control.colnames
        cnames.remove('index')
        control = unique(control, cnames)

        other.data['control_index'] = other.data['control_index'] + self.control['index'].max() + 1
        data = vstack((self.data, other.data))

        data_ind = np.isin(data['control_index'], control['index'])
        data = data[data_ind]

        return type(self)(control, data)

    def __repr__(self):
        return f'<{self.__class__.__name__}\n' \
               f' {self.control.__repr__()}\n' \
               f' {self.data.__repr__()}\n' \
               f'>'

    def get_energies(self):
        obs_time = self.obs_avg
        if 'e_unit' in self.control.colnames:
            #L4 request
            energies=get_energies_from_edges(obs_time=obs_time,
                    e_low=self.control['e_low'][0],
                    e_high=self.control['e_high'][0], e_unit=self.control['e_unit'][0])

        elif 'e_low' in self.control.colnames and 'e_high' in self.control.colnames:
            energies=get_energies_from_edges(obs_time=obs_time,
                    e_low=self.control['e_low'][0],
                    e_high=self.control['e_high'][0], e_unit=None)

        elif 'energy_bin_edge_mask' in self.control.colnames:
            energies = get_energies_from_mask(obs_time=obs_time , 
                    mask=self.control['energy_bin_edge_mask'][0])
        elif 'energy_bin_mask' in self.control.colnames:
            energies = get_energies_from_mask(obs_time=obs_time, 
                    mask=self.control['energy_bin_mask'][0])
        else:
            energies = get_energies_from_mask(obs_time=obs_time)

        return energies

    @staticmethod
    def get_energy_channel_dict(start_time):
        return get_energy_channels(start_time)

    @classmethod
    def from_fits(cls, fitspath):
        header = fits.getheader(fitspath)
        control = QTable.read(fitspath, hdu='CONTROL')
        data = QTable.read(fitspath, hdu='DATA')
        obs_beg = Time(header['DATE_OBS'])
        data['time'] = (data['time'] + obs_beg)
        return cls(control=control, data=data)

    def to_requests(self):
        for ci in unique(self.control, keys=['tc_packet_seq_control', 'request_id'])['index']:
            control = self.control[self.control['index'] == ci]
            data = self.data[self.data['control_index'] == ci]
        # for req_id in self.control['request_id']:
        #     ctrl_inds = np.where(self.control['request_id'] == req_id)
        #     control = self.control[ctrl_inds]
        #     data_index = control['index'][0]
        #     data_inds = np.where(self.data['control_index'] == data_index)
        #     data = self.data[data_inds]

            yield type(self)(control=control, data=data)




class XrayL0(Product):
    """
    X-ray Level 0 data
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'xray-rpd'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        control = Control.from_packets(packets)

        control.remove_column('num_structures')
        control = unique(control)

        if len(control) != 1:
            raise ValueError('Creating a science product form packets from multiple products')

        control['index'] = 0

        data = Data()
        data['start_time'] = (np.array(packets.get('NIX00404'), np.uint16)) * 0.1 * u.s
        try:
            data['rcr'] = np.array(packets.get('NIX00401'), np.ubyte)
        except Exception as e:
            data['rcr'] = np.array(packets.get('NIX00401')[0], np.ubyte)
        data['integration_time'] = (np.array(packets.get('NIX00405')[0], np.int16)) * 0.1 * u.s
        data['pixel_masks'] = _get_pixel_mask(packets, 'NIXD0407')
        data['detector_masks'] = _get_detector_mask(packets)
        data['triggers'] = np.array([packets.get(f'NIX00{i}') for i in range(408, 424)],
                                    np.int64).T
        data['num_samples'] = np.array(packets.get('NIX00406'), np.int16)

        num_detectors = 32
        num_energies = 32
        num_pixels = 12

        # Data
        tmp = dict()
        tmp['pixel_id'] = np.array(packets.get('NIXD0158'), np.ubyte)
        tmp['detector_id'] = np.array(packets.get('NIXD0153'), np.ubyte)
        tmp['channel'] = np.array(packets.get('NIXD0154'), np.ubyte)
        tmp['continuation_bits'] = packets.get('NIXD0159', np.ubyte)

        control['energy_bin_edge_mask'] = np.full((1, 32), False, np.ubyte)
        all_energies = set(tmp['channel'])
        control['energy_bin_edge_mask'][:, list(all_energies)] = True

        # Find contiguous time indices
        unique_times = np.unique(data['start_time'])
        time_indices = np.searchsorted(unique_times, data['start_time'])

        # Create full count array 0s are not send down, if cb = 0 1 count, for cb 1 just extract
        # and for cb 2 extract and sum
        raw_counts = packets.get('NIX00065')
        counts_1d = []
        raw_count_index = 0
        for cb in tmp['continuation_bits']:
            if cb == 0:
                counts_1d.append(1)
            elif cb == 1:
                cur_count = raw_counts[raw_count_index]
                counts_1d.append(cur_count)
                raw_count_index += cb
            elif cb == 2:
                cur_count = raw_counts[raw_count_index:(raw_count_index + cb)]
                combined_count = int.from_bytes((cur_count[0]+1).to_bytes(2, 'big')
                        + cur_count[1].to_bytes(1, 'big'), 'big')
                counts_1d.append(combined_count)
                raw_count_index += cb
            else:
                raise ValueError(f'Continuation bits value of {cb} not allowed (0, 1, 2)')
        counts_1d = np.array(counts_1d, np.uint16)
        # raw_counts = counts_1d

        end_inds = np.cumsum(data['num_samples'])
        start_inds = np.hstack([0, end_inds[:-1]])
        dd = [(tmp['pixel_id'][s:e], tmp['detector_id'][s:e], tmp['channel'][s:e], counts_1d[s:e])
              for s, e in zip(start_inds.astype(int), end_inds)]

        counts = np.zeros((len(unique_times), num_detectors, num_pixels, num_energies), np.uint32)
        for i, (pid, did, cid, cc) in enumerate(dd):
            counts[time_indices[i], did, pid, cid] = cc

     
        sub_index = np.searchsorted(data['start_time'], unique_times)
        data = data[sub_index]
        data['time'] = Time(scet_to_datetime(f'{int(control["time_stamp"][0])}:0'))\
            + data['start_time'] + data['integration_time']/2
        data['timedel'] = data['integration_time']
        data['counts'] = counts * u.ct
        data['control_index'] = control['index'][0]

        data.remove_columns(['start_time', 'integration_time', 'num_samples'])

        return cls(control=control, data=data)


    


class XrayL1(Product):
    """
    X-ray Compression Level 1/2 data
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'xray-cpd'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        # Control

        ssid = packets['SSID'][0]

        control = Control.from_packets(packets)

        control.remove_column('num_structures')
        control = unique(control)

        if len(control) != 1:
            #print(packets[0])
            #raise ValueError('Creating a science product form packets from multiple products')
            print('Control is not unique')
            

        control['index'] = 0

        data = Data()
        try:
            data['delta_time'] = (np.array(packets['NIX00441'], np.int32)) * 0.1 * u.s
        except KeyError:
            #replaced with NIX00404 for versions after asw v180, 
            data['delta_time'] = (np.array(packets['NIX00404'], np.int32)) * 0.1 * u.s

        unique_times = np.unique(data['delta_time'])


        qadd(data, 'rcr',  np.array(packets['NIX00401'], np.ubyte))

        pixel_sets= np.array(packets['NIX00442'], np.ubyte)

        qadd(data,'num_pixel_sets',  pixel_sets)

        pixel_masks = _get_pixel_mask(packets, 'NIXD0407')
        pixel_masks = pixel_masks.reshape(-1, pixel_sets[0], 12)
        if ssid == 21 and pixel_sets[0] != 12:
            pixel_masks = np.pad(pixel_masks, ((0, 0),  (0, 12- pixel_sets[0]), (0, 0)))
        qadd(data,'pixel_masks', pixel_masks)
        qadd(data,'detector_masks', _get_detector_mask(packets))
        qadd(data,'integration_time', (np.array(packets.get('NIX00405'), np.uint16)) * 0.1 * u.s)

        # TODO change once FSW fixed
        #ts, tk, tm = control['compression_scheme_counts_skm'][0]
        
        ts, tk, tm = control['compression_scheme_triggers_skm'][0]

        trigger_from_packets=[packets.get(f'NIX00{i}') for i in range(242, 258)]


        triggers, triggers_var = decompress(trigger_from_packets,
                                            s=ts, k=tk, m=tm,
                                            return_variance=True)

        qadd(data,'triggers', triggers.T)
        qadd(data,'triggers_comp_err', np.sqrt(triggers_var).T)
        qadd(data,'num_energy_groups', np.array(packets['NIX00258'], np.ubyte))

        tmp = dict()
        tmp['e_low'] = np.array(packets['NIXD0016'], np.ubyte)
        tmp['e_high'] = np.array(packets['NIXD0017'], np.ubyte)
        tmp['num_data_elements'] = np.array(packets['NIX00259'])
        unique_energies_low = np.unique(tmp['e_low'])
        unique_energies_high = np.unique(tmp['e_high']) 

        control['e_low']=[unique_energies_low]
        control['e_high']=[unique_energies_high]


        # counts = np.array(eng_packets['NIX00260'], np.uint32)

        cs, ck, cm = control['compression_scheme_counts_skm'][0]
        counts, counts_var = decompress(packets.get('NIX00260'), s=cs, k=ck, m=cm,
                                        return_variance=True)

        #print('unique id:', packets['NIX00037'])
        #print('unique times',len(unique_times), unique_times)
        #print('pixel mask',len(data['pixel_masks']), pixel_masks.size)

        #print(unique_times.size, 
        #                        data['detector_masks'][0].sum(), data['num_pixel_sets'][0].sum(),unique_energies_low.size)
        counts = counts.reshape(unique_times.size, unique_energies_low.size,
                                data['detector_masks'][0].sum(), data['num_pixel_sets'][0].sum())
        #comment from Hualin, 2021, Sept. probably there is a bug here, when the number of energy bins is not 32, it cashes
        #maybe need to replaced to:
        #print(unique_times.size, 
        #                        data['detector_masks'][0].sum(), data['num_pixel_sets'][0].sum(),unique_energies_low.size)

        counts_var = counts_var.reshape(unique_times.size, unique_energies_low.size,
                                        data['detector_masks'][0].sum(),
                                        data['num_pixel_sets'][0].sum())
        # t x e x d x p -> t x d x p x e
        counts = counts.transpose((0, 2, 3, 1))
        counts_var = np.sqrt(counts_var.transpose((0, 2, 3, 1)))
        if ssid == 21:
            out_counts = np.zeros((unique_times.size, 32, 12, 32))
            out_var = np.zeros((unique_times.size, 32, 12, 32))
        elif ssid == 22:
            out_counts = np.zeros((unique_times.size, 32, 4, 32))
            out_var = np.zeros((unique_times.size, 32, 4, 32))


        data['time'] = Time(scet_to_datetime(f'{int(control["time_stamp"][0])}:0')) \
            + data['delta_time'] + data['integration_time'] / 2

        energy_channels=Product.get_energy_channel_dict(data['time'][0])

        dl_energies = np.array([[energy_channels[lch]['e_lower'], energy_channels[hch]['e_upper']]
            for lch, hch in zip(unique_energies_low, unique_energies_high)]).reshape(-1)
        dl_energies = np.unique(dl_energies)
        
        sci_energies = np.hstack([[energy_channels[ch]['e_lower'] for ch in range(32)],
                                  energy_channels[31]['e_upper']])

        # If there is any onboard summing of energy channels rebin back to standard sci channels
        #print(config.ASW_VERSION)
        logger.info("rebining energies..")
        #print(unique_energies_low)
        #print(unique_energies_high)
        if (unique_energies_high - unique_energies_low).sum() > 0:
            # there is a bug here
            #print('Onboard summing rebinned ')
            rebinned_counts = np.zeros((*counts.shape[:-1], 32))
            rebinned_counts_var = np.zeros((*counts_var.shape[:-1], 32))
            e_ch_start = 0
            e_ch_end = counts.shape[-1]
            if dl_energies[0] == 0.0:
                rebinned_counts[..., 0] = counts[..., 0]
                rebinned_counts_var[..., 0] = counts_var[..., 0]
                e_ch_start += 1
            elif dl_energies[-1] == np.inf:
                rebinned_counts[..., -1] = counts[..., -1]
                rebinned_counts_var[..., -1] = counts_var[..., -1]
                e_ch_end -= 1

            torebin = np.where((dl_energies >= 4.0) & (dl_energies <= 150.0))
            rebinned_counts[..., 1:-1] = np.apply_along_axis(rebin_proportional, -1,
                    counts[..., e_ch_start:e_ch_end].reshape(-1, e_ch_end-e_ch_start),
                    dl_energies[torebin],
                    sci_energies[1:-1]).reshape((*counts.shape[:-1], 30))

            rebinned_counts_var[..., 1:-1] = np.apply_along_axis(rebin_proportional, -1,
                    counts_var[..., e_ch_start:e_ch_end].reshape(-1, e_ch_end-e_ch_start),
                    dl_energies[torebin],
                    sci_energies[1:-1]).reshape((*counts_var.shape[:-1], 30))

            energy_indices = np.full(32, True)
            energy_indices[[0,-1]] = False

            ix = np.ix_(np.full(unique_times.size, True), data['detector_masks'][0].astype(bool),
                        np.ones(data['num_pixel_sets'][0], dtype=bool), np.full(32, True))

            out_counts[ix] = rebinned_counts
            out_var[ix] = rebinned_counts_var
        else:
            energy_indices = np.full(32, False)
            energy_indices[unique_energies_low.min():unique_energies_high.max()+1] = True

            ix = np.ix_(np.full(unique_times.size, True),
                              data['detector_masks'][0].astype(bool),
                              np.ones(data['num_pixel_sets'][0], dtype=bool),
                              energy_indices)

            out_counts[ix] = counts
            out_var[ix] = counts_var


        if counts.sum() != out_counts.sum():
            #import ipdb; ipdb.set_trace()
            raise ValueError(f'Original and reformatted count totals do not match: {counts.sum()} vs {out_counts.sum()}')

        control['energy_bin_edge_mask'] = np.full((1, 32), False, np.ubyte)

        #control['energy_science_bins']=np.hstack([unique_energies_low, unique_energies_high])
        #print("Energy bins:", control['energy_science_bins'])
        all_sci_energies = np.hstack([tmp['e_low'], tmp['e_high']])
        all_energies = set(all_sci_energies)
        control['energy_bin_edge_mask'][:, list(all_energies)] = True
        # time x energy x detector x pixel
        # counts = np.array(
        #     eng_packets['NIX00260'], np.uint16).reshape(unique_times.size, num_energies,
        #                                                 num_detectors, num_pixels)
        # time x channel x detector x pixel need to transpose to time x detector x pixel x channel

        sub_index = np.searchsorted(data['delta_time'], unique_times)
        data = data[sub_index]


        logger.info("adding data to qtable..")
        qadd(data,'timedel', data['integration_time'])
        qadd(data,'counts', out_counts * u.ct)
        qadd(data,'counts_comp_err', out_var * u.ct)
        qadd(data,'control_index', control['index'][0])

        data.remove_columns(['delta_time', 'integration_time'])

        data = data['time', 'timedel', 'rcr', 'pixel_masks', 'detector_masks', 'num_pixel_sets',
                    'num_energy_groups', 'triggers', 'triggers_comp_err', 'counts', 'counts_comp_err']
        data['control_index'] = 0
        logger.info("qtable ready")

        return cls(control=control, data=data)


class XrayL2(XrayL1):
    """
    X-ray Compression Level 2 data
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'xray-scpd'
        self.level = 'L1A'


class XrayL3(Product):
    """
    X-ray Compression Level 3 data (visibilities)
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'xray-vis'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        # Control
        control = Control.from_packets(packets)
        control.remove_column('num_structures')
        control = unique(control)
        if len(control) != 1:
            raise ValueError()
        control['index'] = range(len(control))

        data = Data()
        try:
            data['control_index'] = np.full(len(packets['NIX00441']), 0)
            data['delta_time'] = (np.array(packets['NIX00441'], np.uint16)) * 0.1 * u.s
        except KeyError:
            data['control_index'] = np.full(len(packets['NIX00404']), 0)
            data['delta_time'] = (np.array(packets['NIX00404'], np.int32)) * 0.1 * u.s


        unique_times = np.unique(data['delta_time'])

        # time = np.array([])
        # for dt in set(self.delta_time):
        #     i, = np.where(self.delta_time == dt)
        #     nt = sum(np.array(packets['NIX00258'])[i])
        #     time = np.append(time, np.repeat(dt, nt))
        # self.time = time

        data['rcr'] = packets['NIX00401']
        data['pixel_mask1'] = _get_pixel_mask(packets, 'NIXD0407')
        data['pixel_mask2'] = _get_pixel_mask(packets, 'NIXD0444')
        data['pixel_mask3'] = _get_pixel_mask(packets, 'NIXD0445')
        data['pixel_mask4'] = _get_pixel_mask(packets, 'NIXD0446')
        data['pixel_mask5'] = _get_pixel_mask(packets, 'NIXD0447')
        data['detector_masks'] = _get_detector_mask(packets)
        data['integration_time'] = (np.array(packets['NIX00405'])) * 0.1

        ts, tk, tm = control['compression_scheme_triggers_skm'][0]
        triggers, triggers_var = decompress([packets[f'NIX00{i}'] for i in range(242, 258)],
                                            s=ts, k=tk, m=tm, return_variance=True)

        data['triggers'] = triggers.T
        data['triggers_comp_err'] = np.sqrt(triggers_var).T

        tids = np.searchsorted(data['delta_time'], unique_times)
        data = data[tids]

        num_energy_groups = sum(packets['NIX00258'])

        # Data
        vis = np.zeros((unique_times.size, 32, 32), dtype=complex)
        vis_err = np.zeros((unique_times.size, 32, 32), dtype=complex)
        e_low = np.array(packets['NIXD0016'])
        e_high = np.array(packets['NIXD0017'])

        # TODO create energy bin mask
        control['energy_bin_edge_mask'] = np.full((1, 32), False, np.ubyte)
        all_energies = set(np.hstack([e_low,e_high]))
        control['energy_bin_edge_mask'][:, list(all_energies)] = True

        data['flux'] = np.array(packets['NIX00261']).reshape(unique_times.size, -1)
        num_detectors = packets['NIX00262'][0]
        detector_id = np.array(packets['NIX00100']).reshape(unique_times.size, -1,
                                                            num_detectors)

        # vis[:, detector_id[0], e_low.reshape(unique_times.size, -1)[0]] = (
        #         np.array(packets['NIX00263']) + np.array(packets['NIX00264'])
        #         * 1j).reshape(unique_times.size, num_detectors, -1)

        ds, dk, dm = control['compression_scheme_counts_skm'][0]
        real, real_var = decompress(packets['NIX00263'], s=ds, k=dk, m=dm, return_variance=True)
        imaginary, imaginary_var = decompress(packets['NIX00264'], s=ds, k=dk, m=dm,
                                              return_variance=True)

        mesh = np.ix_(np.arange(unique_times.size), detector_id[0][0],
                      e_low.reshape(unique_times.size, -1)[0])
        vis[mesh] = (real + imaginary * 1j).reshape(unique_times.size, num_detectors, -1)

        # TODO this doesn't seem correct prob need combine in a better
        vis_err[mesh] = (np.sqrt(real_var) + np.sqrt(imaginary_var)*1j).reshape(unique_times.size,
                                                                                num_detectors, -1)

        data['visibility'] = vis
        data['visibility_err'] = vis_err

        data['time'] = Time(scet_to_datetime(f'{int(control["time_stamp"][0])}:0')) \
            + data['delta_time'] + data['integration_time'] / 2
        data['timedel'] = data['integration_time']

        return cls(control=control, data=data)


class Spectrogram(Product):
    """
    Spectrogram data
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'xray-spec'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        # Control
        control = Control.from_packets(packets)

        control['pixel_mask'] = np.unique(_get_pixel_mask(packets), axis=0)
        control['detector_mask'] = np.unique(_get_detector_mask(packets), axis=0)

        raw_rcr=np.array(packets.get('NIX00401'), np.ubyte)
        #print(rcr)
        num_times = np.array(packets.get('NIX00089'))
        #    #raise an exception for  data  acquired on Oct 28, 2021
        #except ValueError:
        #print(rcr)
        #print(control['pixel_mask'])
        rcr = np.hstack([np.full(nt, rcr) for rcr, nt in zip(raw_rcr, num_times)]).astype(np.ubyte)

        #control['rcr'] =rcr
        control['index'] = range(len(control))

        e_min = np.array(packets['NIXD0442'])
        e_max = np.array(packets['NIXD0443'])
        energy_unit = np.array(packets['NIXD0019']) + 1


        num_times = np.array(packets['NIX00089'])
        total_num_times = num_times.sum()

        cs, ck, cm = control['compression_scheme_counts_skm'][0]

        counts, counts_var = decompress(packets['NIX00268'], s=cs, k=ck, m=cm, return_variance=True)
        counts = counts.reshape(total_num_times, -1)
        counts_var = counts_var.reshape(total_num_times, -1)
        full_counts = np.zeros((total_num_times, 32))
        full_counts_var = np.zeros((total_num_times, 32))

        cids = [np.arange(emin, emax+1, eunit) for (emin, emax, eunit)
                in zip(e_min, e_max, energy_unit)]

        control['energy_bin_edge_mask'] = np.full((1, 32), False, np.ubyte)
        control['energy_bin_edge_mask'][:, cids] = True
        control['e_low']=[e_min[0]]
        control['e_high']=[e_max[0]]
        control['e_unit']=[energy_unit[0]]
        #print("Emin:",e_min, e_max, energy_unit)
        

        start=scet_to_datetime(f'{int(control["time_stamp"][0])}:0')
        energy_channels=Product.get_energy_channel_dict(start)

        dl_energies = np.array([[energy_channels[ch]['e_lower'] for ch in chs]
                               + [energy_channels[chs[-1]]['e_upper']] for chs in cids][0])

        sci_energies = np.hstack([[energy_channels[ch]['e_lower'] for ch in range(32)],
                                  energy_channels[31]['e_upper']])
        ind = 0
        for nt in num_times:
            e_ch_start = 0
            e_ch_end = counts.shape[1]
            if dl_energies[0] == 0:
                full_counts[ind:ind+nt, 0] = counts[ind:ind+nt, 0]
                full_counts_var[ind:ind+nt, 0] = counts_var[ind:ind+nt, 0]
                e_ch_start = 1
            if dl_energies[-1] == np.inf:
                full_counts[ind:ind+nt, -1] = counts[ind:ind+nt, -1]
                full_counts_var[ind:ind+nt, -1] = counts[ind:ind+nt, -1]
                e_ch_end -= 1

            torebin = np.where((dl_energies >= 4.0) & (dl_energies <= 150.0))
            full_counts[ind:ind+nt, 1:-1] = np.apply_along_axis(rebin_proportional, 1,
                                                                counts[ind:ind+nt,
                                                                    e_ch_start:e_ch_end],
                                                                dl_energies[torebin],
                                                                sci_energies[1:-1])

            full_counts_var[ind:ind+nt, 1:-1] = np.apply_along_axis(rebin_proportional, 1,
                                                                    counts_var[ind:ind+nt,
                                                                        e_ch_start:e_ch_end],
                                                                    dl_energies[torebin],
                                                                    sci_energies[1:-1])

            ind += nt

        if counts.sum() != full_counts.sum():
            raise ValueError('Original and reformatted count totals do not match')

        try:
            delta_time = (np.array(packets['NIX00441'], np.uint32)) * 0.1 * u.s
        except KeyError:
            delta_time = (np.array(packets['NIX00404'], np.uint32)) * 0.1 * u.s



        closing_time_offset = (np.array(packets['NIX00269'], np.uint32)) * 0.1 * u.s

        # TODO incorporate into main loop above
        centers = []
        deltas = []
        last = 0
        for i, nt in enumerate(num_times):
            edge = np.hstack(
                [delta_time[last:last + nt], delta_time[last + nt - 1] + closing_time_offset[i]])
            delta = np.diff(edge)
            center = edge[:-1] + delta / 2
            centers.append(center)
            deltas.append(delta)
            last = last + nt

        centers = np.hstack(centers)
        deltas = np.hstack(deltas)

        # Data
        #print("RCR shape:", rcr.shape)
        #print("detla shape:", deltas.shape)


        data = Data()
        data['time'] = Time(scet_to_datetime(f'{int(control["time_stamp"][0])}:0')) \
            + centers
        data['timedel'] = deltas
        data['rcr']=rcr



        ts, tk, tm = control['compression_scheme_triggers_skm'][0]
        triggers, triggers_var = decompress(packets['NIX00267'], s=ts, k=tk, m=tm,
                                            return_variance=True)

        data['triggers'] = triggers
        data['triggers_comp_err'] = np.sqrt(triggers_var)
        data['counts'] = full_counts * u.ct
        data['counts_comp_err'] = np.sqrt(full_counts_var) * u.ct
        data['control_index'] = 0

        return cls(control=control, data=data)


class Aspect(Product):
    """
    Aspect
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'burst-aspect'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        # Header
        control = Control()
        scet_coarse = packets['NIX00445']
        scet_fine = packets['NIX00446']
        start_times = Time([scet_to_datetime(f'{scet_coarse[i]}:{scet_fine[i]}')
                            for i in range(len(scet_coarse))])

        control['summing_value'] = packets['NIX00088']
        control['averaging_value'] = packets['NIX00490']
        control['index'] = range(len(control))
        if 'NIX00037' in packets:
            control['request_id'] = np.array(packets['NIX00037'], np.uint32)
            #request id introduced after FSW v183, 2021-12-09
        #except KeyError:
        #    pass

        delta_time = ((control['summing_value'] * control['averaging_value']) / 1000.0)
        samples = packets['NIX00089']

        offsets = [delta_time[i]* np.arange(ns) * u.s for i, ns in enumerate(samples)]
        time_st=[start_times[i] + offsets[i] for i in range(len(offsets))]
        try: 
            time = Time(time_st)
        except AttributeError:
            #throw exception after asw 183, updated by Hualin in order to 
            time=np.array(time_st).flatten()

        timedel = np.hstack(offsets)

        # Data
        try:
            data = Data()
            data['time'] = time
            data['timedel'] = timedel
            data['cha_diode0'] = packets['NIX00090']
            data['cha_diode1'] = packets['NIX00091']
            data['chb_diode0'] = packets['NIX00092']
            data['chb_diode1'] = packets['NIX00093']
            data['control_index'] = np.hstack([np.full(ns, i) for i, ns in enumerate(samples)])
        except ValueError as e:
            logger.warning(e)
            #raise
            return None

        return cls(control=control, data=data)

    def to_requests(self):
        return [type(self)(control=self.control, data=self.data)]
        #updated by Hualin Xiao, the lines below created multiple fits as reported by Frederic
        """
        days = set([(t.year, t.month, t.day) for t in self.data['time'].to_datetime()])
        date_ranges = [(datetime(*day), datetime(*day) + timedelta(days=1)) for day in days]
        for dstart, dend in date_ranges:
            i = np.where((self.data['time'] >= dstart) &
                         (self.data['time'] < dend))

            data = self.data[i]
            control_indices = np.unique(data['control_index'])
            control = self.control[np.isin(self.control['index'], control_indices)]
            control_index_min = control_indices.min()

            data['control_index'] = data['control_index'] - control_index_min
            control['index'] = control['index'] - control_index_min
        type(self)(control=control, data=data)
        """


def rebin_proportional(y1, x1, x2):
    x1 = np.asarray(x1)
    y1 = np.asarray(y1)
    x2 = np.asarray(x2)

    # the fractional bin locations of the new bins in the old bins
    i_place = np.interp(x2, x1, np.arange(len(x1)))

    cum_sum = np.r_[[0], np.cumsum(y1)]

    # calculate bins where lower and upper bin edges span
    # greater than or equal to one original bin.
    # This is the contribution from the 'intact' bins (not including the
    # fractional start and end parts.
    whole_bins = np.floor(i_place[1:]) - np.ceil(i_place[:-1]) >= 1.
    start = cum_sum[np.ceil(i_place[:-1]).astype(int)]
    finish = cum_sum[np.floor(i_place[1:]).astype(int)]

    y2 = np.where(whole_bins, finish - start, 0.)

    bin_loc = np.clip(np.floor(i_place).astype(int), 0, len(y1) - 1)

    # fractional contribution for bins where the new bin edges are in the same
    # original bin.
    same_cell = np.floor(i_place[1:]) == np.floor(i_place[:-1])
    frac = i_place[1:] - i_place[:-1]
    contrib = (frac * y1[bin_loc[:-1]])
    y2 += np.where(same_cell, contrib, 0.)

    # fractional contribution for bins where the left and right bin edges are in
    # different original bins.
    different_cell = np.floor(i_place[1:]) > np.floor(i_place[:-1])
    frac_left = np.ceil(i_place[:-1]) - i_place[:-1]
    contrib = (frac_left * y1[bin_loc[:-1]])

    frac_right = i_place[1:] - np.floor(i_place[1:])
    contrib += (frac_right * y1[bin_loc[1:]])

    y2 += np.where(different_cell, contrib, 0.)

    return y2
