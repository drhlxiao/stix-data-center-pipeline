"""
High level STIX data products created from single stand alone packets or a sequence of packects.
"""
from datetime import timedelta, datetime
from itertools import chain

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import QTable, vstack, unique
from astropy.time import Time

from stix.fits.calibration.integer_compression import decompress
from stix.spice.time_utils import scet_to_datetime
from stix.fits.products.common import _get_detector_mask, _get_pixel_mask, _get_energy_bins, \
    _get_compression_scheme, _get_num_energies, _get_sub_spectrum_mask

# __all__ = ['LightCurve', 'Background', 'Spectra', 'Variance', 'CalibrationSpectra',
#            'FlareFlagAndLocation', 'TMManagementAndFlareList', 'get_energies_from_mask']

ENERGY_CHANNELS = {
    0: {'channel_edge': 0, 'energy_edge': 0, 'e_lower': 0.0, 'e_upper': 4.0, 'bin_width': 4.0,
        'dE_E': 2.000, 'ql_channel': None},
    1: {'channel_edge': 1, 'energy_edge': 4, 'e_lower': 4.0, 'e_upper': 5.0, 'bin_width': 1.0,
        'dE_E': 0.222, 'ql_channel': 0},
    2: {'channel_edge': 2, 'energy_edge': 5, 'e_lower': 5.0, 'e_upper': 6.0, 'bin_width': 1.0,
        'dE_E': 0.182, 'ql_channel': 0},
    3: {'channel_edge': 3, 'energy_edge': 6, 'e_lower': 6.0, 'e_upper': 7.0, 'bin_width': 1.0,
        'dE_E': 0.154, 'ql_channel': 0},
    4: {'channel_edge': 4, 'energy_edge': 7, 'e_lower': 7.0, 'e_upper': 8.0, 'bin_width': 1.0,
        'dE_E': 0.133, 'ql_channel': 0},
    5: {'channel_edge': 5, 'energy_edge': 8, 'e_lower': 8.0, 'e_upper': 9.0, 'bin_width': 1.0,
        'dE_E': 0.118, 'ql_channel': 0},
    6: {'channel_edge': 6, 'energy_edge': 9, 'e_lower': 9.0, 'e_upper': 10.0, 'bin_width': 1.0,
        'dE_E': 0.105, 'ql_channel': 0},
    7: {'channel_edge': 7, 'energy_edge': 10, 'e_lower': 10.0, 'e_upper': 11.0, 'bin_width': 1.0,
        'dE_E': 0.095, 'ql_channel': 1},
    8: {'channel_edge': 8, 'energy_edge': 11, 'e_lower': 11.0, 'e_upper': 12.0, 'bin_width': 1.0,
        'dE_E': 0.087, 'ql_channel': 1},
    9: {'channel_edge': 9, 'energy_edge': 12, 'e_lower': 12.0, 'e_upper': 13.0, 'bin_width': 1.0,
        'dE_E': 0.080, 'ql_channel': 1},
    10: {'channel_edge': 10, 'energy_edge': 13, 'e_lower': 13.0, 'e_upper': 14.0, 'bin_width': 1.0,
         'dE_E': 0.074, 'ql_channel': 1},
    11: {'channel_edge': 11, 'energy_edge': 14, 'e_lower': 14.0, 'e_upper': 15.0, 'bin_width': 1.0,
         'dE_E': 0.069, 'ql_channel': 1},
    12: {'channel_edge': 12, 'energy_edge': 15, 'e_lower': 15.0, 'e_upper': 16.0, 'bin_width': 1.0,
         'dE_E': 0.065, 'ql_channel': 2},
    13: {'channel_edge': 13, 'energy_edge': 16, 'e_lower': 16.0, 'e_upper': 18.0, 'bin_width': 1.0,
         'dE_E': 0.061, 'ql_channel': 2},
    14: {'channel_edge': 14, 'energy_edge': 18, 'e_lower': 18.0, 'e_upper': 20.0, 'bin_width': 2.0,
         'dE_E': 0.105, 'ql_channel': 2},
    15: {'channel_edge': 15, 'energy_edge': 20, 'e_lower': 20.0, 'e_upper': 22.0, 'bin_width': 2.0,
         'dE_E': 0.095, 'ql_channel': 2},
    16: {'channel_edge': 16, 'energy_edge': 22, 'e_lower': 22.0, 'e_upper': 25.0, 'bin_width': 3.0,
         'dE_E': 0.128, 'ql_channel': 2},
    17: {'channel_edge': 17, 'energy_edge': 25, 'e_lower': 25.0, 'e_upper': 28.0, 'bin_width': 3.0,
         'dE_E': 0.113, 'ql_channel': 3},
    18: {'channel_edge': 18, 'energy_edge': 28, 'e_lower': 28.0, 'e_upper': 32.0, 'bin_width': 4.0,
         'dE_E': 0.133, 'ql_channel': 3},
    19: {'channel_edge': 19, 'energy_edge': 32, 'e_lower': 32.0, 'e_upper': 36.0, 'bin_width': 4.0,
         'dE_E': 0.118, 'ql_channel': 3},
    20: {'channel_edge': 20, 'energy_edge': 36, 'e_lower': 36.0, 'e_upper': 40.0, 'bin_width': 4.0,
         'dE_E': 0.105, 'ql_channel': 3},
    21: {'channel_edge': 21, 'energy_edge': 40, 'e_lower': 40.0, 'e_upper': 45.0, 'bin_width': 5.0,
         'dE_E': 0.118, 'ql_channel': 3},
    22: {'channel_edge': 22, 'energy_edge': 45, 'e_lower': 45.0, 'e_upper': 50.0, 'bin_width': 5.0,
         'dE_E': 0.105, 'ql_channel': 3},
    23: {'channel_edge': 23, 'energy_edge': 50, 'e_lower': 50.0, 'e_upper': 56.0, 'bin_width': 6.0,
         'dE_E': 0.113, 'ql_channel': 4},
    24: {'channel_edge': 24, 'energy_edge': 56, 'e_lower': 56.0, 'e_upper': 63.0, 'bin_width': 7.0,
         'dE_E': 0.118, 'ql_channel': 4},
    25: {'channel_edge': 25, 'energy_edge': 63, 'e_lower': 63.0, 'e_upper': 70.0, 'bin_width': 7.0,
         'dE_E': 0.105, 'ql_channel': 4},
    26: {'channel_edge': 26, 'energy_edge': 70, 'e_lower': 70.0, 'e_upper': 76.0, 'bin_width': 6.0,
         'dE_E': 0.082, 'ql_channel': 4},
    27: {'channel_edge': 27, 'energy_edge': 76, 'e_lower': 76.0, 'e_upper': 84.0, 'bin_width': 8.0,
         'dE_E': 0.100, 'ql_channel': 4},
    28: {'channel_edge': 28, 'energy_edge': 84, 'e_lower': 84.0, 'e_upper': 100.0,
         'bin_width': 16.0, 'dE_El': 0.174, 'ql_channel': 4},
    29: {'channel_edge': 29, 'energy_edge': 100, 'e_lower': 100.0, 'e_upper': 120.0,
         'bin_width': 20.0, 'dE_El': 0.182, 'ql_channel': 4},
    30: {'channel_edge': 30, 'energy_edge': 120, 'e_lower': 120.0, 'e_upper': 150.0,
         'bin_width': 30.0, 'dE_El': 0.222, 'ql_channel': 4},
    31: {'channel_edge': 31, 'energy_edge': 150, 'e_lower': 150.0, 'e_upper': np.inf,
         'bin_width': np.inf, 'dE_E': np.inf, 'ql_channel': None}
}


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
        # Header
        control = cls()
        # self.energy_bin_mask = None
        # self.samples = None
        control['scet_coarse'] = np.array(packets['NIX00445'], np.uint32)
        # Not all QL data have fine time in TM default to 0 if no present
        scet_fine = packets.get('NIX00446')
        if scet_fine:
            control['scet_fine'] = np.array(scet_fine, np.uint32)
        else:
            control['scet_fine'] = np.zeros_like(control['scet_coarse'], np.uint32)

        integration_time = packets.get('NIX00405')
        if integration_time:
            control['integration_time'] = (np.array(integration_time, np.float) + 1) * 0.1 * u.s
        else:
            control['integration_time'] = np.zeros_like(control['scet_coarse'], np.float) * u.s

        # control = unique(control)
        control['index'] = np.arange(len(control))

        return control


class Data(QTable):
    def __repr__(self):
        return f'<{self.__class__.__name__} \n {super().__repr__()}>'

    @classmethod
    def from_packets(cls, packets):
        raise NotImplementedError


class Product:
    def __init__(self, control, data):
        """
        Generic product compose of control and data
        Parameters
        ----------
        control : parser.products.quicklook.Control
            Table containing control information
        data : parser.products.quicklook.Data
            Table containing data
        """
        self.type = 'ql'
        self.control = control
        self.data = data

        self.obs_beg = self.data['time'][0] - self.control['integration_time'][0] / 2
        self.obs_end = self.data['time'][-1] + self.control['integration_time'][-1] / 2
        self.obs_avg = self.obs_beg + (self.obs_end - self.obs_beg) / 2

    def __add__(self, other):
        """
        Combine two products stacking data along columns and removing duplicated data using time as
        the primary key.

        Parameters
        ----------
        other : A subclass of parser.products.quicklook.Product

        Returns
        -------
        A subclass of parser.products.quicklook.Product
            The combined data product
        """
        if not isinstance(other, type(self)):
            raise TypeError(f'Products must of same type not {type(self)} and {type(other)}')

        # TODO reindex and update data control_index
        other.control['index'] = other.control['index'] + self.control['index'].max() + 1
        control = vstack((self.control, other.control))
        # control = unique(control, keys=['scet_coarse', 'scet_fine'])
        # control = control.group_by(['scet_coarse', 'scet_fine'])

        other.data['control_index'] = other.data['control_index'] + self.control['index'].max() + 1
        data = vstack((self.data, other.data))
        data = unique(data, keys='time')
        # data = data.group_by('time')
        unique_control_inds = np.unique(data['control_index'])
        control = control[np.isin(control['index'], unique_control_inds)]

        return type(self)(control, data)

    def __repr__(self):
        return f'<{self.__class__.__name__}\n' \
               f' {self.control.__repr__()}\n' \
               f' {self.data.__repr__()}\n' \
               f'>'

    def to_days(self):
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
            yield type(self)(control=control, data=data)

    @classmethod
    def from_packets(cls, packets, eng_packets):
        control = Control.from_packets(packets)
        data = Data.from_packets(packets)
        return cls(control, data)

    @classmethod
    def from_fits(cls, fitspath):
        header = fits.getheader(fitspath)
        control = QTable.read(fitspath, hdu='CONTROL')
        data = QTable.read(fitspath, hdu='DATA')
        obs_beg = Time(header['DATE_OBS'])
        data['time'] = (data['time'] + obs_beg)
        return cls(control=control, data=data)

    def get_energies(self):
        if 'energy_bin_edge_mask' in self.control.colnames:
            energies = get_energies_from_mask(self.control['energy_bin_edge_mask'][0])
        elif 'energy_bin_mask' in self.control.colnames:
            energies = get_energies_from_mask(self.control['energy_bin_mask'][0])
        else:
            energies = get_energies_from_mask()

        return energies


class LightCurve(Product):
    """
    Quick Look Light Curve data product.
    """
    def __init__(self, control=None, data=None):
        super().__init__(control=control, data=data)
        self.name = 'lightcurve'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        control = Control.from_packets(packets)
        control['detector_mask'] = _get_detector_mask(packets)
        control['pixel_mask'] = _get_pixel_mask(packets)
        control['energy_bin_edge_mask'] = _get_energy_bins(packets, 'NIX00266', 'NIXD0107')
        control['compression_scheme_counts_skm'] = \
            _get_compression_scheme(packets, 'NIXD0101', 'NIXD0102', 'NIXD0103')
        control['compression_scheme_triggers_skm'] = \
            _get_compression_scheme(packets, 'NIXD0104', 'NIXD0105', 'NIXD0106')
        control['num_energies'] = _get_num_energies(packets)
        control['num_samples'] = np.array(packets['NIX00271'])[
            np.cumsum(control['num_energies']) - 1]

        time, duration = control._get_time()
        # Map a given entry back to the control info through index
        control_indices = np.hstack([np.full(ns, cind) for ns, cind in
                                     control[['num_samples', 'index']]])

        cs, ck, cm = control['compression_scheme_counts_skm'][0]
        counts, counts_var = decompress(packets['NIX00272'], s=cs, k=ck, m=cm, return_variance=True)

        ts, tk, tm = control['compression_scheme_triggers_skm'][0]
        triggers, triggers_var = decompress(packets['NIX00274'], s=ts, k=tk, m=tm,
                                            return_variance=True)

        flat_indices = np.hstack((0, np.cumsum([*control['num_samples']]) *
                                  control['num_energies'])).astype(int)
        counts_reformed = [
            np.array(counts[flat_indices[i]:flat_indices[i + 1]]).reshape(n_eng, n_sam)
            for i, (n_sam, n_eng) in enumerate(control[['num_samples', 'num_energies']])]

        counts_var_reformed = [
            np.array(counts_var[flat_indices[i]:flat_indices[i + 1]]).reshape(n_eng, n_sam)
            for i, (n_sam, n_eng) in enumerate(control[['num_samples', 'num_energies']])]

        counts = np.hstack(counts_reformed).T
        counts_var = np.hstack(counts_var_reformed).T

        data = Data()
        data['control_index'] = control_indices
        data['time'] = time
        data['timedel'] = duration
        data['triggers'] = triggers
        data['triggers_err'] = np.sqrt(triggers_var)
        data['rcr'] = packets['NIX00276']
        data['counts'] = counts * u.ct
        data['counts_err'] = np.sqrt(counts_var) * u.ct

        return cls(control=control, data=data)

    def __repr__(self):
        return f'{self.name}, {self.level}\n' \
               f'{self.obs_beg.fits}, {self.obs_end}\n ' \
               f'{len(self.control)}, {len(self.data)}'


class Background(Product):
    """
    Background product.
    """

    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'background'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        control = Control.from_packets(packets)

        # Control
        control['energy_bin_mask'] = _get_energy_bins(packets, 'NIX00266', 'NIXD0111')
        control['compression_scheme_background_skm'] = _get_compression_scheme(packets, 'NIXD0108',
                                                                               'NIXD0109',
                                                                               'NIXD0110')
        control['compression_scheme_triggers_skm'] = _get_compression_scheme(packets, 'NIXD0112',
                                                                             'NIXD0113', 'NIXD0114')

        control['num_energies'] = _get_num_energies(packets)
        control['num_samples'] = np.array(packets['NIX00277'])[
            np.cumsum(control['num_energies']) - 1]

        time, duration = control._get_time()
        # Map a given entry back to the control info through index
        control_indices = np.hstack([np.full(ns, cind) for ns, cind in
                                     control[['num_samples', 'index']]])

        # Data
        bs, bk, bm = control['compression_scheme_background_skm'][0]
        counts, counts_var = decompress(packets['NIX00278'], s=bs, k=bk, m=bm, return_variance=True)

        flat_indices = np.hstack((0, np.cumsum([*control['num_samples']]) *
                                  control['num_energies'])).astype(int)

        counts_reformed = [
            np.array(counts[flat_indices[i]:flat_indices[i + 1]]).reshape(n_eng, n_sam)
            for i, (n_sam, n_eng) in enumerate(control[['num_samples', 'num_energies']])]

        counts_var_reformed = [
            np.array(counts_var[flat_indices[i]:flat_indices[i + 1]]).reshape(n_eng, n_sam)
            for i, (n_sam, n_eng) in enumerate(control[['num_samples', 'num_energies']])]

        counts = np.hstack(counts_reformed).T
        counts_var = np.hstack(counts_var_reformed).T

        ts, tk, tm = control['compression_scheme_triggers_skm'][0]
        triggers, triggers_var = decompress(packets['NIX00274'], s=ts, k=tk, m=tm,
                                            return_variance=True)

        data = Data()
        data['control_index'] = control_indices
        data['time'] = time
        data['timedel'] = duration
        data['background'] = counts * u.ct
        data['background_err'] = np.sqrt(counts_var) * u.ct
        data['triggers'] = triggers
        data['triggers_err'] = np.sqrt(triggers_var)

        return cls(control=control, data=data)


class Spectra(Product):
    """
    Spectra product.
    """

    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'spectrum'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        # Header
        control = Control.from_packets(packets)

        # Control
        control['pixel_mask'] = _get_pixel_mask(packets)
        control['compression_scheme_spectra_skm'] = _get_compression_scheme(packets, 'NIXD0115',
                                                                            'NIXD0116', 'NIXD0117')
        control['compression_scheme_triggers_skm'] = _get_compression_scheme(packets, 'NIXD0112',
                                                                             'NIXD0113', 'NIXD0114')
        # Fixed for spectra
        num_energies = 32
        control['num_energies'] = num_energies
        control['num_samples'] = np.array(packets['NIX00089'])

        # Due to the way packets are split up full contiguous block of detector 1-32 are not always
        # down-linked to the ground so need to pad the array to write to table and later fits
        total_samples = control['num_samples'].sum()
        full, partial = divmod(total_samples, 32)
        pad_after = 0
        if partial != 0:
            pad_after = 32 - partial

        control_indices = np.pad(np.hstack([np.full(ns, cind) for ns, cind in
                                            control[['num_samples', 'index']]]), (0, pad_after),
                                 constant_values=-1)
        control_indices = control_indices.reshape(-1, 32)

        duration, time = cls._get_time(control, num_energies, packets, pad_after)

        # sample x detector x energy
        # counts = np.array([eng_packets.get('NIX00{}'.format(i)) for i in range(452, 484)],
        #                   np.uint32).T * u.ct
        ss, sk, sm = control['compression_scheme_spectra_skm'][0]
        counts, counts_var = zip(
            *[decompress(packets.get('NIX00{}'.format(i)), s=ss, k=sk, m=sm, return_variance=True)
              for i in range(452, 484)])
        counts = np.vstack(counts).T
        counts_var = np.vstack(counts_var).T

        counts = np.pad(counts, ((0, pad_after), (0, 0)), constant_values=0)
        counts_var = np.pad(counts_var, ((0, pad_after), (0, 0)), constant_values=0)

        ts, tk, tm = control['compression_scheme_triggers_skm'][0]
        triggers, triggers_var = decompress(packets.get('NIX00484'), s=ts, k=tk, m=tm,
                                     return_variance=True)

        triggers = np.pad(triggers, (0, pad_after), constant_values=0)
        triggers_var = np.pad(triggers_var, (0, pad_after), constant_values=0)

        detector_index = np.pad(np.array(packets.get('NIX00100'), np.int16), (0, pad_after),
                                constant_values=-1)
        num_integrations = np.pad(np.array(packets.get('NIX00485'), np.uint16), (0, pad_after),
                                  constant_values=0)

        # Data
        data = Data()
        data['control_index'] = control_indices[:, 0]
        data['time'] = time[:, 0]
        data['timedel'] = duration[:, 0]
        data['detector_index'] = detector_index.reshape(-1, 32) * u.ct
        data['spectra'] = counts.reshape(-1, 32, num_energies) * u.ct
        data['spectra_err'] = np.sqrt(counts_var.reshape(-1, 32, num_energies))
        data['triggers'] = triggers.reshape(-1, num_energies)
        data['triggers_err'] = np.sqrt(triggers_var.reshape(-1, num_energies))
        data['num_integrations'] = num_integrations.reshape(-1, num_energies)

        return cls(control=control, data=data)

    @classmethod
    def _get_time(cls, control, num_energies, packets, pad_after):
        times = []
        durations = []
        start = 0
        for i, (ns, it) in enumerate(control['num_samples', 'integration_time']):
            off_sets = np.array(packets.get('NIX00485')[start:start + ns]) * it
            base_time = Time(scet_to_datetime(
                f'{control["scet_coarse"][i]}:{control["scet_fine"][i]}'))
            start_times = base_time + off_sets
            end_times = base_time + off_sets + it
            cur_time = start_times + (end_times - start_times) / 2
            times.extend(cur_time)
            durations.extend([it]*ns)
            start += ns
        time = Time(times)
        time = Time(np.pad(time.datetime64, (0, pad_after), constant_values=time[-1].datetime64))
        time = time.reshape(-1, num_energies)
        duration = np.pad(np.hstack(durations), (0, pad_after)).reshape(-1, num_energies) * it.unit
        return duration, time


class Variance(Product):
    """
    Variance product.
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'variance'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        # Header
        control = Control.from_packets(packets)

        # Control
        control['samples_per_variance'] = np.array(packets.get('NIX00279'), np.ubyte)
        control['pixel_mask'] = _get_pixel_mask(packets)
        control['detector_mask'] = _get_detector_mask(packets)
        control['compression_scheme_variance_skm'] = _get_compression_scheme(packets, 'NIXD0118',
                                                                             'NIXD0119', 'NIXD0120')
        energy_masks = np.array([
            [bool(int(x)) for x in format(packets.get('NIX00282')[i], '032b')]
            for i in range(len(packets.get('NIX00282')))])

        control['energy_bin_mask'] = energy_masks
        control['num_energies'] = 1
        control['num_samples'] = packets.get('NIX00280')

        time, duration = control._get_time()
        # Map a given entry back to the control info through index
        control_indices = np.hstack([np.full(ns, cind) for ns, cind in
                                     control[['num_samples', 'index']]])

        vs, vk, vm = control['compression_scheme_variance_skm'][0]
        variance, variance_var = decompress(packets.get('NIX00281'), s=vs, k=vk, m=vm,
                                            return_variance=True)

        # Data
        data = Data()
        data['time'] = time
        data['timedel'] = duration
        data['control_index'] = control_indices
        data['variance'] = variance
        data['variance_err'] = np.sqrt(variance_var)

        return cls(control=control, data=data)


class CalibrationSpectra(Product):
    """
    Calibration Spectra data product.
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'calibration_spectrum'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        control = Control.from_packets(packets)

        control['integration_time'] = (np.array(packets['NIX00122'], np.uint32) + 1) * 0.1 * u.s
        # control['obs_beg'] = control['obs_utc']
        # control['.obs_end'] = control['obs_beg'] + timedelta(seconds=control['duration'].astype('float'))
        # control['.obs_avg'] = control['obs_beg'] + (control['obs_end'] - control['obs_beg']) / 2

        # Control
        control['quiet_time'] = np.array(packets['NIX00123'], np.uint16)
        control['live_time'] = np.array(packets['NIX00124'], np.uint32)
        control['average_temperature'] = np.array(packets['NIX00125'], np.uint16)
        control['detector_mask'] = _get_detector_mask(packets)
        control['pixel_mask'] = _get_pixel_mask(packets)
        control['subspectrum_mask'] = _get_sub_spectrum_mask(packets)
        control['compression_scheme_counts_skm'] = _get_compression_scheme(packets, 'NIXD0126',
                                                                           'NIXD0127', 'NIXD0128')
        subspec_data = {}
        j = 129
        for subspec, i in enumerate(range(300, 308)):
            subspec_data[subspec + 1] = {'num_points': packets.get(f'NIXD0{j}')[0],
                                         'num_summed_channel': packets.get(f'NIXD0{j + 1}')[0],
                                         'lowest_channel': packets.get(f'NIXD0{j + 2}')[0]}
            j += 3

        control['num_samples'] = np.array(packets.get('NIX00159'), np.uint16)
        # control.remove_column('index')
        # control = unique(control)
        # control['index'] = np.arange(len(control))

        control['subspec_num_points'] = np.array(
            [v['num_points'] for v in subspec_data.values()]).reshape(1, -1)
        control['subspec_num_summed_channel'] = np.array(
            [v['num_summed_channel'] for v in subspec_data.values()]).reshape(1, -1)
        control['subspec_lowest_channel'] = np.array(
            [v['lowest_channel'] for v in subspec_data.values()]).reshape(1, -1)

        subspec_index = np.argwhere(control['subspectrum_mask'][0].flatten() == 1)
        num_sub_spectra = control['subspectrum_mask'].sum(axis=1)
        sub_channels = [np.arange(control['subspec_num_points'][0, index] + 1)
                        * (control['subspec_num_summed_channel'][0, index] + 1)
                        + control['subspec_lowest_channel'][0, index] for index in subspec_index]
        channels = list(chain(*[ch.tolist() for ch in sub_channels]))
        control['num_channels'] = len(channels)

        # Data
        data = Data()
        data['control_index'] = [0]
        data['time'] = (Time(scet_to_datetime(f"{control['scet_coarse'][0]}"
                                              f":{control['scet_fine'][0]}"))
                        + control['integration_time'][0]/2).reshape(1)
        data['timedel'] = control['integration_time'][0]
        # data['detector_id'] = np.array(packets.get('NIXD0155'), np.ubyte)
        # data['pixel_id'] = np.array(packets.get('NIXD0156'), np.ubyte)
        # data['subspec_id'] = np.array(packets.get('NIXD0157'), np.ubyte)
        num_spec_points = np.array(packets.get('NIX00146'))

        cs, ck, cm = control['compression_scheme_counts_skm'][0]
        counts, counts_var = decompress(packets.get('NIX00158'), s=cs, k=ck, m=cm,
                                        return_variance=True)


        counts_rebinned = np.apply_along_axis(rebin_proportional, 1,
                                              counts.reshape(-1, len(channels)), channels,
                                              np.arange(1025))

        counts_var_rebinned = np.apply_along_axis(rebin_proportional, 1,
                                                  counts_var.reshape(-1, len(channels)), channels,
                                                  np.arange(1025))

        dids = np.array(packets.get('NIXD0155'), np.ubyte).reshape(-1, num_sub_spectra[0])[:, 0]
        pids = np.array(packets.get('NIXD0156'), np.ubyte).reshape(-1, num_sub_spectra[0])[:, 0]

        full_counts = np.zeros((32, 12, 1024))
        full_counts[dids, pids] = counts_rebinned
        full_counts_var = np.zeros((32, 12, 1024))
        full_counts_var[dids, pids] = counts_var_rebinned
        data['counts'] = full_counts.reshape((1, *full_counts.shape))
        data['counts_err'] = np.sqrt(full_counts_var).reshape((1, *full_counts_var.shape))

        return cls(control=control, data=data)


class FlareFlagAndLocation(Product):
    """
    Flare flag and location product
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'flare_flag'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets):
        control = Control.from_packets(packets)
        control['num_samples'] = packets.get('NIX00089')

        control_indices = np.hstack([np.full(ns, cind) for ns, cind in
                                     control[['num_samples', 'index']]])

        time, duration = control._get_time()

        # DATA
        data = Data()
        data['control_index'] = control_indices
        data['time'] = time
        data['duration'] = duration
        data['loc_z'] = np.array(packets['NIX00283'], np.int16)
        data['loc_y'] = np.array(packets['NIX00284'], np.int16)
        data['thermal_index'] = np.array(packets['NIXD0061'], np.uint16)
        data['non_thermal_index'] = np.array(packets['NIXD0060'], np.uint16)
        data['location_status'] = np.array(packets['NIXD0059'], np.uint16)

        return cls(control=control, data=data)


class TMManagementAndFlareList(Product):
    """
    TM Management and Flare list product.
    """
    def __init__(self, control, data):
        super().__init__(control=control, data=data)
        self.name = 'flare_list'
        self.level = 'L1A'

    @classmethod
    def from_packets(cls, packets, eng_packets):
        tmp = QTable()
        tmp['scet_coarse'] = packets['coarse_time']
        tmp['scet_fine'] = packets['coarse_time']
        control = Control(tmp)
        data = Data()
        if 'parameters' in packets:

            control['ubsd_counter'] = packets.get('NIX00285')[0]
            control['pald_counter'] = packets.get('NIX00286')[0]
            control['num_samples'] = packets.get('NIX00286')[0]

            # DATA
            data['start_scet_coarse'] = packets.get('NIX00287')
            data['end_scet_coarse'] = packets.get('NIX00287')
            data['obs_utc'] = scet_to_datetime(f"{data['start_scet_coarse']}:0")
            data['highest_flareflag'] = packets.get('NIX00289')[0]
            data['tm_byte_volume'] = packets.get('NIX00290')[0]
            data['average_z_loc'] = packets.get('NIX00291')[0]
            data['average_y_loc'] = packets.get('NIX00292')[0]
            data['processing_mask'] = packets.get('NIX00293')[0]

            return cls(control=control, data=data)
        else:
            return None


def get_energies_from_mask(mask=None):
    """
    Return energy channels for
    Parameters
    ----------
    mask : list or array
        Energy bin mask

    Returns
    -------
    tuple
        Lower and high energy edges
    """

    if mask is None:
        low = [ENERGY_CHANNELS[edge]['e_lower'] for edge in range(32)]
        high = [ENERGY_CHANNELS[edge]['e_upper'] for edge in range(32)]
    elif len(mask) == 33:
        edges = np.where(np.array(mask) == 1)[0]
        channel_edges = [edges[i:i + 2].tolist() for i in range(len(edges) - 1)]
        low = []
        high = []
        for edge in channel_edges:
            l, h = edge
            low.append(ENERGY_CHANNELS[l]['e_lower'])
            high.append(ENERGY_CHANNELS[h - 1]['e_upper'])
    elif len(mask) == 32:
        edges = np.where(np.array(mask) == 1)
        low_ind = np.min(edges)
        high_ind = np.max(edges)
        low = [ENERGY_CHANNELS[low_ind]['e_lower']]
        high = [ENERGY_CHANNELS[high_ind]['e_upper']]
    else:
        raise ValueError(f'Energy mask or edges must have a length of 32 or 33 not {len(mask)}')

    return low, high


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
