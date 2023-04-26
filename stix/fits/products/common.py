import numpy as np

from stix.core.logger import get_logger
from stix.spice import time_utils

logger = get_logger()

DEFAULT_ENERGY_CHANNELS = {
    0: {
        'channel_edge': 0,
        'energy_edge': 0,
        'e_lower': 0.0,
        'e_upper': 4.0,
        'bin_width': 4.0,
        'dE_E': 2.000,
        'ql_channel': None
    },
    1: {
        'channel_edge': 1,
        'energy_edge': 4,
        'e_lower': 4.0,
        'e_upper': 5.0,
        'bin_width': 1.0,
        'dE_E': 0.222,
        'ql_channel': 0
    },
    2: {
        'channel_edge': 2,
        'energy_edge': 5,
        'e_lower': 5.0,
        'e_upper': 6.0,
        'bin_width': 1.0,
        'dE_E': 0.182,
        'ql_channel': 0
    },
    3: {
        'channel_edge': 3,
        'energy_edge': 6,
        'e_lower': 6.0,
        'e_upper': 7.0,
        'bin_width': 1.0,
        'dE_E': 0.154,
        'ql_channel': 0
    },
    4: {
        'channel_edge': 4,
        'energy_edge': 7,
        'e_lower': 7.0,
        'e_upper': 8.0,
        'bin_width': 1.0,
        'dE_E': 0.133,
        'ql_channel': 0
    },
    5: {
        'channel_edge': 5,
        'energy_edge': 8,
        'e_lower': 8.0,
        'e_upper': 9.0,
        'bin_width': 1.0,
        'dE_E': 0.118,
        'ql_channel': 0
    },
    6: {
        'channel_edge': 6,
        'energy_edge': 9,
        'e_lower': 9.0,
        'e_upper': 10.0,
        'bin_width': 1.0,
        'dE_E': 0.105,
        'ql_channel': 0
    },
    7: {
        'channel_edge': 7,
        'energy_edge': 10,
        'e_lower': 10.0,
        'e_upper': 11.0,
        'bin_width': 1.0,
        'dE_E': 0.095,
        'ql_channel': 1
    },
    8: {
        'channel_edge': 8,
        'energy_edge': 11,
        'e_lower': 11.0,
        'e_upper': 12.0,
        'bin_width': 1.0,
        'dE_E': 0.087,
        'ql_channel': 1
    },
    9: {
        'channel_edge': 9,
        'energy_edge': 12,
        'e_lower': 12.0,
        'e_upper': 13.0,
        'bin_width': 1.0,
        'dE_E': 0.080,
        'ql_channel': 1
    },
    10: {
        'channel_edge': 10,
        'energy_edge': 13,
        'e_lower': 13.0,
        'e_upper': 14.0,
        'bin_width': 1.0,
        'dE_E': 0.074,
        'ql_channel': 1
    },
    11: {
        'channel_edge': 11,
        'energy_edge': 14,
        'e_lower': 14.0,
        'e_upper': 15.0,
        'bin_width': 1.0,
        'dE_E': 0.069,
        'ql_channel': 1
    },
    12: {
        'channel_edge': 12,
        'energy_edge': 15,
        'e_lower': 15.0,
        'e_upper': 16.0,
        'bin_width': 1.0,
        'dE_E': 0.065,
        'ql_channel': 2
    },
    13: {
        'channel_edge': 13,
        'energy_edge': 16,
        'e_lower': 16.0,
        'e_upper': 18.0,
        'bin_width': 1.0,
        'dE_E': 0.061,
        'ql_channel': 2
    },
    14: {
        'channel_edge': 14,
        'energy_edge': 18,
        'e_lower': 18.0,
        'e_upper': 20.0,
        'bin_width': 2.0,
        'dE_E': 0.105,
        'ql_channel': 2
    },
    15: {
        'channel_edge': 15,
        'energy_edge': 20,
        'e_lower': 20.0,
        'e_upper': 22.0,
        'bin_width': 2.0,
        'dE_E': 0.095,
        'ql_channel': 2
    },
    16: {
        'channel_edge': 16,
        'energy_edge': 22,
        'e_lower': 22.0,
        'e_upper': 25.0,
        'bin_width': 3.0,
        'dE_E': 0.128,
        'ql_channel': 2
    },
    17: {
        'channel_edge': 17,
        'energy_edge': 25,
        'e_lower': 25.0,
        'e_upper': 28.0,
        'bin_width': 3.0,
        'dE_E': 0.113,
        'ql_channel': 3
    },
    18: {
        'channel_edge': 18,
        'energy_edge': 28,
        'e_lower': 28.0,
        'e_upper': 32.0,
        'bin_width': 4.0,
        'dE_E': 0.133,
        'ql_channel': 3
    },
    19: {
        'channel_edge': 19,
        'energy_edge': 32,
        'e_lower': 32.0,
        'e_upper': 36.0,
        'bin_width': 4.0,
        'dE_E': 0.118,
        'ql_channel': 3
    },
    20: {
        'channel_edge': 20,
        'energy_edge': 36,
        'e_lower': 36.0,
        'e_upper': 40.0,
        'bin_width': 4.0,
        'dE_E': 0.105,
        'ql_channel': 3
    },
    21: {
        'channel_edge': 21,
        'energy_edge': 40,
        'e_lower': 40.0,
        'e_upper': 45.0,
        'bin_width': 5.0,
        'dE_E': 0.118,
        'ql_channel': 3
    },
    22: {
        'channel_edge': 22,
        'energy_edge': 45,
        'e_lower': 45.0,
        'e_upper': 50.0,
        'bin_width': 5.0,
        'dE_E': 0.105,
        'ql_channel': 3
    },
    23: {
        'channel_edge': 23,
        'energy_edge': 50,
        'e_lower': 50.0,
        'e_upper': 56.0,
        'bin_width': 6.0,
        'dE_E': 0.113,
        'ql_channel': 4
    },
    24: {
        'channel_edge': 24,
        'energy_edge': 56,
        'e_lower': 56.0,
        'e_upper': 63.0,
        'bin_width': 7.0,
        'dE_E': 0.118,
        'ql_channel': 4
    },
    25: {
        'channel_edge': 25,
        'energy_edge': 63,
        'e_lower': 63.0,
        'e_upper': 70.0,
        'bin_width': 7.0,
        'dE_E': 0.105,
        'ql_channel': 4
    },
    26: {
        'channel_edge': 26,
        'energy_edge': 70,
        'e_lower': 70.0,
        'e_upper': 76.0,
        'bin_width': 6.0,
        'dE_E': 0.082,
        'ql_channel': 4
    },
    27: {
        'channel_edge': 27,
        'energy_edge': 76,
        'e_lower': 76.0,
        'e_upper': 84.0,
        'bin_width': 8.0,
        'dE_E': 0.100,
        'ql_channel': 4
    },
    28: {
        'channel_edge': 28,
        'energy_edge': 84,
        'e_lower': 84.0,
        'e_upper': 100.0,
        'bin_width': 16.0,
        'dE_El': 0.174,
        'ql_channel': 4
    },
    29: {
        'channel_edge': 29,
        'energy_edge': 100,
        'e_lower': 100.0,
        'e_upper': 120.0,
        'bin_width': 20.0,
        'dE_El': 0.182,
        'ql_channel': 4
    },
    30: {
        'channel_edge': 30,
        'energy_edge': 120,
        'e_lower': 120.0,
        'e_upper': 150.0,
        'bin_width': 30.0,
        'dE_El': 0.222,
        'ql_channel': 4
    },
    31: {
        'channel_edge': 31,
        'energy_edge': 150,
        'e_lower': 150.0,
        'e_upper': np.inf,
        'bin_width': np.inf,
        'dE_E': np.inf,
        'ql_channel': None
    }
}

FINE_TIME_BIN_TEST = {
    0: {
        'channel_edge': 0.0,
        'energy_edge': 0.0,
        'e_lower': 0.0,
        'e_upper': 4.0,
        'bin_width': 4.0,
        'dE_E': 4.0,
        'ql_channel': None
    },
    1: {
        'channel_edge': 1.0,
        'energy_edge': 4.0,
        'e_lower': 4.0,
        'e_upper': 4.45,
        'bin_width': 0.45,
        'dE_E': 0.4500000000000002,
        'ql_channel': 0
    },
    2: {
        'channel_edge': 2.0,
        'energy_edge': 4.45,
        'e_lower': 4.45,
        'e_upper': 4.95,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 0
    },
    3: {
        'channel_edge': 3.0,
        'energy_edge': 4.95,
        'e_lower': 4.95,
        'e_upper': 5.45,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 0
    },
    4: {
        'channel_edge': 4.0,
        'energy_edge': 5.45,
        'e_lower': 5.45,
        'e_upper': 5.95,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 0
    },
    5: {
        'channel_edge': 5.0,
        'energy_edge': 5.95,
        'e_lower': 5.95,
        'e_upper': 6.45,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 0
    },
    6: {
        'channel_edge': 6.0,
        'energy_edge': 6.45,
        'e_lower': 6.45,
        'e_upper': 6.95,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 0
    },
    7: {
        'channel_edge': 7.0,
        'energy_edge': 6.95,
        'e_lower': 6.95,
        'e_upper': 7.35,
        'bin_width': 0.4,
        'dE_E': 0.39999999999999947,
        'ql_channel': 1
    },
    8: {
        'channel_edge': 8.0,
        'energy_edge': 7.35,
        'e_lower': 7.35,
        'e_upper': 7.75,
        'bin_width': 0.4,
        'dE_E': 0.40000000000000036,
        'ql_channel': 1
    },
    9: {
        'channel_edge': 9.0,
        'energy_edge': 7.75,
        'e_lower': 7.75,
        'e_upper': 8.25,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 1
    },
    10: {
        'channel_edge': 10.0,
        'energy_edge': 8.25,
        'e_lower': 8.25,
        'e_upper': 8.75,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 1
    },
    11: {
        'channel_edge': 11.0,
        'energy_edge': 8.75,
        'e_lower': 8.75,
        'e_upper': 9.25,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 1
    },
    12: {
        'channel_edge': 12.0,
        'energy_edge': 9.25,
        'e_lower': 9.25,
        'e_upper': 10.0,
        'bin_width': 0.75,
        'dE_E': 0.75,
        'ql_channel': 2
    },
    13: {
        'channel_edge': 13.0,
        'energy_edge': 10.0,
        'e_lower': 10.0,
        'e_upper': 10.5,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 2
    },
    14: {
        'channel_edge': 14.0,
        'energy_edge': 10.5,
        'e_lower': 10.5,
        'e_upper': 11.0,
        'bin_width': 0.5,
        'dE_E': 0.5,
        'ql_channel': 2
    },
    15: {
        'channel_edge': 15.0,
        'energy_edge': 11.0,
        'e_lower': 11.0,
        'e_upper': 12.0,
        'bin_width': 1.0,
        'dE_E': 1.0,
        'ql_channel': 2
    },
    16: {
        'channel_edge': 16.0,
        'energy_edge': 12.0,
        'e_lower': 12.0,
        'e_upper': 13.0,
        'bin_width': 1.0,
        'dE_E': 1.0,
        'ql_channel': 2
    },
    17: {
        'channel_edge': 17.0,
        'energy_edge': 13.0,
        'e_lower': 13.0,
        'e_upper': 15.0,
        'bin_width': 2.0,
        'dE_E': 2.0,
        'ql_channel': 3
    },
    18: {
        'channel_edge': 18.0,
        'energy_edge': 15.0,
        'e_lower': 15.0,
        'e_upper': 18.0,
        'bin_width': 3.0,
        'dE_E': 3.0,
        'ql_channel': 3
    },
    19: {
        'channel_edge': 19.0,
        'energy_edge': 18.0,
        'e_lower': 18.0,
        'e_upper': 21.0,
        'bin_width': 3.0,
        'dE_E': 3.0,
        'ql_channel': 3
    },
    20: {
        'channel_edge': 20.0,
        'energy_edge': 21.0,
        'e_lower': 21.0,
        'e_upper': 25.0,
        'bin_width': 4.0,
        'dE_E': 4.0,
        'ql_channel': 3
    },
    21: {
        'channel_edge': 21.0,
        'energy_edge': 25.0,
        'e_lower': 25.0,
        'e_upper': 28.0,
        'bin_width': 3.0,
        'dE_E': 3.0,
        'ql_channel': 3
    },
    22: {
        'channel_edge': 22.0,
        'energy_edge': 28.0,
        'e_lower': 28.0,
        'e_upper': 32.0,
        'bin_width': 4.0,
        'dE_E': 4.0,
        'ql_channel': 3
    },
    23: {
        'channel_edge': 23.0,
        'energy_edge': 32.0,
        'e_lower': 32.0,
        'e_upper': 36.0,
        'bin_width': 4.0,
        'dE_E': 4.0,
        'ql_channel': 4
    },
    24: {
        'channel_edge': 24.0,
        'energy_edge': 36.0,
        'e_lower': 36.0,
        'e_upper': 43.0,
        'bin_width': 7.0,
        'dE_E': 7.0,
        'ql_channel': 4
    },
    25: {
        'channel_edge': 25.0,
        'energy_edge': 43.0,
        'e_lower': 43.0,
        'e_upper': 50.0,
        'bin_width': 7.0,
        'dE_E': 7.0,
        'ql_channel': 4
    },
    26: {
        'channel_edge': 26.0,
        'energy_edge': 50.0,
        'e_lower': 50.0,
        'e_upper': 59.0,
        'bin_width': 9.0,
        'dE_E': 9.0,
        'ql_channel': 4
    },
    27: {
        'channel_edge': 27.0,
        'energy_edge': 59.0,
        'e_lower': 59.0,
        'e_upper': 70.0,
        'bin_width': 11.0,
        'dE_E': 11.0,
        'ql_channel': 4
    },
    28: {
        'channel_edge': 28.0,
        'energy_edge': 70.0,
        'e_lower': 70.0,
        'e_upper': 84.0,
        'bin_width': 14.0,
        'dE_E': 14.0,
        'ql_channel': 4
    },
    29: {
        'channel_edge': 29.0,
        'energy_edge': 84.0,
        'e_lower': 84.0,
        'e_upper': 110.0,
        'bin_width': 26.0,
        'dE_E': 26.0,
        'ql_channel': 4
    },
    30: {
        'channel_edge': 30.0,
        'energy_edge': 110.0,
        'e_lower': 110.0,
        'e_upper': 150.0,
        'bin_width': 40.0,
        'dE_E': 40.0,
        'ql_channel': 4
    },
    31: {
        'channel_edge': 31.0,
        'energy_edge': 150.0,
        'e_lower': 150.0,
        'e_upper': np.inf,
        'bin_width': np.inf,
        'dE_E': np.inf,
        'ql_channel': None
    }
}


def get_energy_channels(obs_time):
    obs_datetime = time_utils.anytime(obs_time)
    #print(obs_datetime)
    if time_utils.anytime(
            '2023-01-31T00:13:00Z') <= obs_datetime <= time_utils.anytime(
                '2023-02-02T16:00:00Z'):
        return FINE_TIME_BIN_TEST
    return DEFAULT_ENERGY_CHANNELS


def _get_compression_scheme(packets, nix1, nix2, nix3):
    """
    Get the compression scheme parameters.

    Parameters
    ----------
    packets : dict
        Packets
    nix1 : str
        Parameter name for S value
    nix2 : str
        Parameter name for K value
    nix3 : str
        Parameter name for M value

    Returns
    -------
    np.ndarray
        S,K,M compression scheme parameters
    """
    comp_counts = np.array((packets[nix1], packets[nix2], packets[nix3]),
                           np.ubyte).T

    return comp_counts


def _get_energy_bins(packets, nixlower, nixuppper):
    """
    Get energy bin mask from packets

    Parameters
    ----------
    packets : dict
        Packets
    nixlower : str
        Parameter name of lower 32 bins
    nixuppper : str
        Parameters name of the upper bin

    Returns
    -------
    np.ndarray
        Full energy mask of len 33
    """
    energy_bin_mask = np.array(packets[nixlower], np.uint32)
    energy_bin_mask_uppper = np.array(packets[nixuppper], np.bool8)

    full_energy_mask = [
        format(mask, 'b').zfill(32)[::-1] + format(upper, 'b')
        for mask, upper in zip(energy_bin_mask, energy_bin_mask_uppper)
    ]
    full_energy_mask = [list(map(int, m)) for m in full_energy_mask]
    full_energy_mask = np.array(full_energy_mask).astype(np.ubyte)
    return full_energy_mask


def _get_detector_mask(packets):
    """
    Get the detector mask.
    Parameters
    ----------
    packets : dict
        Packets

    Returns
    -------
    np.ndarray
        Detector mask
    """
    detector_masks = np.array(
        [
            [bool(int(x)) for x in format(packets.get('NIX00407')[i], '032b')
             ][::-1]  # reverse ind
            for i in range(len(packets.get('NIX00407')))
        ],
        np.ubyte)

    return detector_masks


def _get_pixel_mask(packets, param_name='NIXD0407'):
    """
    Get pixel mask.

    Parameters
    ----------
    packets : dict
        Packets

    Returns
    -------
    np.ndarray
        Pixel mask
    """
    pixel_masks = np.array(
        [
            [bool(int(x)) for x in format(packets.get(param_name)[i], '012b')
             ][::-1]  # reverse ind
            for i in range(len(packets.get(param_name)))
        ],
        np.ubyte)

    return pixel_masks


def _get_num_energies(packets):
    """
    Get number of energies.

    Parameters
    ----------
    packets : dict
        Packets

    Returns
    -------
    int
        Number of energies
    """
    return packets['NIX00270']


def _get_unique(packets, param_name, dtype):
    """
    Get a unique parameter raise warning if not unique.

    Parameters
    ----------
    param_name : str
        STIX parameter name eg NIX00001
    dtype : np.dtype
        Dtype to cast to eg. np.uint16/np.uint32

    Returns
    -------
    np.ndarray
        First value even if not unique
    """
    param = np.array(packets[param_name], dtype)
    if not np.all(param == param[0]):
        logger.warning('%s has changed in complete packet sequence',
                       param_name)
    return param[0]


def _get_sub_spectrum_mask(packets):
    """
    Get subspectrum mask as bool array

    Parameters
    ----------
    packets : dict
        Merged packets

    Returns
    -------
    numpy.ndarray
        Bool array of mask
    """
    sub_spectrum_masks = np.array(
        [[bool(int(x))
          for x in format(packets.get('NIX00160')[i], '08b')][::-1]
         for i in range(len(packets.get('NIX00160')))], np.ubyte)

    return sub_spectrum_masks


def get_energies_from_mask(obs_time, mask=None):
    """
    Return energy channels for
    Parameters
    ----------
    obs_time: None
        Observation time
    mask : list or array
        Energy bin mask

    Returns
    -------
    tuple
        Lower and high energy edges
    """


    energy_channels = get_energy_channels(obs_time)

    if mask is None:
        low = [energy_channels[edge]['e_lower'] for edge in range(32)]
        high = [energy_channels[edge]['e_upper'] for edge in range(32)]
    else:
        edges = np.where(np.array(mask) == 1)[0]
        channel_edges = [
            edges[i:i + 2].tolist() for i in range(len(edges) - 1)
        ]
        low = []
        high = []
        for edge in channel_edges:
            l, h = edge
            low.append(energy_channels[l]['e_lower'])
            high.append(energy_channels[h - 1]['e_upper'])
    #elif len(mask) == 32:
    #    edges = np.where(np.array(mask) == 1)
    #    low_ind = np.min(edges)
    #    high_ind = np.max(edges)
    #    low = [energy_channels[low_ind]['e_lower']]
    #    high = [energy_channels[high_ind]['e_upper']]
    #else:
    #    raise ValueError(f'Energy mask or edges must have a length of 32 or 33 not {len(mask)}')
    return low, high


def get_energies_from_edges(obs_time, e_low, e_high, e_unit=None):
    """
    get energy bins from two lists which define science channel edges
    for L1 data eunit = None and L4 request eunit >1
    """
    energy_channels = get_energy_channels(obs_time)
    if e_unit is None:
        #L1 request
        low = [energy_channels[edge]['e_lower'] for edge in e_low]
        high = [energy_channels[edge]['e_upper'] for edge in e_high]
    else:
        #L4 request
        low, high = [], []
        ibin = 0
        ibins = int((e_high - e_low + 1) / e_unit)
        for i in range(ibins):
            start = e_low + i * e_unit
            end = start + e_unit
            print('start, end, ibin:', ibin, start, end, e_unit)
            low.append(energy_channels[start]['e_lower'])
            high.append(energy_channels[end - 1]['e_upper'])
            if start > 31 or end > 31:
                break

    return low, high
