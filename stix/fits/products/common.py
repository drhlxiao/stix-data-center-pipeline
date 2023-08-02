import numpy as np

from stix.core.logger import get_logger
from stix.core.energy_bins import StixEnergyBins
from stix.spice import time_utils

logger = get_logger()




def get_energy_channels(obs_time):
    return StixEnergyBins.get_energy_bins_dict(obs_time)


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
