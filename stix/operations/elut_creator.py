"""
ELUT Creator

This script generates Energy Look-Up Tables (ELUT) from a FITS file containing ECC results.
It creates two CSV files, one for the SSW (Solar SoftWare) and another for creating IORs.

Usage:
    python elut_creator.py <ECC_filename>
    - <ECC_filename>: The filename of the FITS file containing data.

Example:
    python elut_creator.py data.fits

After running the script, ELUT files will be generated in the current directory.

Author:
    Hualin Xiao (hualin.xiao@fhnw.ch)
History:
    Apr 16, 2024, v1.0 
"""

import os
import sys
from datetime import datetime

import numpy as np
from astropy.io import fits

DEFAULT_EBIN_EDGES = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70,
                      76, 84, 100, 120, 150, np.inf]


def create_elut(ecc_fname, energy_edges=DEFAULT_EBIN_EDGES):
    """
    Create ELUT (Energy Look-Up Table) from a given FITS file.

    Parameters:
    - ecc_fname (str): The filename of the FITS file containing data.
    - energy_edges (list): List of energy bin edges.

    Returns:
    - None: Saves ELUT data to CSV files.
    """
    # Check if the length of energy bins is valid
    if len(energy_edges) != 33:
        raise ValueError('Invalid length of the energy bins')

    # Open the FITS file
    with fits.open(ecc_fname) as f:
        values = f['BestValues'].data

    # Initialize lists to store ELUT data
    elut_ssw = [['Offset (ADC)', 'Gain (ADC/keV)', 'Pixel', 'Detector'] +
                [f'ADC Edge #{ibin} - {ehigh} keV' for ibin, ehigh in enumerate(energy_edges[1:-1])]]
    elut_ops = []

    # Process each row in the FITS data
    for row in values:
        pix, offset, gain, _, _ = row
        offset *= 4  # Convert offset to ADC units
        gain = 1 / (4 * gain)  # Convert gain to ADC/keV
        to_adc = lambda x: round(offset + x / gain)  # Function to convert energy to ADC units
        ipx = pix % 12  # Pixel number
        idet = pix // 12  # Detector number
        # Generate ELUT data for SSW and OPS
        pix_elut_ssw = [offset, gain, ipx, idet] + [to_adc(x) for x in energy_edges[1:-1]]
        pix_elut_ops = [idet, ipx] + [to_adc(x) for x in energy_edges[1:-1]]
        # Append ELUT data to lists
        elut_ssw.append(pix_elut_ssw)
        elut_ops.append(pix_elut_ops)

    # Generate filenames for CSV files
    date_str = datetime.now().strftime('%Y%m%d')
    fname_ssw = f'elut_table_{date_str}.csv'
    fname_ops = f'elut_table_starlet_{date_str}.csv'

    # Create header content for SSW file
    ecc_basename = os.path.basename(ecc_fname)
    ssw_content = f'Created from {ecc_basename} at {datetime.now()}\n Channel energy edges obtained from elut_creator.py\n'
    ops_content = ''

    # Write ELUT data to CSV files
    for row in elut_ssw:
        ssw_content += ','.join(map(str, row)) + '\n'
    for row in elut_ops:
        ops_content += ','.join(map(str, row)) + '\n'

    # Write to files
    with open(fname_ssw, 'w') as fssw:
        fssw.write(ssw_content)
    with open(fname_ops, 'w') as fops:
        fops.write(ops_content)

    # Print confirmation message
    print(f'ELUT has been saved to {fname_ssw} and {fname_ops}.')


if __name__ == '__main__':
    # Check if filename is provided as command-line argument
    if len(sys.argv) == 1:
        print('Elut creator\n')
        print('Usage:\n')
        print('elut_creator <ECC_filename>\n')
    else:
        # Generate ELUT
        create_elut(sys.argv[1])

