"""
Housekeeping fits structures.
"""
import numpy as np

from astropy.io import fits


def mini(num_samples):
    """
    Generate housekeeping mini report fits structure.

    Returns
    -------
    astropy.io.fits.HDUList
        HDU list, primary and binary extensions data, energy, control.
    """
    data_columns = (
        fits.Column(name='TIME', format='E', array=np.zeros(num_samples)),
        fits.Column(name='SW_RUNNING', format='L', array=np.zeros(num_samples)),
        fits.Column(name='INSTRUMENT_NUMBER', format='B', array=np.zeros(num_samples)),
        fits.Column(name='INSTRUMENT_MODE', format='B', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_PCB_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_FPGA_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_3V3_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_2V5_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_1V5_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_SPW_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_SPW0_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_SPW1_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='SW_VERSION', format='B', array=np.zeros(num_samples)),
        fits.Column(name='CPU_LOAD', format='B', array=np.zeros(num_samples)),
        fits.Column(name='ARCHIVE_MEMORY_USAGE', format='B', array=np.zeros(num_samples)),
        fits.Column(name='AUTONOMOUS_ASW_BOOT_STAT', format='L', array=np.zeros(num_samples)),
        fits.Column(name='MEMORY_LOAD_ENA_FLAG', format='L', array=np.zeros(num_samples)),
        fits.Column(name='IDPU_IDENTIFIER', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ACTIVE_SPW_LINK', format='L', array=np.zeros(num_samples)),
        fits.Column(name='OVERRUNS_FOR_TASKS', format='B', array=np.zeros(num_samples)),
        fits.Column(name='WATCHDOG_STATE', format='L', array=np.zeros(num_samples)),
        fits.Column(name='RECEIVED_SPW_PACKETS', format='I', array=np.zeros(num_samples)),
        fits.Column(name='REJECTED_SPW_PACKETS', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_1V5_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_REF_2V5_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_2V9_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_PSU_TEMP_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='FDIR_STATUS', format='J', array=np.zeros(num_samples)),
        fits.Column(name='FDIR_STATUS_MASK_OF_HK_TEMPERATURE', format='I',
                    array=np.zeros(num_samples)),
        fits.Column(name='FDIR_STATUS_MASK_OF_HK_VOLTAGE', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_SELFTEST_STATUS_FLAG', format='L', array=np.zeros(num_samples)),
        fits.Column(name='MEMORY_STATUS_FLAG', format='L', array=np.zeros(num_samples)),
        fits.Column(name='FDIR_STATUS_MASK_OF_HK_CURRENT', format='I', array=np.zeros(num_samples)),
        fits.Column(name='NUMBER_EXECUTED_TC', format='I', array=np.zeros(num_samples)),
        fits.Column(name='NUMBER_SENT_TM', format='I', array=np.zeros(num_samples)),
        fits.Column(name='NUMBER_FAILED_TM_GEN', format='I', array=np.zeros(num_samples))
    )

    data_coldefs = fits.ColDefs(data_columns)
    data_hdu = fits.BinTableHDU.from_columns(data_coldefs)
    primary = fits.PrimaryHDU()
    return fits.HDUList([primary, data_hdu])


def maxi(num_samples):
    """
    Generate housekeeping maxi report fits structure.

    Parameters
    ----------


    Returns
    -------
    astropy.io.fits.HDUList
        HDU list, primary and binary extensions data, energy, control.
    """

    data_columns = (
        fits.Column(name='TIME', format='E', array=np.zeros(num_samples)),
        fits.Column(name='SW_RUNNING', format='L', array=np.zeros(num_samples)),
        fits.Column(name='INSTRUMENT_NUMBER', format='B', array=np.zeros(num_samples)),
        fits.Column(name='INSTRUMENT_MODE', format='B', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_PCB_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_FPGA_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_3V3_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_2V5_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_1V5_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_SPW_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_SPW0_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_SPW1_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_REF_2V5A_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_REF_2V5B_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM01_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM02_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM03_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM04_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM05_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM06_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM07_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_TIM08_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_VSENSA_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_VSENSB_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ATT_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ATT_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_HV_01_16_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_HV_17_32_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='DET_Q1_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='DET_Q2_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='DET_Q3_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='DET_Q4_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_1V5_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_REF_2V5_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DPU_2V9_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_PSU_TEMP_T', format='I', array=np.zeros(num_samples)),
        fits.Column(name='SW_VERSION', format='B', array=np.zeros(num_samples)),
        fits.Column(name='CPU_LOAD', format='B', array=np.zeros(num_samples)),
        fits.Column(name='ARCHIVE_MEMORY_USAGE', format='B', array=np.zeros(num_samples)),
        fits.Column(name='AUTONOMOUS_ASW_BOOT_STAT', format='L', array=np.zeros(num_samples)),
        fits.Column(name='MEMORY_LOAD_ENA_FLAG', format='L', array=np.zeros(num_samples)),
        fits.Column(name='IDPU_IDENTIFIER', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ACTIVE_SPW_LINK', format='L', array=np.zeros(num_samples)),
        fits.Column(name='OVERRUNS_FOR_TASKS', format='B', array=np.zeros(num_samples)),
        fits.Column(name='WATCHDOG_STATE', format='L', array=np.zeros(num_samples)),
        fits.Column(name='RECEIVED_SPW_PACKETS', format='I', array=np.zeros(num_samples)),
        fits.Column(name='REJECTED_SPW_PACKETS', format='I', array=np.zeros(num_samples)),
        fits.Column(name='ENDIS_DETECTOR_STATUS', format='J', array=np.zeros(num_samples)),
        fits.Column(name='SPW1_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='SPW0_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='Q4_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='Q3_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='Q2_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='Q1_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ASPECT_B_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ASPECT_A_POWER_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ATT_M2_MOVING', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ATT_M1_MOVING', format='L', array=np.zeros(num_samples)),
        fits.Column(name='HV17_32_ENABLED_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='HV01_16_ENABLED_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='LV_ENABLED_STATUS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='HV1_DEPOLAR_IN_PROGRESS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='HV2_DEPOLAR_IN_PROGRESS', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ATT_AB_FLAG_OPEN', format='L', array=np.zeros(num_samples)),
        fits.Column(name='ATT_BC_FLAG_CLOSED', format='L', array=np.zeros(num_samples)),
        fits.Column(name='MED_VALUE_TRG_ACC', format='E', array=np.zeros(num_samples)),
        fits.Column(name='MAX_VALUE_OF_TRIG_ACC', format='E', array=np.zeros(num_samples)),
        fits.Column(name='HV_REGULATORS_MASK', format='B', array=np.zeros(num_samples)),
        fits.Column(name='TC_20_128_SEQ_CNT', format='I', array=np.zeros(num_samples)),
        fits.Column(name='ATTENUATOR_MOTIONS', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_PHOTOA0_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_PHOTOA1_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_PHOTOB0_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ASP_PHOTOB1_V', format='I', array=np.zeros(num_samples)),
        fits.Column(name='ATTENUATOR_CURRENTS', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_ATT_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='HK_DET_C', format='I', array=np.zeros(num_samples)),
        fits.Column(name='FDIR_FUNCTION_STATUS', format='I', array=np.zeros(num_samples))
    )

    data_coldefs = fits.ColDefs(data_columns)
    data_hdu = fits.BinTableHDU.from_columns(data_coldefs, name='DATA')
    primary = fits.PrimaryHDU()
    return fits.HDUList([primary, data_hdu])
