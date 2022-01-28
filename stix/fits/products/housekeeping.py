"""
House Keeping data products
"""
import numpy as np

from stix.spice.datetime import scet_to_datetime
from stix.fits.io.housekeeping import mini, maxi
SKIP_ATTRS = {'scet_coarse', 'scet_fine', 'obs_utc', 'obs_beg', 'period', 'obs_avg', 'obs_end',
              'num_samples'}

class MiniReport:
    """
    Mini house keeping reported during start up of the flight software.
    """
    def __init__(self, packets):

        self.num_samples = len(packets['coarse_time'])
        # Header
        self.scet_coarse = packets['coarse_time']
        self.scet_fine = packets['fine_time']
        self.obs_utc = scet_to_datetime(f'{self.scet_coarse[0]}:{self.scet_fine[0]}')
        self.obs_beg = self.obs_utc
        self.obs_end = scet_to_datetime(f'{self.scet_coarse[-1]}:{self.scet_fine[-1]}')
        self.obs_avg = self.obs_beg + (self.obs_end - self.obs_beg) / 2.0

        # Create array of times as dt from date_obs
        times = [scet_to_datetime(f"{packets['coarse_time'][i]}:"
                                  f"{packets['fine_time'][i]}")
                 for i in range(self.num_samples)]
        time = np.array(times) - times[0]

        # Data
        self.time = [t.total_seconds() for t in time]
        self.sw_running = packets.get('NIXD0021')
        self.instrument_number = packets.get('NIXD0022')
        self.instrument_mode = packets.get('NIXD0023')
        self.hk_dpu_pcb_t = packets.get('NIXD0025')
        self.hk_dpu_fpga_t = packets.get('NIXD0026')
        self.hk_dpu_3v3_c = packets.get('NIXD0027')
        self.hk_dpu_2v5_c = packets.get('NIXD0028')
        self.hk_dpu_1v5_c = packets.get('NIXD0029')
        self.hk_dpu_spw_c = packets.get('NIXD0030')
        self.hk_dpu_spw0_v = packets.get('NIXD0031')
        self.hk_dpu_spw1_v = packets.get('NIXD0032')
        self.sw_version = packets.get('NIXD0001')
        self.cpu_load = packets.get('NIXD0002')
        self.archive_memory_usage = packets.get('NIXD0003')
        self.autonomous_asw_boot_stat = packets.get('NIXD0166')
        self.memory_load_ena_flag = packets.get('NIXD0167')
        self.idpu_identifier = packets.get('NIXD0004')
        self.active_spw_link = packets.get('NIXD0005')
        self.overruns_for_tasks = packets.get('NIXD0168')
        self.watchdog_state = packets.get('NIXD0169')
        self.received_spw_packets = packets.get('NIXD0079')
        self.rejected_spw_packets = packets.get('NIXD0079')
        self.hk_dpu_1v5_v = packets.get('NIXD0035')
        self.hk_ref_2v5_v = packets.get('NIXD0036')
        self.hk_dpu_2v9_v = packets.get('NIXD0037')
        self.hk_psu_temp_t = packets.get('NIXD0024')
        self.fdir_status = packets.get('NIX00085')
        self.fdir_status_mask_of_hk_temperature = packets.get('NIX00161')
        self.fdir_status_mask_of_hk_voltage = packets.get('NIX00162')
        self.hk_selftest_status_flag = packets.get('NIXD0163')
        self.memory_status_flag = packets.get('NIXD0164')
        self.fdir_status_mask_of_hk_current = packets.get('NIXD0165')
        self.number_executed_tc = packets.get('NIX00166')
        self.number_sent_tm = packets.get('NIX00167')
        self.number_failed_tm_gen = packets.get('NIX00168')

    def to_hdul(self):
        """
        Create a housekeeping mini report HDUL based on the number of samples

        Returns
        -------
        `astropy.io.fits.HUDList`

        """
        hdul = mini(self.num_samples)
        for key, value in self.__dict__.items():
            if key not in SKIP_ATTRS:
                hdul[1].data[key.upper()] = getattr(self, key)

        return hdul


class MaxiReport:
    """
    Maxi house keeping reported in all modes while the flight software is running.
    """
    def __init__(self, packets):

        self.num_samples = len(packets['coarse_time'])
        # Header
        self.scet_coarse = packets['coarse_time']
        self.scet_fine = packets['fine_time']
        self.obs_utc = scet_to_datetime(f'{self.scet_coarse[0]}:{self.scet_fine[0]}')
        self.obs_beg = self.obs_utc
        self.obs_end = scet_to_datetime(f'{self.scet_coarse[-1]}:{self.scet_fine[-1]}')
        self.obs_avg = self.obs_beg + (self.obs_end - self.obs_beg) / 2.0

        # Create array of times as dt from date_obs
        times = [scet_to_datetime(f"{packets['coarse_time'][i]}:"
                                  f"{packets['fine_time'][i]}")
                 for i in range(self.num_samples)]
        time = np.array(times) - times[0]

        # Data
        self.time = [t.total_seconds() for t in time]
        self.sw_running = packets.get('NIXD0021')
        self.instrument_number = packets.get('NIXD0022')
        self.instrument_mode = packets.get('NIXD0023')
        self.hk_dpu_pcb_t = packets.get('NIXD0025')
        self.hk_dpu_fpga_t = packets.get('NIXD0026')
        self.hk_dpu_3v3_c = packets.get('NIXD0027')
        self.hk_dpu_2v5_c = packets.get('NIXD0028')
        self.hk_dpu_1v5_c = packets.get('NIXD0029')
        self.hk_dpu_spw_c = packets.get('NIXD0030')
        self.hk_dpu_spw0_v = packets.get('NIXD0031')
        self.hk_dpu_spw1_v = packets.get('NIXD0032')
        self.hk_asp_ref_2v5a_v = packets.get('NIXD0038')
        self.hk_asp_ref_2v5b_v = packets.get('NIXD0039')
        self.hk_asp_tim01_t = packets.get('NIXD0040')
        self.hk_asp_tim02_t = packets.get('NIXD0041')
        self.hk_asp_tim03_t = packets.get('NIXD0042')
        self.hk_asp_tim04_t = packets.get('NIXD0043')
        self.hk_asp_tim05_t = packets.get('NIXD0044')
        self.hk_asp_tim06_t = packets.get('NIXD0045')
        self.hk_asp_tim07_t = packets.get('NIXD0046')
        self.hk_asp_tim08_t = packets.get('NIXD0047')
        self.hk_asp_vsensa_v = packets.get('NIXD0048')
        self.hk_asp_vsensb_v = packets.get('NIXD0049')
        self.hk_att_v = packets.get('NIXD0050')
        self.hk_att_t = packets.get('NIXD0051')
        self.hk_hv_01_16_v = packets.get('NIXD0052')
        self.hk_hv_17_32_v = packets.get('NIXD0053')
        self.det_q1_t = packets.get('NIXD0054')
        self.det_q2_t = packets.get('NIXD0055')
        self.det_q3_t = packets.get('NIXD0056')
        self.det_q4_t = packets.get('NIXD0057')
        self.hk_dpu_1v5_v = packets.get('NIXD0035')
        self.hk_ref_2v5_v = packets.get('NIXD0036')
        self.hk_dpu_2v9_v = packets.get('NIXD0037')
        self.hk_psu_temp_t = packets.get('NIXD0024')
        self.sw_version = packets.get('NIXD0001')
        self.cpu_load = packets.get('NIXD0002')
        self.archive_memory_usage = packets.get('NIXD0003')
        self.autonomous_asw_boot_stat = packets.get('NIXD0166')
        self.memory_load_ena_flag = packets.get('NIXD0167')
        self.idpu_identifier = packets.get('NIXD0004')
        self.active_spw_link = packets.get('NIXD0005')
        self.overruns_for_tasks = packets.get('NIXD0168')
        self.watchdog_state = packets.get('NIXD0169')
        self.received_spw_packets = packets.get('NIXD0079')
        self.rejected_spw_packets = packets.get('NIXD0078')
        self.endis_detector_status = packets.get('NIXD0070')
        self.spw1_power_status = packets.get('NIXD0080')
        self.spw0_power_status = packets.get('NIXD0081')
        self.q4_power_status = packets.get('NIXD0082')
        self.q3_power_status = packets.get('NIXD0083')
        self.q2_power_status = packets.get('NIXD0084')
        self.q1_power_status = packets.get('NIXD0085')
        self.aspect_b_power_status = packets.get('NIXD0086')
        self.aspect_a_power_status = packets.get('NIXD0087')
        self.att_m2_moving = packets.get('NIXD0088')
        self.att_m1_moving = packets.get('NIXD0089')
        self.hv17_32_enabled_status = packets.get('NIXD0090')
        self.hv01_16_enabled_status = packets.get('NIXD0091')
        self.lv_enabled_status = packets.get('NIXD0092')
        self.hv1_depolar_in_progress = packets.get('NIXD0066')
        self.hv2_depolar_in_progress = packets.get('NIXD0067')
        self.att_ab_flag_open = packets.get('NIXD0068')
        self.att_bc_flag_closed = packets.get('NIXD0069')
        self.med_value_trg_acc = packets.get('NIX00072')
        self.max_value_of_trig_acc = packets.get('NIX00073')
        self.hv_regulators_mask = packets.get('NIXD0074')
        self.tc_20_128_seq_cnt = packets.get('NIXD0077')
        self.attenuator_motions = packets.get('NIX00076')
        self.hk_asp_photoa0_v = packets.get('NIX00078')
        self.hk_asp_photoa1_v = packets.get('NIX00079')
        self.hk_asp_photob0_v = packets.get('NIX00080')
        self.hk_asp_photob1_v = packets.get('NIX00081')
        self.attenuator_currents = packets.get('NIXD0075')
        self.hk_att_c = packets.get('NIXD0075')
        self.hk_det_c = packets.get('NIXD0058')
        self.fdir_function_status = packets.get('NIX00085')

    def to_hdul(self):
        """
        Create a housekeeping maxi report HDUL based on the number of samples

        Returns
        -------
        `astropy.io.fits.HUDList`
        """
        hdul = maxi(self.num_samples)
        for key, value in self.__dict__.items():
            if key not in SKIP_ATTRS:
                hdul[1].data[key.upper()] = getattr(self, key)

        return hdul
