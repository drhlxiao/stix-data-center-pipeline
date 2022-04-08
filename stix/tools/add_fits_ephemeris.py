
            try:
                eph_header= spice.get_fits_headers(start_time=start_dt,
                                                average_time=center_dt)
                primary_hdu.header.update(eph_header)
                #added on Feb 03, 2022
            except Exception as e:
                pass

            primary_hdu.header.update(primary_header)
