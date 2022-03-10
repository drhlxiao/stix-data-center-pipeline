import pymongo
try:
    connect = pymongo.MongoClient()
    db = connect["stix"]
    collection_packets = db['packets']
    collection_raw_files = db['raw_files']
    collection_calibration = db['calibration_runs']
    collection_ql= db['quick_look']
    collection_data_request_forms= db['data_requests']
    collection_fits= db['fits']
    collection_goes= db['goes_fluxes']
    collection_qllc= db['ql_lightcurves']
    collection_qlspec= db['ql_spectra']
    collection_bsd= db['bsd']
    collection_flares= db['flares']
    collection_iors= db['iors']

    print('creating indexes for goes')
    ior_indexes=[("startTime", 1),
                ("stopTime",1),
                  ("data_request_unique_ids",1),
                  ("md5",1),
                  ]
    for index in ior_indexes:
        collection_iors.create_index(index)


    indexes=[[('unix_time',1)]]
    for index in indexes:
        collection_goes.create_index(index)

    indexes=[[('start_unix_time',1)], [('unique_id',1)], [('name',1),('SPID',1),('header_unix_time',1)]]
    for index in indexes:
        collection_bsd.create_index(index)

    indexes=[[('goes.peak_flux',1)], 
            [('peak_counts',1)], 
            [('peak_unix_time',1)], 
            [('LC_statistics.lc1.signal_max',1)], 
            [('LC_statistics.lc2.signal_max',1)], 
            [('LC_statistics.lc3.signal_max',1)], 
            [('LC_statistics.lc4.signal_max',1)], 
            ]
    for index in indexes:
        collection_flares.create_index(index)

    print('creating indexes for runs')
    if collection_raw_files:
        indexes=[[('file',1)],[('date',1)]]
        for index in indexes:
            collection_raw_files.create_index(index)
    print('creating indexes for calibration')
    if collection_calibration:
        indexes=[[('duration',1)],[('start_unix_time',1)],[('run_id',1)]]
        for index in indexes:
            collection_calibration.create_index(index)

    print('creating indexes for quicklook')
    if collection_ql:
        indexes=[[('stop_unix_time',1)],[('start_unix_time',1)],[('SPID',1)]]
        for index in indexes:
            collection_ql.create_index(index)

    if collection_qllc:
        indexes=[[('stop_unix_time',1)],[('start_unix_time',1)],[('run_id',1)]]
        for index in indexes:
            collection_qllc.create_index(index)
    if collection_qlspec:
        indexes=[[('obs_time',1)]]
        for index in indexes:
            collection_qlspec.create_index(index)


    print('creating indexes for packets')
    if collection_packets:
        indexes=[[('header.unix_time',1)],[('header.SPID',1)],[('header.service_type',1)],
                [('header.service_subtype',1)],[('header.name',1)],
                [('run_id',1)],[('header.TMTC',1)], [('hash',1)]]
        for index in indexes:
            collection_packets.create_index(index)
    if collection_data_request_forms:
        print('creating indexes for user data requests')
        indexes=[[('request_type',1), ('detector_mask',1),('pixel_mask',1)],[('request_type',1)], [('detector_mask',1)],[('pixel_mask',1),('ior_id',1),('hidden',1)]]
        for index in indexes:
            print(index)
            collection_data_request_forms.create_index(index)
    if collection_fits:
        print('creating indexes for fits')
        indexes=[[('request_id',1), ('request_id',1)]#,('pixel_mask',1)],[('detector_mask',1)], [('detector_mask',1)],[('pixel_mask',1)]]
        for index in indexes:
            print(index)
            collection_fits.create_index(index)



    connect.close()

except Exception as e:
    print(e)
