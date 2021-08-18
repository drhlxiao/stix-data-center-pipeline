import pymongo
try:
    connect = pymongo.MongoClient()
    db = connect["stix"]
    collection_packets = db['packets']
    collection_raw_files = db['raw_files']
    collection_calibration = db['calibration_runs']
    collection_ql= db['quick_look']
    collection_data_request_forms= db['bsd_req_forms']
    collection_fits= db['fits']
    collection_goes= db['goes_fluxes']
    collection_qllc= db['ql_lightcurves']
    collection_qlspec= db['ql_spectra']
    print('creating indexes for goes')

    indexes=[[('unix_time',1)]]

    for index in indexes:
        collection_goes.create_index(index)

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
        indexes=[[('request_type',1), ('detector_mask',1),('pixel_mask',1)],[('request_type',1)], [('detector_mask',1)],[('pixel_mask',1)]]
        for index in indexes:
            print(index)
            collection_data_request_forms.create_index(index)
    #if collection_fits:
    #    print('creating indexes for fits')
    #    indexes=[[('request_type',1), ('detector_mask',1),('pixel_mask',1)],[('request_type',1)], [('detector_mask',1)],[('pixel_mask',1)]]
    #    for index in indexes:
    #        print(index)
    #        collection_fits.create_index(index)



    connect.close()

except Exception as e:
    print(e)
