from stix.flare_pipeline import plot_idl as p

import sys
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage plot_idl [both|aia|stix|stixall] <doc_id>')
        #p.create_all_for_all()
    else:
        opt=sys.argv[1]

        if opt=='stixall':
            p.create_stix_for_all()
        else:
            create_stix_images=bool(opt !='aia')
            create_aia =bool(opt !='stix')
            p.plot_idl(int(sys.argv[2]), create_aia, create_stix_images)
