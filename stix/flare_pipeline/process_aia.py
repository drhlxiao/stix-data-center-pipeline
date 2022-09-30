from stix.flare_pipeline import plot_idl as p

import sys
if __name__ == '__main__':
    if len(sys.argv) not in [3,4]: 
        print('Usage plot_idl [all|aia_stix_eui|stixall] <doc_id> <end_id>')
        #p.create_all_for_all()
    else:
        opt=sys.argv[1]

        if opt=='stixall':
            start_id=int(sys.argv[2])
            end_id=int(sys.argv[3])
            p.create_stix_for_all(start_id, end_id)
        else:
            create_stix_images= bool('stix' in opt)
            create_eui= bool('eui' in opt)
            create_aia= bool('aia' in opt)
            p.plot_idl(int(sys.argv[2]),  create_stix_images, create_eui, create_aia)
