from stix.flare_pipeline import plot_idl as p

import sys
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage plot_idl <doc_id>')
        p.create_all_for_all()
    else:
        p.plot_idl(int(sys.argv[1]), True)
