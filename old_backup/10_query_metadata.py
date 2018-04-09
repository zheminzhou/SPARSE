# msh_barcoding.py

import sys, os
import utils, pandas as pd, numpy as np

if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    fout = open(params['seqlist'], 'w') if 'seqlist' in params else sys.stdout
    data = utils.load_database(**params)
    
    db_columns = params['db_columns'] + params['metadata_columns'] + params['taxa_columns']
    filter = { key:params[key] for key in db_columns + ['name', 'tag', 'min', 'max', 'group'] if key in params }
    for fld, value in filter.iteritems() :
        if fld in db_columns :
            data = data[ data[fld].isin(value.split(',')) ]
        elif fld == 'min' :
            data = data[ data['size'].astype(int) >= int(value) ]
        elif fld == 'max' :
            data = data[ data['size'].astype(int) <= int(value) ]
        elif fld == 'group' :
            data = data[ data['barcode'].str.contains(value) ]
        elif fld == 'tag' :
            data = data.reset_index(drop=True)
            barcodes = pd.DataFrame(data['barcode'].apply(lambda barcode:[int(b[1:]) for b in barcode.split('.')]).tolist(), columns=params['barcode_tag'])
            
            for f in value.split(';') :
                f = f.strip()
                g1, g2 = f[0], f[-1]
                if f.find('==') > 0 :
                    barcodes = barcodes[barcodes[g1] == barcodes[g2]]
                else :
                    barcodes = barcodes[barcodes[g1] != barcodes[g2]]
            data = data.loc[barcodes.index].reset_index(drop=True)

    data.to_csv(fout, index=False, sep='\t')
