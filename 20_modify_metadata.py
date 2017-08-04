# msh_barcoding.py

import sys, os
import utils, pandas as pd, numpy as np

def update_data(exist_db, new_data, **params) :
    db_columns = params['db_columns'] + params['metadata_columns'] + params['taxa_columns']
    alterable = {col:1 for col in params['metadata_columns'] + params['taxa_columns']}
    if os.path.isfile(exist_db) :
        existing = pd.DataFrame(pd.read_msgpack(exist_db), columns = db_columns).fillna('-').astype(str).set_index('index', False)
    else :
        existing = pd.DataFrame(columns = db_columns, dtype=str).set_index('index', False)
    with open(new_data, 'r') as fin :
        cur_loc = 0
        for line in fin :
            if line.find('\t') < 0 :
                cur_loc += len(line)
            else :
                break
        fin.seek(cur_loc)
        entries = pd.read_csv(fin, delimiter='\t', dtype=str, na_values=['na'], skip_blank_lines=True).fillna('')
    if 'url_path' not in entries.columns and 'ftp_path' in entries.columns :
        entries['url_path'] = [ '-' if p == '-' else ( p if p.endswith('.fna.gz') else p + '/' + p.rsplit('/', 1)[-1] + '_genomic.fna.gz') for p in entries['ftp_path']]

    for col in entries.columns :
        if col.startswith('#') :
            entries[col[1:].strip()] = entries[col]

    entries = pd.DataFrame(entries, columns=db_columns).fillna('-').set_index('index', False)
    
    updated = []
    for (idx, entry) in entries.iterrows() :
        if idx in existing.index :
            old = existing.loc[idx]
            for fld, value in entry.iteritems() :
                if fld in alterable and old[fld] != value and value != '' :
                    old[fld] = value
                    updated.append([index, fld, old[fld], value])

    return updated


if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    assert 'seqlist' in params, 'Please feed in a tab-delimited table in "seqlist="'
    
    exist_db = os.path.join(params['dbname'], 'db_metadata.msg')
    modified = update_data(exist_db, params['seqlist'], **params)
    pd.DataFrame(modified, columns=['#index', 'field', 'oldValue', 'newValue']).to_csv(sys.stdout, index=False, sep='\t')
