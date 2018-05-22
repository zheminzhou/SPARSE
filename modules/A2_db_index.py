# msh_barcoding.py

import sys, os, capnp, subprocess
import shutil, utils, pandas as pd, numpy as np
from multiprocessing import Pool
import time, signal

capnp.remove_import_hook()

kill_signal = False
def signal_handler(signal, frame):
    global kill_signal
    print('Received a kill signal.')
    kill_signal = True


def save2mash(inputs, codes, **params) :
    mash_db, sub_db_level = params['mash_db'], params['representative_level']    
    to_merge = {}
    
    outputs = []
    for idx, _, fmsh, _, idx2 in inputs :
        c = '.'.join(['{0}{1}'.format(t,k) for t,k in zip(params['barcode_tag'], codes[idx] + [idx2])])

        minihash = capnp.load( os.path.join(os.path.dirname(os.path.realpath(__file__)), 'minihash.capnp') )
        msh = minihash.MinHash.read(open(fmsh,'r')).to_dict()
        msh['referenceList']['references'][0]['name'] = c
        size = msh['referenceList']['references'][0]['length64']
        msh = minihash.MinHash.new_message(**msh)
        msh.write(open(fmsh, 'w'))
        outputs.append([idx, size, c, fmsh, idx2])
    
        if codes[idx][sub_db_level] == idx2 :
            for default_db in params['default_mash'] :
                if default_db['ref'] >= size :
                    filename = os.path.join(mash_db, default_db['name'])
                    if filename not in to_merge :
                        to_merge[filename] = [fmsh]
                    else :
                        to_merge[filename].append(fmsh)
                    break
        elif codes[idx][-1] == idx2 :
            msh_db = os.path.join(mash_db, c.split('.')[params['representative_level']] + '.msh')
            if msh_db not in to_merge :
                to_merge[msh_db] = [fmsh]
            else :
                to_merge[msh_db].append(fmsh)
    for db, infiles in to_merge.iteritems() :
        if os.path.exists(db) :
            subprocess.Popen('mash paste {0} {1} {2}'.format(db[:-4] + '.2', db, ' '.join(infiles)).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
            shutil.move(db[:-4] + '.2.msh', db)
        else :
            subprocess.Popen('mash paste {0} {1}'.format(db[:-4], ' '.join(infiles)).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    return outputs

def genotype_and_saving(inputs, pool, **params) :
    input_idx = {i[3]:i[0] for i in inputs}

    barcode_dist = params['barcode_dist']
    codes = {i[0]:[i[4] for b in params['barcode_dist']] for i in inputs}
    res = {}
    for r in pool.map(utils.run_mash, [[i[2], None, 1, params] for i in inputs]) :
        if len(r) > 0 :
            res[input_idx[r[0][0]]] = r[0]
    
    merged_input = os.path.join(params['dbname'], 'merged_input.msh')
    subprocess.Popen('{mash} paste {0} {1}'.format(merged_input, ' '.join([i[2] for i in inputs]), **params).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    for r in utils.run_mash([merged_input, merged_input, params['n_thread'], params]) :
        if input_idx[r[0]] > input_idx[r[1]] :
            if input_idx[r[0]] not in res or res[input_idx[r[0]]][2] > r[2] :
                r[1] = input_idx[r[1]]
                res[input_idx[r[0]]] = r
    os.unlink(merged_input)

    for idx, code in sorted(codes.iteritems()) :
        if res.get(idx, None) is not None :
            best_sim, best_hit = res[idx][2], res[idx][1]
            if isinstance(best_hit, basestring) :
                best_hit = [int(r[1:]) for r in best_hit.split('.')]
            else :
                best_hit = codes[best_hit]
            for i, d in enumerate(params['barcode_dist']) :
                if d >= best_sim :
                    code[i] = int(best_hit[i])
                else :
                    break
    return save2mash(inputs, codes, **params)

def load_data(exist_db, new_data, **params) :
    db_columns = params['db_columns'] + params['metadata_columns'] + params['taxa_columns']
    if os.path.isfile(exist_db) :
        existing = pd.DataFrame(pd.read_msgpack(exist_db), columns = db_columns).fillna('-').astype(str)
    else :
        existing = pd.DataFrame(columns = db_columns, dtype=str)
    with open(new_data, 'r') as fin :
        cur_loc = 0
        for line in fin :
            if line.find('\t') < 0 :
                cur_loc += len(line)
            else :
                break
        fin.seek(cur_loc)
        entries = pd.read_csv(fin, delimiter='\t', dtype=str, na_values=['na'], skip_blank_lines=True).fillna('-')
    if 'url_path' not in entries.columns and 'ftp_path' in entries.columns :
        entries['url_path'] = [ '-' if p == '-' else ( p if p.endswith('.fna.gz') else p + '/' + p.rsplit('/', 1)[-1] + '_genomic.fna.gz') for p in entries['ftp_path']]

    entries = entries.loc[((entries.get('url_path', '-') != '-') | (entries.get('file_path', '-') != '-')) & \
                          (entries.get('excluded_from_refseq', '-') == '-') &\
                          (entries.get('version_status', 'latest') == 'latest') & \
                          (entries.get('genome_rep', 'Full') == 'Full') ]

    for col in entries.columns :
        if col.startswith('#') :
            entries[col[1:].strip()] = entries[col]
    
    entries['version'] = entries['assembly_accession'].apply(lambda x:x.rsplit('.', 1)[-1])

    entries = pd.DataFrame(entries, columns=db_columns).fillna('-').reset_index(drop=True)
    entries['refseq_category'] = pd.Categorical(entries['refseq_category'], ['reference genome', 'representative genome', '-'])
    entries['assembly_level'] = pd.Categorical(entries['assembly_level'], ['Complete Genome', 'Chromosome', 'Draft']).fillna('Draft')
    entries.version = pd.to_numeric(entries.version, errors='coerce').fillna(1).astype(np.int64)
    
    update_rec = { r[0].rsplit('.', 1)[0]:r for r in entries[['assembly_accession', 'version', 'file_path', 'url_path']].as_matrix() }
    for row_id, r in existing.iterrows() :
        acc = r['assembly_accession'].rsplit('.', 1)[0]
        if acc in update_rec :
            ur = update_rec[acc]
            for id, tag in enumerate(['assembly_accession', 'version', 'file_path', 'url_path']) :
                if str(ur[id]) != r[tag] and ur[id] != '-':
                    r[tag] = str(ur[id])
    existing.to_msgpack(exist_db)
    cur_rec = {rr:id \
               for id, r in enumerate(existing[['assembly_accession', 'file_path', 'url_path']].as_matrix()) \
               for rr in r if rr != '-'}
    novel_ids = [id \
                 for id, r in enumerate(entries[['assembly_accession', 'file_path', 'url_path']].as_matrix()) \
                 if r[0] not in cur_rec and r[1] not in cur_rec and r[2] not in cur_rec ]
    
    entries = entries.loc[novel_ids].sort_values(by=['refseq_category', 'assembly_level', 'version',  'assembly_accession'], ascending=[True, True, False,  True]).reset_index(drop=True)
    entries['refseq_category'] = entries['refseq_category'].astype(str)
    entries['assembly_level'] = entries['assembly_level'].astype(str)
    entries = pd.DataFrame(entries, columns=db_columns, dtype=str)
    return existing, entries

def add_taxa_col(data, **params) :
    taxids = np.unique(data['taxid'].as_matrix())
    if np.sum(taxids[taxids != '-'].astype(int) > 0) == 0 :
        return data
    
    nodes = utils.get_taxonomy(**params)
    taxa = {r:{} for r in params['taxa_columns'] }
    for t in taxids :
        for r, i in nodes.get(t, []) :
            taxa[r][t] = i

    for r in params['taxa_columns'] :
        data[r] = [ (taxa[r].get(tid, ori) if ori == '-' else ori) for tid, ori in data[['taxid', r]].as_matrix() ]
    return data

def mash_proc(data) :
    idx, file_link, url_link, params = data
    if not os.path.isfile(file_link) :
        file_link = '-'
    if file_link == '-' :
        fname = url_link.rsplit('/',1)[-1]
        try:
            utils.get_file(url_link, fname)
        except :
            return idx, '', '', ''
    else :
        fname = file_link
    try :
        sha = utils.fileSHA(fname)
        if os.path.isfile(fname.rsplit('/',1)[-1] + '.msh') :
            os.unlink(fname.rsplit('/',1)[-1] + '.msh')
        msh_file = utils.get_mash(fname, fname.rsplit('/',1)[-1], is_read=False, **params)
        if file_link == '-' :
            os.unlink(fname)
        return idx, sha, msh_file, fname
    except :
        return idx, '', '', ''

def db_index(params) : 
    params = utils.load_paramDict(params)
    if params.get('update', False) :
        summary_link = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt'
        summary_file = 'assembly_summary_refseq.txt'
        utils.get_file(summary_link, summary_file)
        params['seqlist'] = summary_file
    
    if params.get('update', False) or not os.path.isfile(os.path.join(params['taxonomy_db'], 'names.dmp')) or not os.path.isfile(os.path.join(params['taxonomy_db'], 'nodes.dmp')):
        taxdump = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
        import tarfile
        utils.get_file(taxdump, 'taxdump.tar.gz')
        if not os.path.isdir(params['taxonomy_db']) :
            os.makedirs(params['taxonomy_db'])
        
        with tarfile.open('taxdump.tar.gz', 'r') as tf :
            tf.extractall(path=params['taxonomy_db'])
        os.unlink('taxdump.tar.gz')
    
    assert 'seqlist' in params, 'use seqlist to bring in a list of genomes.'
    
    exist_db = os.path.join(params['dbname'], 'db_metadata.msg')
    existing, entries = load_data(exist_db, params['seqlist'], **params)
    entries = add_taxa_col(entries, **params)
    phylum_order = [(m[0] not in ('Archaea', 'Bacteria') )*100 + \
                    (m[1] in ('Metazoa', 'Viridiplantae','nan'))*10 + \
                    (m[2] in ('nan', 'Chordata', 'Arthropoda', 'Streptophyta', 'Echinodermata', 'Platyhelminthes', 'Mollusca')) \
                    for m in entries[ params['taxa_columns'][-1:-4:-1] ].as_matrix()]
    entries = entries.loc[np.argsort(phylum_order, kind='mergesort')].reset_index(drop=True)
    
    index_id = max(existing['index'].as_matrix().astype(int))+1 if existing.shape[0] > 0 else 0

    pool, batches = Pool(params['n_thread']), params['n_thread']*3
    
    sha_dict = {c:1 for c in existing['sha256'].as_matrix()}
    sha_dict[''] = 1
    
    for group_id in np.arange(0, entries.shape[0], batches) :
        inputs2 = pool.map(mash_proc, [ [idx, record['file_path'], record['url_path'], params] for idx, record in entries.loc[group_id:(group_id+batches-1)].iterrows() ])
        inputs = []
        for i in inputs2 :
            entries.loc[i[0], 'sha256'] = i[1]
            if i[1] not in sha_dict :
                sha_dict[i[1]] = 1
                inputs.append(list(i) + [index_id])
                index_id += 1
            elif i[2] != '' :
                os.unlink(i[2])
        if not len(inputs) :
            continue
        
        if kill_signal :
            sys.exit(0)
        results = genotype_and_saving(inputs, pool, **params)
        
        for idx, size, c, fmsh, index_id2 in results :
            genome = entries.loc[idx]
            genome['index'] = str(index_id2)
            genome['barcode'] = c
            genome['size'] = str(size)
            existing = existing.append(genome)
            os.unlink(fmsh)
            print time.strftime('%X %x %Z'),':', genome['organism_name'], c
        existing.to_msgpack(exist_db)
        if kill_signal :
            sys.exit(0)

if __name__ == '__main__' :
    db_index(dict([[ k.strip() for k in arg.split('=', 1)] for arg in sys.argv[1:]]))