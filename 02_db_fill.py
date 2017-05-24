# msh_barcoding.py

import sys, os, subprocess, capnp
import shutil, utils, pandas as pd, numpy as np
from multiprocessing import Pool
import time, signal

kill_signal = False
def signal_handler(signal, frame):
    global kill_signal
    print('Received a kill signal.')
    kill_signal = True



def get_taxonomy(**param) :
    names = {}
    with open(os.path.join(param['taxonomy_db'], 'names.dmp')) as fin :
        for line in fin :
            part = line.strip().split('\t')
            if part[6] == 'scientific name' :
                names[part[0]] = part[2]
    children = {}
    nodes = {'1':[]}
    category = { t:1 for t in param['taxa_columns'] }
    with open(os.path.join(param['taxonomy_db'], 'nodes.dmp')) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            nodes[part[0]] = [[part[4], names.get(part[0], '')]] if part[4] in category else []
            if part[2] not in children :
                children[part[2]] = [part[0]]
            else :
                children[part[2]].append(part[0])
    q = ['1']
    while len(q) :
        nq = []
        for qq in q :
            nq.extend(children.get(qq, []))
            for c in children.get(qq, []) :
                nodes[c].extend(nodes[qq])
        q = nq
    return nodes


def save2mash(inputs, codes, **params) :
    capnp.remove_import_hook()
    mash_db, sub_db_level = params['mash_db'], params['representative_level']    
    default_db = os.path.join(mash_db, 'default.msh')
    to_merge = {}
    
    outputs = []
    for idx, _, fmsh, _, idx2 in inputs :
        c = '.'.join(['{0}{1}'.format(t,k) for t,k in zip(params['barcode_tag'], codes[idx] + [idx2])])

        minihash = capnp.load('minihash.capnp')
        msh = minihash.MinHash.read(open(fmsh,'r')).to_dict()
        msh['referenceList']['references'][0]['name'] = c
        size = msh['referenceList']['references'][0]['length64']
        msh = minihash.MinHash.new_message(**msh)
        msh.write(open(fmsh, 'w'))
        outputs.append([idx, size, c, fmsh, idx2])
    
        if codes[idx][sub_db_level] == idx2 :
            if default_db not in to_merge :
                to_merge[default_db] = [fmsh]
            else :
                to_merge[default_db].append(fmsh)
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
        genome_list = pd.DataFrame(pd.read_msgpack(exist_db), columns = db_columns, dtype=str)
    else :
        genome_list = pd.DataFrame(columns = db_columns, dtype=str)

    new_list = pd.read_csv(new_data, delimiter='\t', dtype=str)
    if 'url_path' not in new_list.columns and 'ftp_path' in new_list.columns :
        new_list['url_path'] = [ 'nan' if p == 'na' else ( p if p.endswith('.fna.gz') else p + '/' + p.rsplit('/', 1)[-1] + '_genomic.fna.gz') for p in new_list['ftp_path']]
    for col in new_list.columns :
        if col.startswith('#') :
            new_list[col.split(' ', 1)[1]] = new_list[col]
    new_list = pd.DataFrame(new_list, columns=db_columns, dtype=str)
    new_list = new_list.loc[(new_list['file_path'] != 'nan') | (new_list['url_path'] != 'nan')].reset_index(drop=True)

    cur_rec = {rr:id \
               for id, r in enumerate(genome_list[['assembly_accession', 'file_path', 'url_path']].as_matrix()) \
               for rr in r if rr != 'nan'}
    novel_ids = [id \
                 for id, r in enumerate(new_list[['assembly_accession', 'file_path', 'url_path']].as_matrix()) \
                 if r[0] not in cur_rec and r[1] not in cur_rec and r[2] not in cur_rec ]
    new_list = new_list.loc[novel_ids].reset_index(drop=True)
    return genome_list, new_list

def add_taxa_col(data, **params) :
    taxids = np.unique(data['taxid'].as_matrix())
    if np.sum(taxids[taxids != 'nan'].astype(int) > 0) == 0 :
        return data
    
    nodes = get_taxonomy(**params)
    taxa = {r:{} for r in params['taxa_columns'] }
    for t in taxids :
        for r, i in nodes.get(t, []) :
            taxa[r][t] = i

    for r in params['taxa_columns'] :
        data[r] = [ taxa[r].get(tid, ori) for tid, ori in data[['taxid', r]].as_matrix() ]
    return data


def reorder_inputs(new_list, **params) :
    group_ranks = params['taxa_columns'][3:]
    group = new_list[group_ranks[0]]
    for rank in group_ranks[1:] :
        group.loc[group == 'nan'] = new_list.loc[group == 'nan', rank]
    group = group.as_matrix()
    order_ids = np.zeros(shape=group.shape, dtype=int)
    group_cnt = dict([[n,0] for n in np.unique(group)]+[['nan', 0]])
    for idx, record in enumerate(new_list[group_ranks].as_matrix()) :
        order_id = np.max([group_cnt.get(n, 0) for n in record] + [group_cnt['nan']])
        order_ids[idx] = order_id
        group_cnt[group[idx]] += 1
    new_list = new_list.loc[np.argsort(order_ids,kind='mergesort')].reset_index(drop=True)
    return new_list

def mash_proc(data) :
    idx, file_link, url_link, params = data
    if not os.path.isfile(file_link) :
        file_link = 'nan'
    if file_link == 'nan' :
        fname = url_link.rsplit('/',1)[-1]
        subprocess.Popen('wget -O {1} {0}'.format(url_link, fname).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    else :
        fname = file_link
    try :
        sha = utils.fileSHA(fname)
        if os.path.isfile(fname.rsplit('/',1)[-1] + '.msh') :
            os.unlink(fname.rsplit('/',1)[-1] + '.msh')
        msh_file = utils.get_mash(fname, fname.rsplit('/',1)[-1], is_read=False, **params)
        if file_link == 'nan' :
            os.unlink(fname)
        return idx, sha, msh_file, fname
    except :
        return idx, '', '', ''

if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    
    assert 'seqlist' in params, 'use seqlist to bring in a list of genomes.'
    
    exist_db = os.path.join(params['dbname'], 'db_metadata.msg')
    genome_list, new_list = load_data(exist_db, params['seqlist'], **params)
    new_list = add_taxa_col(new_list, **params)
    phylum_order = [(m[0] not in ('Archaea', 'Bacteria') )*100 + (m[1] in ('Metazoa', 'Viridiplantae','nan'))*10 + (m[2] in ('nan', 'Chordata', 'Arthropoda', 'Streptophyta', 'Echinodermata', 'Platyhelminthes', 'Mollusca')) for m in new_list[['superkingdom', 'kingdom', 'phylum']].as_matrix()]
    new_list = new_list.loc[np.argsort(phylum_order, kind='mergesort')].reset_index(drop=True)
    
    index_id = max(genome_list['index'].as_matrix().astype(int))+1 if genome_list.shape[0] > 0 else 0

    pool, batches = Pool(params['n_thread']), params['n_thread']*3
    
    sha_dict = {c:1 for c in genome_list['sha256'].as_matrix()}
    sha_dict[''] = 1
    
    for group_id in np.arange(0, new_list.shape[0], batches) :
        inputs2 = pool.map(mash_proc, [ [idx, record['file_path'], record['url_path'], params] for idx, record in new_list.loc[group_id:(group_id+batches-1)].iterrows() ])
        inputs = []
        for i in inputs2 :
            new_list.loc[i[0], 'sha256'] = i[1]
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
            genome = new_list.loc[idx]
            genome['index'] = str(index_id2)
            genome['barcode'] = c
            genome['size'] = str(size)
            genome_list = genome_list.append(genome)
            os.unlink(fmsh)
            print time.strftime('%X %x %Z'),':', genome['organism_name'], c
        genome_list.to_msgpack(exist_db)
        if kill_signal :
            sys.exit(0)
