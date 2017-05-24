# msh_barcoding.py

import sys, os, subprocess, capnp, gzip, glob
import shutil, utils, pandas as pd, numpy as np
from multiprocessing import Pool

def create_db(data) :
    bt_build, prefix, id, genomes = data
    dbname = '{0}.{1}'.format(prefix, id)
    seq_tax = []
    for g in genomes :
        if ':' in g[1] :
            if g[1].upper().endswith('.GZ') :
                subprocess.Popen('wget -O {0}.tmp.gz {1}'.format(dbname, g[1]).split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE).communicate()
                fname = '{0}.tmp.gz'.format(dbname)
            else :
                subprocess.Popen('wget -O {0}.tmp {1}'.format(dbname, g[1]).split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE).communicate()
                fname = '{0}.tmp'.format(dbname)
        else :
            fname = g[1]
        show_cmd = 'zcat' if fname.upper().endswith('.GZ') else 'cat'
        show_run = subprocess.Popen("{0} {1}|tee -a '{2}'|grep '^>'".format(show_cmd, fname, dbname), shell=True, stdout=subprocess.PIPE)
        for line in iter(show_run.stdout.readline, r'') :
            seqname = line[1:].strip().split()[0]
            seq_tax.append([seqname, g[0]])
        if g[1] != fname :
            os.unlink(fname)
    with gzip.open('{0}.taxa.gz'.format(dbname), 'w') as fout :
        for s, b in seq_tax :
            fout.write('{0}\t{1}\n'.format(s, b))
    r = subprocess.Popen('{0} -o 3 {1} {1}'.format(bt_build, dbname).split(), stdout=subprocess.PIPE)
    r.communicate()
    if r.returncode == 0 :
        subprocess.Popen('gzip {0}'.format(dbname).split()).communicate()
    else :
        for fname in glob.glob(dbname + '.*') :
            os.unlink(fname)
    return [prefix, id, r.returncode]
    

if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    db_columns = [c for c in params['db_columns'] + params['metadata_columns'] + params['taxa_columns'] if c not in ('file_path', 'url_path', 'sha256')]
    
    existing_data = os.path.join(params['dbname'], 'db_metadata.msg')
    assert existing_data, 'no data in the database.'
    data = pd.read_msgpack(existing_data)
    
    if 'ref_size' not in params :
        params['ref_size'] = params['default_ref_size']
    
    if 'seqlist' not in params  :
        # No seqlist. Will generate a table of candidates. Modify and feed them back to create the database.
        representative_level = params['representative_level']
        
        dd = [d.split('.') for d in data['barcode'].tolist()]
        size = data['size'].as_matrix().astype(int)
        data['For_MapDB'] = ['T' if s < int(params['ref_size']) and d[representative_level][1:] == d[-1][1:] else 'F' for d, s in zip(dd, size)]
        pd.DataFrame(data, columns=['For_MapDB'] + db_columns).to_csv(sys.stdout, sep='\t')
    else :
        mapdb = params['default_MapDB'] if 'MapDB' not in params else params['MapDB']
        mapdb = os.path.join(params['bowtie_db'], mapdb.rsplit('/', 1)[-1])
        start_id = 0
        
        if params['seqlist'] == 'stdin' :
            fin = sys.stdin
        else :
            fin = open(params['seqlist'])
        glist = pd.read_csv(fin, delimiter='\t', dtype='str')
        fin.close()
        indices = {i:1 for i in glist.query("For_MapDB == 'T'")['index'].tolist()}
        
        if len(glob.glob(mapdb + '.*')) > 0 :
            assert params.get('mode', '') in ('overwrite', 'append'), 'Old database with same name present. You have to use a new name with "MapDB=", or choose between "mode=overwrite" and "mode=append".'
            if params.get('mode', '') == 'overwrite' :
                for fname in glob.glob(mapdb + '.*') :
                    os.unlink(fname)
            elif params.get('mode', '') == 'append' :
                for fname in glob.glob(mapdb + '.*.taxa.gz') :
                    i = int(fname.rsplit('.',3)[1])
                    if i >= start_id :
                        start_id = i + 1
                    with gzip.open(fname) as fin :
                        for line in fin :
                            indices[line.strip().split()[1]] = 2
        
        data = data[['index', 'file_path', 'size']].loc[data['index'].isin([i for i, t in indices.iteritems() if t == 1])]
        data = data.as_matrix()[np.argsort(-data['size'].as_matrix().astype(int))]
        min_file_num = int(np.ceil(np.sum(data.T[2].astype(float))/3800000000))
        
        buckets = [[0, []] for n in xrange(min_file_num)]
        id = -1
        for index, file_link, size in data :
            size, done = int(size), 0
            for id in range(id+1, len(buckets)) + range(id) :
                b = buckets[id]
                if b[0] + size <= 3800000000 :
                    b[0] += size
                    b[1].append([index, file_link, size])
                    done = 1
                    break
            if done == 0 :
                buckets.append([size, [[index, file_link, size]]])
        with open(mapdb + '.info', 'w') as fout :
            for id, bucket in enumerate(buckets) :
                for b,f,s in bucket[1] :
                    fout.write('{0}\t{1}\n'.format(b, id+start_id))
        pool = Pool(params['n_thread'])
        result = pool.imap_unordered(create_db, [[params['bowtie2_build'], mapdb, start_id + id, bucket[1]] for id, bucket in enumerate(buckets)])
        for r in result :
            if r[2] != 0 :
                print 'Database {0}.{1} FAILED with code {2}!'.format(*r)
        print 'Done'
