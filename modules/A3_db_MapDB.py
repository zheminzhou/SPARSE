import sys, os, subprocess, capnp, gzip, glob
import shutil, utils, pandas as pd, numpy as np
from multiprocessing import Pool

def create_db(data) :
    build, prefix, id, genomes, dbtype = data
    dbname = '{0}.{1}'.format(prefix, id)
    seq_tax = []
    for g in genomes :
        if os.path.isfile(g[2]) :
            fname = g[2]
        else :
            fname = g[3].rsplit('/', 1)[-1]
            try:
                utils.get_file(g[3], fname)
            except :
                pass
        show_cmd = 'gzip -cd' if fname.upper().endswith('.GZ') else 'cat'
        try:
            show_run = subprocess.Popen("{0} {1}|tee -a '{2}'|grep '^>'".format(show_cmd, fname, dbname), shell=True, stdout=subprocess.PIPE)
            for line in iter(show_run.stdout.readline, r'') :
                seqname = line[1:].strip().split()[0]
                seq_tax.append([seqname, g[0]])
            if g[2] != fname :
                os.unlink(fname)
        except :
            pass
    if dbtype == 'bowtie2' :
        r = subprocess.Popen('{0} -o 3 {1} {1}'.format(build, dbname).split(), stdout=subprocess.PIPE)
    elif dbtype == 'malt' :
        r = subprocess.Popen('{0} -i {1} -s DNA -d {1}.malt -t 20'.format(build, dbname).split(), stdout=subprocess.PIPE)
    r.communicate()
    if r.returncode == 0 :
        with gzip.open('{0}.taxa.gz'.format(dbname), 'w') as fout :
            for s, b in seq_tax :
                fout.write('{0}\t{1}\n'.format(s, b))
        subprocess.Popen('gzip {0}'.format(dbname).split()).communicate()
    else :
        for fname in glob.glob(dbname + '.*') :
            os.unlink(fname)
    return [prefix, id, r.returncode]
    
def db_MapDB(params) : 
    params = utils.load_paramDict(params)
    params['dbtype'] = params.get('dbtype', 'bowtie2')
    db_columns = [c for c in params['db_columns'] + params['metadata_columns'] + params['taxa_columns'] if c not in ('sha256')]

    assert params.get('seqlist', None) is not None, 'seqlist is required. '
    
    data = utils.load_database(**params)
    
    
    if params['seqlist'] in ('stdin', '-', '') :
        fin = sys.stdin
    else :
        fin = open(params['seqlist'])
    glist = pd.read_csv(fin, delimiter='\t', dtype='str')
    fin.close()
    
    mapdb = params['MapDB']
    mapdb = os.path.join(params['bowtie_db'], mapdb)
    start_id = 0
    
    indices = {i:1 for i in glist['index'].tolist()}
    
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
    data = data.set_index('index', drop=False)
    data['size'] = data['size'].astype(int)
    data = data.loc[[i for i, t in indices.iteritems() if t == 1]].sort_values(by=['size'], ascending=[False])
    min_file_num = int(np.ceil(np.sum(data['size']).astype(float)/3800000000))
    
    buckets = [[0, []] for n in xrange(min_file_num)]
    id = -1
    for index, size, file_path, url_path in data[['index', 'size', 'file_path', 'url_path']].as_matrix() :
        size, done = int(size), 0
        for id in range(id+1, len(buckets)) + range(id+1) :
            b = buckets[id]
            if b[0] + size <= 3800000000 :
                b[0] += size
                b[1].append([index, size, file_path, url_path])
                done = 1
                break
        if done == 0 :
            buckets.append([size, [[index, size, file_path, url_path]]])
    if params['dbtype'] == 'bowtie2' :
        pool = Pool(min(params['n_thread'], len(buckets)))
        result = pool.imap_unordered(create_db, [[params['bowtie2_build'], mapdb, start_id + id, bucket[1], params['dbtype']] for id, bucket in enumerate(buckets)])
    else :
        result = map(create_db, [[params['malt_build'], mapdb, start_id + id, bucket[1], params['dbtype']] for id, bucket in enumerate(buckets)])
    for r in result :
        if r[2] != 0 :
            print 'Database {0}.{1} FAILED with code {2}!'.format(*r)

    with open(mapdb + '.info', 'w') as fout :
        for id, bucket in enumerate(buckets) :
            for b,_,_,_ in bucket[1] :
                fout.write('{0}\t{1}\n'.format(b, id+start_id))
    print 'Done'


    if __name__ == '__main__' :
        db_MapDB( dict([[ k.strip() for k in arg.split('=', 1)] for arg in sys.argv[1:]]) )

