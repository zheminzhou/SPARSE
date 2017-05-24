# utils.py

import os, sys, subprocess, ujson as json, numpy as np

def load_params(argv) :
    c2 = dict([[ k.strip() for k in arg.split('=', 1)] for arg in argv[1:]])
    assert 'dbname' in c2, 'need "dbname".'
    assert os.path.join(c2['dbname'], 'dbsetting.cfg'), 'Not configured. Run db_create.py to build the framework for a database.'
    
    config = json.load(open(os.path.join(c2['dbname'], 'dbsetting.cfg')))
    config.update(c2)
    for k in config :
        if isinstance(config[k], basestring) and '{' in config[k] :
            config[k] = config[k].format(**config)
    return config

def get_genome_file(parse, data, folder='.') :
    outputs = []
    for parse in parse.split(',') :
        if parse.isdigit() :
            record = data.loc[data['index'] == parse]
        elif parse[0] in 'usrpnmeca' :
            s, e = (parse.split('-') + [None])[:2]
            if e is None :
                record = data.loc[data['barcode'].str.find(s + '.') >= 0]
            else :
                e = {d:id for id, d in enumerate('usrpnmeca')}[e]
                records = data.loc[data['barcode'].str.find(s + '.') >= 0]
                get_unique = np.vectorize(lambda d, i: d.split('.')[i][1:] == d.split('.')[-1][1:])
                record = records.loc[get_unique(records['barcode'].as_matrix(), e)]
        for name, fname, lname in record.loc[:, ['index', 'file_path', 'url_path']].as_matrix() :
            if not os.path.isfile(fname) :
                fname = os.path.join(folder, lname.rsplit('/', 1)[-1])
                if fname[-3:].upper() == '.GZ' :
                    fname = fname[:-3]
                    d = subprocess.Popen('curl {0}|gzip -d > {1}'.format(lname, fname), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                else :
                    d = subprocess.Popen('curl {0} > {1}'.format(lname, fname), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                d.communicate()
                assert d.returncode == 0, 'failed during downloading {0}'.format(lname)
            outputs.append([name, fname])
    return outputs

def get_mash(fname, prefix=None, is_read=False, **param) :
    if prefix is None :
        prefix = fname
    if os.path.exists(prefix + '.msh') :
        return prefix + '.msh'
    else :
        if is_read :
            subprocess.Popen('{mash} sketch -r -b 2 -p {n_thread} {msh_param} -o {1} {0}'.format(fname, prefix, **param).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
        else :
            subprocess.Popen('{mash} sketch -p {n_thread} {msh_param} -o {1} {0}'.format(fname, prefix, **param).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
        assert prefix + '.msh'
        return prefix + '.msh'

def run_mash(data) :
    msh_file, ref_msh, n_thread, params = data
    neighbors = []
    if ref_msh is None :
        mash_db, sub_db_level = params['mash_db'], params['representative_level']
        report_neighbor, deep_search = params['barcode_dist'][0]*1.2, params['barcode_dist'][max(sub_db_level-1, 0)]
        
        default_db = os.path.join(mash_db, 'default.msh')
        
        msh_dbs = []
        if os.path.exists(default_db) :
            msh_run = subprocess.Popen('{mash} dist -p {2} {0} {1}'.format(default_db, msh_file, n_thread, **params).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            for line in iter(msh_run.stdout.readline, r'') :
                part = line.strip().split()
                dist = float(part[2])
                if dist <= deep_search :
                    msh_dbs.append(os.path.join(mash_db, part[0].split('.')[sub_db_level] +'.msh'))
                if dist <= report_neighbor :
                    neighbors.append([part[1], part[0], dist])
            for msh_db in msh_dbs :
                if os.path.exists(msh_db) :
                    msh_run = subprocess.Popen('{mash} dist -p {2} {0} {1}'.format(msh_db, msh_file, n_thread, **params).split(), stdout=subprocess.PIPE)
                    for line in iter(msh_run.stdout.readline, r'') :
                        part = line.strip().split()
                        dist = float(part[2])
                        if dist <= report_neighbor :
                            neighbors.append([part[1], part[0], dist])
    elif os.path.exists(ref_msh) :
        msh_run = subprocess.Popen('{mash} dist -p {2} {0} {1}'.format(ref_msh, msh_file, n_thread, **params).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        for line in iter(msh_run.stdout.readline, r'') :
            part = line.strip().split()
            dist = float(part[2])
            if dist <= 1 :
                neighbors.append([part[1], part[0], dist])

    return sorted(neighbors, key=lambda x:(x[0], x[2]))

import hashlib

def fileSHA(fname) :
    def hash_bytestr_iter(bytesiter, hasher, ashexstr=True):
        for block in bytesiter:
            hasher.update(block)
        return (hasher.hexdigest() if ashexstr else hasher.digest())
    def file_as_blockiter(afile, blocksize=65536):
        with afile:
            block = afile.read(blocksize)
            while len(block) > 0:
                yield block
                block = afile.read(blocksize)

    return hash_bytestr_iter(file_as_blockiter(open(fname, 'rb')), hashlib.sha256())

def retrieve_info(group, data=None, **params) :
    if data is None :
        data = pd.read_msgpack(os.path.join(params['dbname'], 'db_metadata.msg'))

    g = {'group' : group}
    tag = group if 'a' in group else group + '.'
    d = data.loc[data.barcode.str.startswith(tag)]
    g['n_member'] = d.shape[0]
    rep = next(d.itertuples())
    g['representative'] = {'assembly_accession':rep.assembly_accession, 'organism_name':rep.organism_name}
    g['taxonomy'] = []
    for tlevel in params['taxa_columns'] :
        r = np.unique(d[tlevel].as_matrix(), return_counts=True)
        g['taxonomy'].append( dict(rank=tlevel, freq=[dict(count=c, taxon=t) for c, t in sorted(zip(r[1], r[0]))] ))
    return g
