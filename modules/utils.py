# utils.py

import os, sys, subprocess, numpy as np, urllib2, pandas as pd
import contextlib, time
try :
    import ujson as json
except :
    import json

def load_database(**params) :
    exist_db = os.path.join(params['dbname'], 'db_metadata.msg')
    db_columns = params['db_columns'] + params['metadata_columns'] + params['taxa_columns']
    data = pd.DataFrame(pd.read_msgpack(exist_db), columns = db_columns).fillna('-').astype(str)
    if 'deleted' in params :
        return data
    else :
        return data[data.get('deleted', '-') == '-']


try :
    import capnp
    capnp.remove_import_hook()
    
    def run_mash(data, is_read=False) :
        msh_file, ref_msh, n_thread, params = data

        minihash = capnp.load( os.path.join(os.path.dirname(os.path.realpath(__file__)), 'minihash.capnp') )
        msh = minihash.MinHash.read(open(msh_file,'r')).to_dict()
        size = msh['referenceList']['references'][0]['length64']

        neighbors = []
        if ref_msh is None :
            mash_db, sub_db_level = params['mash_db'], params['representative_level']
            report_neighbor, deep_search = params['barcode_dist'][0]*1.2, params['barcode_dist'][max(sub_db_level-1, 0)]
            
            if not is_read :
                init_dbs = [ os.path.join(mash_db, msh_db['name']) for msh_db in params['default_mash'] if size >= msh_db['min'] and size <= msh_db['max'] ]
            else :
                init_dbs = [ os.path.join(mash_db, msh_db['name']) for msh_db in params['default_mash'] ]

            msh_dbs = []
            for default_db in init_dbs :
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
except :
    pass

def get_file(url, fname) :
    for ite in xrange(30) :
        try :
            req = urllib2.urlopen(url, timeout=2)
            CHUNK = 16 * 1024
            with open(fname, 'wb') as fout :
                while True :
                    chunk = req.read(CHUNK)
                    if not chunk:
                        break
                    fout.write(chunk)
            return fname
        except :
            time.sleep(2)
    raise Exception('Downloading failed')

def load_params(argv) :
    c2 = dict([[ k.strip() for k in arg.split('=', 1)] for arg in argv[1:]])
    assert 'dbname' in c2, 'need "dbname".'
    assert os.path.join(c2['dbname'], 'dbsetting.cfg'), 'Not configured. Run db_create.py to build the framework for a database.'
    
    config = json.load(open(os.path.join(c2['dbname'], 'dbsetting.cfg')))
    if config['BIN'] != '' and not config['BIN'].endswith('/') :
        config['BIN'] = config['BIN'].strip() + '/'
    config['SPARSE'] = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    for k in c2 :
        if k in config :
            if isinstance(config[k], int) :
                c2[k] = int(c2[k])
            elif isinstance(config[k], float) :
                c2[k] = float(c2[k])
            elif isinstance(config[k], list) or isinstance(config[k], dict) :
                c2[k] = json.loads(c2[k])
    config.update(c2)
    for k in config :
        if isinstance(config[k], basestring) and '{' in config[k] :
            config[k] = config[k].format(**config)
    return config

def load_paramDict(c2) :
    assert 'dbname' in c2, 'need "dbname".'
    assert os.path.join(c2['dbname'], 'dbsetting.cfg'), 'Not configured. Run db_create.py to build the framework for a database.'
    
    config = json.load(open(os.path.join(c2['dbname'], 'dbsetting.cfg')))
    if config['BIN'] != '' and not config['BIN'].endswith('/') :
        config['BIN'] = config['BIN'].strip() + '/'
    config['SPARSE'] = os.path.realpath(__file__).replace('\\', '/').rsplit('/', 2)[0]
    for k in c2 :
        if k in config :
            if isinstance(config[k], int) :
                c2[k] = int(c2[k])
            elif isinstance(config[k], float) :
                c2[k] = float(c2[k])
            elif isinstance(config[k], list) or isinstance(config[k], dict) :
                c2[k] = json.loads(c2[k])
    config.update(c2)
    for k in config :
        if isinstance(config[k], basestring) and '{' in config[k] :
            config[k] = config[k].format(**config)
    for key in ('mash', 'minimap2', 'samtools') :
        if not os.path.isfile(config[key]) :
            config[key] = key
    return config


@contextlib.contextmanager
def get_genome_file(parse, data, folder='.') :
    outputs = []
    to_del = []
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
                record = records.loc[get_unique(records['barcode'].values, e)]
        for name, fname, lname in record.loc[:, ['index', 'file_path', 'url_path']].values :
            if not os.path.isfile(fname) :
                fname = os.path.join(folder, lname.rsplit('/', 1)[-1])
                get_file(lname, fname)
                if fname[-3:].upper() == '.GZ' :
                    subprocess.Popen('gzip -d {0}'.format(fname).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
                    fname = fname[:-3]
                to_del.append(fname)

            outputs.append([name, fname])
    yield outputs
    for fname in to_del :
        if os.path.isfile(fname) :
            os.unlink(fname)
    

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
        r = list(np.unique(d[tlevel].values, return_counts=True))
        r[0], r[1] = r[0][np.argsort(-r[1])], r[1][np.argsort(-r[1])]
        r[0], r[1] = r[0][r[1] > r[1][0]*0.1], r[1][r[1] > r[1][0]*0.1]
        g['taxonomy'].append( dict(rank=tlevel, freq=[dict(count=c, taxon=t) for c, t in zip(*r)] ))
    return g


def get_taxonomy(**param) :
    names, authority = {}, {}
    with open(os.path.join(param['taxonomy_db'], 'names.dmp')) as fin :
        for line in fin :
            part = line.strip().split('\t')
            if part[6] == 'scientific name' :
                names[part[0]] = part[2]
            elif part[6] in ('authority', 'type material', 'genbank common name', 'common name') :
                authority[part[0]] = 1
    children = {}
    nodes = {'1':[]}
    category = { t:1 for t in param['taxa_columns'] }
    with open(os.path.join(param['taxonomy_db'], 'nodes.dmp')) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            if part[4] in category :
                if part[4] != 'species' or part[0] in authority :
                    nodes[part[0]] = [[part[4], names.get(part[0], '')]]
                else :
                    nodes[part[0]] = [[part[4], '*' + names.get(part[0], '')]]
            else :
                nodes[part[0]] = []
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
