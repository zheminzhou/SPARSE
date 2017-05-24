# sigma integrated
import sys, subprocess, os, gzip, re, shutil, numpy as np, pandas as pd
import utils, time, glob
from multiprocessing import Pool
from scipy.stats import norm
from scipy.special import factorial
from scipy.stats.mstats import gmean
'''
0	read_id
1	tax_id
2	seq_id
3	site_id
4	mutation
5	unvariable
6	lk_similarity
7	weight_summary
8	weight_dense
9	weight_core
10	prop_similarity
'''
ln_fact = np.log(factorial(np.arange(171)))
def dnorm(freq, n) :
    if len(freq) == 0 :
        return []
    res = norm.cdf((freq-n+0.5)/np.sqrt(freq)) - norm.cdf((freq-n-0.5)/np.sqrt(freq))
    res[res < 1e-100] = -230.2585
    res[res > 0] = np.log(res[res>0])
    return res

def dpois(freq, n) :
    freq, n = np.array(freq), np.array(n)
    res = np.zeros(n.shape)
    x = ((n - np.ceil(freq)) >= 20*np.sqrt(np.ceil(freq)))
    res[x], freq[x] = -230.2585, -1
    x = (freq <= 20) & (freq > 0)
    res[x] = np.log(freq[x])*n[x] - ln_fact[n[x].astype(int)] - freq[x]
    freq[x] = -1
    res[freq > 0] = dnorm(freq[freq > 0], n[freq > 0])
    return res
    
def bam2matrix(data, MapDB, workspace, bowtie_db, mismatch=0.05, n_thread=10, **params) :
    db_list = searchDBs(bowtie_db, MapDB)

    o_t = time.time()
    # # get total numbers of reads
    if os.path.isfile(os.path.join(workspace, 'r1.fastq')) :
        with open(os.path.join(workspace, 'r1.fastq')) as fin :
            fin.seek(-10000, 2)
            line = fin.readlines()[-4]
            n_read = int(line.strip()[1:])+1
    else :        
        with gzip.open(os.path.join(workspace, 'r1.fastq.gz')) as fin :
            line = fin.readlines()[-4]
            n_read = int(line.strip()[1:])+1

    # # prepare taxonomy info for sequences
    seq_files, map_files = [], []
    for dbname in db_list :
        seq_files.append(os.path.join(workspace, dbname.rsplit('/', 1)[-1] + '.seqinfo.npy'))
        map_files.extend(glob.glob(os.path.join(workspace, dbname.rsplit('/', 1)[-1] + '.map.*.npy')))
    print time.time()-o_t, "info prepared."
    sys.stdout.flush()
    
    maps, seqinfo = np.vstack([np.load(mfile) for mfile in map_files]), np.vstack([np.load(sfile) for sfile in seq_files])

    m2 = np.zeros(shape=[maps.shape[0], maps.shape[1] + 4])
    m2[:, :-4] = maps
    maps = m2
    print time.time()-o_t, "Got mappings."
    sys.stdout.flush()

    # # prepare cluster information from source DB
    b2 = np.zeros(shape=[np.max(data['index'].astype(int))+1, len(next(data.itertuples()).barcode.split('.'))], dtype=int)
    barcodes = np.zeros(shape=[np.max(data['index'].astype(int))+1, len(next(data.itertuples()).barcode.split('.'))], dtype=int)
    b2[data['index'].as_matrix().astype(int)] = [ [int(bb[1:]) for bb in b.split('.')] for b in data['barcode'].as_matrix()]
    barcodes[:, :] = -1
    barcodes[np.unique(maps.T[1]).astype(int)] = b2[np.unique(maps.T[1]).astype(int)]
    
    # # PER match conservation weighting
    conservation = np.zeros(shape=maps.shape[0])
    for id in xrange(4) :
        taxa_cnt = barcodes[np.unique(barcodes.T[4], return_index=True)[1], id]
        taxa_cnt = np.bincount(*np.unique(taxa_cnt[taxa_cnt >= 0], return_counts=True))

        if np.max(taxa_cnt) == 1 : break
        taxa_cnt[taxa_cnt == 1] = (10.0)/(8.0-id)
        
        maps[:, 10] = taxa_cnt[barcodes[maps.T[1].astype(int), id]]
        maps[:, 9] = barcodes[maps.T[1].astype(int), id] + maps.T[0].astype(int)*(taxa_cnt.size + 1)
        c = np.unique(maps.T[9], return_inverse=True, return_counts=True)
        maps[:, 9] = c[2][c[1]]
        c = (conservation < maps.T[9]/maps.T[10])
        conservation[c] = maps[c, 9]/maps[c, 10]
    maps[:, 9] = norm.cdf((conservation-0.86)*2/0.14)
    print time.time()-o_t, "conservation score."
    sys.stdout.flush()
    
    # # similarity based weighting
    maps = maps[np.argsort(maps.T[6])[::-1]]
    read_ids = np.unique(maps.T[0], return_index=True, return_inverse=True)
    maps[:, 10] = np.exp(maps[:, 6] - maps[read_ids[1][read_ids[2]], 6])
    # remove a match if it is 4 SNP more divergent than the best match
    maps = maps[ maps.T[10] >= np.power(mismatch/(1-mismatch), 3) ]
    
    # # generalized weights for conservation
    map_score = np.bincount(maps.T[0].astype(int), weights=maps.T[10])
    map_dense = np.bincount(maps.T[0].astype(int), weights=maps.T[10]*maps.T[9])
    map_dense[map_score >0] = map_dense[map_score > 0]/map_score[map_score>0]    
    maps.T[9] = map_dense[maps.T[0].astype(int)]
    print time.time()-o_t, "similarity score."
    sys.stdout.flush()
    
    # # PER match depth weighting
    # encode every 0.5 / 2 kb 
    seqinfo = seqinfo[np.in1d(seqinfo.T[2], maps.T[1])]
    seqlen = np.bincount(seqinfo.T[0], seqinfo.T[1])
    maps.T[8] = 0
    maps[seqlen[maps.T[2].astype(int)] < 1000, 8] = -400
    for fid, fragment_len in enumerate([500, 2000]) :
        cov_score = np.zeros(shape=maps.shape[0], dtype=float)
        seq_info = np.zeros(shape=[np.max(seqinfo.T[0])+1, 5], dtype=int)
        seq_info[seqinfo.T[0], 0] = seqinfo.T[1]
        seq_info[seqinfo.T[0], 4] = seqinfo.T[2]
        seq_info.T[1] = 1+(seq_info.T[0]/fragment_len)
        prev = 0
        for s in seq_info :
            s[2:4] = [prev+1, prev+s[1]+1]
            prev = s[3]
        #
        # get depth for every region
        regions = np.zeros(shape=seq_info[-1, 3]+1)
        cov_score = seq_info[maps.T[2].astype(int), 2] + (maps.T[3]/fragment_len).astype(int)
        regions[:(np.max(cov_score).astype(int)+1)] = np.bincount(cov_score.astype(int), weights=maps.T[10])
        # -- correct for continuity
        regions += np.max(np.vstack([np.concatenate([regions[1:], [0]]), np.concatenate([[0], regions[:-1]])]), 0)
        regions[np.concatenate([seq_info.T[3], [0]])] = 0
        
        cov_score = regions[cov_score.astype(int)]
        print time.time()-o_t, "Start to run depth weighting."
        sys.stdout.flush()
        
        # tax_list = [[tax_id, region_start, region_end]]
        tax_ids = np.unique(seq_info.T[4], return_index=True)
        tax_list = np.vstack([tax_ids[0],tax_ids[1]])[:, np.argsort(tax_ids[1])].astype(float)
        tax_list = np.vstack([tax_list, np.concatenate([tax_list[1][1:], [seq_info.shape[0]]])  ]).T
        
        # get mean depth
        for tid, (tax_id, s, e) in enumerate(tax_list.astype(int)) :
            rid = np.concatenate([ np.arange(seq_info[ss][2], seq_info[ss][3]) for ss in np.arange(s, e) ])
            geometric_cov = 1.0
            if np.max(regions[rid]) > 1 :
                for ite in xrange(10) :
                    ngcov = gmean(regions[rid]+geometric_cov) - geometric_cov
                    if np.abs(ngcov-geometric_cov) < 1e-4 :
                        geometric_cov = ngcov
                        break
                    else :
                        geometric_cov = ngcov
            default_lk = dpois(geometric_cov, np.ceil(geometric_cov))
            tax_list[tid][1:] = geometric_cov, default_lk
            
        geometric_covs = np.bincount(tax_list.T[0].astype(int), tax_list.T[1])
        default_lks = np.bincount(tax_list.T[0].astype(int), tax_list.T[2])
        cov_score[np.ceil(cov_score) <= np.ceil(geometric_covs[maps.T[1].astype(int)])] = 0.0
        cov_score[cov_score > 0] = dpois(geometric_covs[maps[cov_score > 0, 1].astype(int)], np.ceil(cov_score[cov_score > 0])) - default_lks[maps[cov_score > 0, 1].astype(int)]
        maps[maps.T[8] > cov_score, 8] = cov_score[maps.T[8] > cov_score]
    
        # # generalized weights for mapping depth
    map_dense = np.bincount(maps.T[0].astype(int), weights=maps.T[10]*np.exp(maps.T[8]))
    map_dense[map_score >0] = map_dense[map_score > 0]/map_score[map_score>0]
    maps.T[8] = map_dense[maps.T[0].astype(int)]
    
    # # summarise maps.T[6] (similarity weighting),    maps.T[8] (depth weighting),    maps.T[9] (conservation weighting)
    maps.T[6] = np.exp(maps.T[6])
    maps.T[7] = maps.T[8]*maps.T[9]
    maps = maps[(maps.T[7] > 0) & (maps.T[6] > 0)]
    print time.time()-o_t, "Get all weighting."
    sys.stdout.flush()

    for mid in xrange(0, maps.shape[0], 20000000) :
        maps[mid:(mid+20000000)].dump( os.path.join(workspace, 'read_scores.{0}.npy'.format(mid/20000000)) )
    
    maps = maps[(maps.T[9] >= 0.01) & (maps.T[8] >= 0.01) & (maps.T[7]>= 0.001) ]
    picked = []
    workon = np.copy(maps)
    r_limit = max(n_read*0.000001, 3)
    while True :
        tax_match = np.bincount(workon.T[1].astype(int), weights=workon.T[10])
        if tax_match.size == 0 : break
        tax_id = np.argmax(tax_match)
        if tax_match[tax_id] < r_limit :
            break
        picked.append([tax_id, tax_match[tax_id]])
        workon = workon[(~np.in1d(workon.T[0], workon[(workon.T[1] == tax_id) & (workon.T[10] >= 0.5), 0])) & (workon.T[1] != tax_id)]
    picked = np.array(picked)
    if picked.shape[0] > 5000 :
        r_limit = max(n_read*0.00001, 10)
        picked = picked[ picked.T[1] >= r_limit ]

    print time.time()-o_t, "got candidates."
    sys.stdout.flush()

    maps = maps[np.in1d(maps.T[1], picked.T[0])]
    maps = maps[np.lexsort([-maps.T[10], maps.T[0]])]

    def write_down(fout, reads, w) :
        if len(reads) > 0 :
            fout.write('*\t{0}\t{2}\t{1}\n'.format(int(reads.keys()[0]), '\t'.join(['{0}={1}'.format(t, p) for t, p in sorted(reads.values()[0].items())]), w))
    
    with open(os.path.join(workspace, 'ipopt.qmatrix.txt'), 'w') as fout :
        fout.write('''#\t+\tMatrixName\tTotalNumberReads\tMatchedReads\tUnmatchedReads\ReadLimit
#\t@\tGenomeIndex\tGenomeName\tMatchedNumber
#\t*\tReadID\tWeight\tGenomeIndex=QValue
+\t{0}\t{1}\t{2}\t{3}\t{4}
'''.format(os.path.join(workspace, 'ipopt.qmatrix.txt'), n_read, (np.unique(maps.T[0])).size, n_read-(np.unique(maps.T[0])).size, r_limit))
        oid = {}
        for id, (t, c) in enumerate(picked) :
            fout.write('@\t{0}\t{1}\t{2}\n'.format(id, int(t), c))
            oid[t] = id
        reads, w = {}, 0.0
        for m in maps :
            if m[0] not in reads :
                write_down(fout, reads, w)
                reads = {m[0]:{}}
                w = m[7]
            if oid[m[1]] not in reads[m[0]] or m[6] > reads[m[0]][oid[m[1]]] :
                reads[m[0]][oid[m[1]]] = m[6]
        write_down(fout, reads, w)
    return workspace

def searchDBs(bowtie_db, MapDB) :
    db_list = []
    for mdb in MapDB.split(',') :
        for id in xrange(500) :
            if os.path.exists(os.path.join(bowtie_db, '{0}.{1}.rev.1.bt2'.format(mdb, id))) or os.path.exists(os.path.join(bowtie_db, '{0}.{1}.rev.1.bt2l'.format(mdb, id))) :
                db_list.append(os.path.join(bowtie_db, '{0}.{1}'.format(mdb, id)))
            else :
                break
    return db_list

def align(bowtie_db, MapDB, workspace, mismatch, r1, r2=None, n_thread=10, **params) :
    db_list = searchDBs(bowtie_db, MapDB)
    # prepare databases
    if not os.path.isdir(workspace) :
        os.mkdir(workspace)

    seqs, seqID = {}, -1
    for dbname in db_list :
        with gzip.open('{0}.taxa.gz'.format(dbname)) as fin :
            for line in fin :
                seqID += 1
                seqname, tax_id = line.strip().split()
                seqs[seqname] = [int(tax_id), seqID]

    with open(os.path.join(workspace, 'r1.fastq'), 'w') as fout :
        fin = gzip.open(r1) if r1.endswith('.gz') else open(r1)
        for id, line in enumerate(fin) :
            fout.write(line if id%4 else '@{0}\n'.format(id/4))
        fin.close()
    if r2 is None :
        read_info = '-U {0}'.format(os.path.join(workspace, 'r1.fastq'))
    else :
        read_info = '-1 {0} -2 {1}'.format(os.path.join(workspace, 'r1.fastq'), os.path.join(workspace, 'r2.fastq'))
        with open(os.path.join(workspace, 'r2.fastq'), 'w') as fout :
            fin = gzip.open(r2) if r2.endswith('.gz') else open(r2)
            for id, line in enumerate(fin) :
                fout.write(line if id%4 else '@{0}\n'.format(id/4))
            fin.close()
            
    # run bowtie2 in multiple threads
    for db_prefix in db_list :
        out_prefix = os.path.join(workspace, db_prefix.rsplit('/', 1)[-1])
        
        cmd_bt2 = '{0} -k 300 --no-unal --ignore-quals --score-min L,20,1.1 --mp 6,6 --np 6 --sensitive-local -q -p {2} -I 25 -X 800 -x {3} {1}'.format(
            params['bowtie2'], read_info, n_thread, db_prefix)
        run_bt2 = subprocess.Popen(cmd_bt2.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        prev, res, seq_info = {}, [], []
        used_names = {}
        for line in iter(run_bt2.stdout.readline, r'') :
            if line.startswith('@') :
                if line.startswith('@SQ') :
                    name, seq_len = [ p[3:] for p in line.strip().split()[1:3] ]
                    if name in seqs :
                        used_names[name] = 1
                        seq_info.append([seqs[name][1], int(seq_len), seqs[name][0]] )
            else :
                part = line.strip().split()
                if part[2] in seqs :
                    key = (part[0], int(part[1]) & 128)
                    if key not in prev :
                        s = [int(l) for l in re.findall('(\d+)S', part[5])]
                        if sum(s) > 100 or sum(s) * 2 >= len(part[9]) or (len(s) > 1 and min(s) > 3) :
                            continue
                        else :
                            prev = {key:{seqs[part[2]][0]:1}}
                    elif seqs[part[2]][0] not in prev[key] :
                        s = [int(l) for l in re.findall('(\d+)S', part[5])]
                        if sum(s) > 100 or sum(s) * 2 >= len(part[9]) or (len(s) > 1 and min(s) > 3) :
                            continue
                        else :
                            prev[key][seqs[part[2]][0]] = 1
                    else :
                        continue
                    score = int(part[11][5:]) + (0 if len(s) == 0 else max(s))
                    mut = (len(part[9])*2 - score)/8.0
                    nom = len(part[9]) - mut
                    res.append([ int(part[0]), 
                                 seqs[part[2]][0], seqs[part[2]][1], int(part[3]), 
                                 mut, nom, 0.0])
        for name in used_names: 
            seqs.pop(name, None)
        sfile = out_prefix + '.seqinfo.npy'
        np.array(seq_info).dump(sfile)
        res = np.array(res)
        if res.size > 0 :
            res[:, 6] = ( res[:, 4]*np.log(mismatch) + res[:, 5]*np.log(1-mismatch) ) - np.sum(res[:, 4:6], 1)*(mismatch*np.log(mismatch)+(1-mismatch)*np.log(1-mismatch))
            for mid in xrange(0, res.shape[0], 20000000) :
                mfile = out_prefix + '.map.{0}.npy'.format(mid/20000000)
                res[mid:(mid+20000000)].dump( mfile )
    return True

def ipopt(workspace, **params) :
    # run sigma
    subprocess.Popen('{ipopt} -t {n_thread} -i {0}'.format(os.path.join(workspace, 'ipopt.qmatrix.txt'), **params).split()).communicate()
    return os.path.join(workspace, 'ipopt.gvector.txt')

def assign_reads(data, qvector, workspace, **params) :
    maps = []
    for id in xrange(99999) :
        if os.path.isfile(os.path.join(workspace, 'read_scores.{0}.npy'.format(id))) :
            maps.append(np.load(os.path.join(workspace, 'read_scores.{0}.npy'.format(id))))
        else :
            break
    maps = np.vstack(maps)

    taxon = []
    with open(qvector) as fin :
        for line in fin :
            if line.startswith('+') :
                part = line.strip().split('\t')
                n_read, c_read, r_limit = int(part[2]), int(part[3]), float(part[-1])
            elif line.startswith('@') :
                part = line.strip().split('\t')
                if c_read * float(part[4])/100 >= r_limit :
                    taxon.append([int(part[2]), float(part[4])/100.0])
    taxon = np.array(taxon)
    taxa = np.bincount(taxon.T[0].astype(int)+1, taxon.T[1])
    maps[maps.T[1] >= taxa.shape[0]-1, 1] = -1
    maps[taxa[maps.T[1].astype(int)+1] == 0, 1] = -1

    # 4: read distance;   5: mean distance
    significant_aln = (maps.T[10] >= 0.05)
    aln_len = np.bincount(maps[significant_aln, 1].astype(int)+1, (maps.T[7]*(maps.T[4] + maps.T[5]))[significant_aln])
    aln_mut = np.bincount(maps[significant_aln, 1].astype(int)+1, (maps.T[7]*maps.T[4])[significant_aln])
    
    maps.T[4] = 100.0*maps.T[4]/(maps.T[4] + maps.T[5])
    maps.T[5] = 100.0*(aln_mut/aln_len)[maps.T[1].astype(int)+1]
    
    # justified similarity score
    maps.T[10] = maps.T[6]*taxa[maps.T[1].astype(int)+1]
    # remove low level hits
    maps = maps[np.argsort(-maps.T[10])]
    r_max_score = np.unique(maps.T[0], return_index=True, return_inverse=True)
    mean_weight = np.sum(maps[r_max_score[1], 7])/r_max_score[1].size
    
    maps.T[10] = (maps.T[10]+1e-301) / (maps[ r_max_score[1][r_max_score[2]] , 10] + 1e-300)
    maps = maps[maps.T[10] >= 0.1]

    read = np.bincount(maps.T[0].astype(int), maps.T[10])
    maps.T[10] = maps.T[10]/read[maps.T[0].astype(int)]
    
    taxa = np.bincount(maps.T[1].astype(int)+1, maps.T[10])
    maps = maps[(maps.T[1] < 0) | (taxa[maps.T[1].astype(int)+1] >= r_limit)]
    
    maps = maps[np.lexsort([maps.T[4], -maps.T[10], maps.T[0]])]
    maps[maps.T[1] < 0, 10] = 0.0
    
    read = {}
    for m in maps :
        if m[0] not in read :
            read[m[0]] = [m[8], m[9]]
        elif m[1] < 0 and read[m[0]][-1][0] < 0 :
            read[m[0]][-1][1] +=m[10]
            continue
        read[m[0]].append([m[1].astype(int), m[10], m[4], m[5]])
    with gzip.open(os.path.join(workspace, 'read_assignment.gz'), 'w') as fout :
        fout.write('#\tReadID\tReadName\tDepth_Weight\tCore_Weight\tReferenceID:Proportion:Distance:MeanDistance\n')
        fin = gzip.open(params['r1']) if params['r1'].endswith('.gz') else open(params['r1'])
        for id, line in enumerate(fin) :
            if id % 4 == 0 :
                rid = id / 4
                if rid in read :
                    fout.write('{0}\t{1}\t{2:.5f}\t{3:.5f}\t{4}\n'.format(rid, line.strip().split()[0][1:], read[rid][0], read[rid][1], '\t'.join(['{0}:{1:.2f}:{2:.2f}:{3:.2f}'.format(*x) for x in read[rid][2:]])))
                else :
                    fout.write('{0}\t{1}\n'.format(rid, line.strip().split()[0][1:]))
        

def group_major(mat, match_group, **params) :
    groups = {}
    for part in mat[['index', 'barcode'] + list(reversed(params['taxa_columns']))].as_matrix() :
        index, barcode, taxa = part[0], part[1], part[2:]
        barcodes = barcode.split('.')
        for b in barcodes[:4] :
            if b not in match_group :
                continue
            if b not in groups :
                groups[b] = [{'':0} for c in params['taxa_columns'] ]
            for t, g in zip(taxa, groups[b]) :
                if t != '' :
                    g[t] = g.get(t, 0) + 1
    for g, taxa in groups.iteritems() :
        if g[0] in 'p' :
            taxa[:] = [max(t.items(), key=lambda x:x[1])[0] for t in taxa]+['{0[0]}: {0[1]}'.format(mat.loc[mat['index'] == g[1:], ['organism_name', 'assembly_accession']].as_matrix().tolist()[0])]
        elif g[0] == 'r' :
            taxa[:] = [max(t.items(), key=lambda x:x[1])[0] for t in taxa]
        elif g[0] == 's' :
            taxa[:] = [max(t.items(), key=lambda x:x[1])[0] for t in taxa][:-1]
        elif g[0] == 'u' :
            taxa[:] = [max(t.items(), key=lambda x:x[1])[0] for t in taxa][:-2]
    return groups

def profiling(data, assign, **params) :
    basic_aln = [0, [0, 0.], [0, 0.]]
    match_ref = {}
    with gzip.open(assign) as fin :
        for line in fin :
            if line.startswith('#') : continue
            part = line.strip().split('\t')
            if len(part) < 3 :
                basic_aln[0] += 1
            else :
                dweight, cweight = float(part[2]), float(part[3])
                basic_aln[1][0] += 1
                basic_aln[1][1] += dweight
                for match in part[4:] :
                    ref, prop, d2, d1 = match.split(':')
                    if ref == '-1' :
                        basic_aln[2][0] += 1
                        basic_aln[2][1] += dweight
                    else :
                        d = float(d1) if float(d2) < float(d1)*2 else float(d2)
                        if ref not in match_ref :
                            match_ref[ref] = []
                        match_ref[ref].append([d, dweight*cweight*float(prop), dweight*float(prop)])
    match_ref = { idx:np.array(reads) for idx, reads in match_ref.iteritems() }
    barcodes = data[['index', 'barcode']].as_matrix()
    barcodes = { idx:barcode.split('.') for idx, barcode in barcodes[np.in1d(barcodes.T[0], match_ref.keys())] }
    top_cluster = {}
    for idx, clusters in barcodes.iteritems() :
        ani90 = clusters[0]
        if ani90 not in top_cluster :
            top_cluster[ani90] = np.sum(match_ref[idx].T[1:], 1)
        else :
            top_cluster[ani90][0] += np.sum(match_ref[idx].T[1])
            top_cluster[ani90][1] += np.sum(match_ref[idx].T[2])
    match_group = {}
    for ref, match in match_ref.iteritems() :
        match.T[1] *= top_cluster[barcodes[ref][0]][1]/top_cluster[barcodes[ref][0]][0]
        for dc, gg in zip(params['barcode_dist'], barcodes[ref][:-1])[3::-1] :
            if gg in match_group :
                match_group[gg][0].append(ref)
                match_group[gg][1] += np.sum(match[match.T[0] <= dc*120, 1])
            else :
                match_group[gg] = [[ref], np.sum(match[match.T[0] <= dc*120, 1])]
    groups = group_major(data, match_group, **params)
    basic_aln[0] = [basic_aln[0], basic_aln[0]*basic_aln[1][1]/basic_aln[1][0]]
    tot_weight = basic_aln[0][1] + basic_aln[1][1]
    with open(os.path.join(params['workspace'], 'profile.txt'), 'w') as fout :
        fout.write('Total\t{0}\t{1}\n'.format(basic_aln[0][0]+basic_aln[1][0], basic_aln[1][0]))
        fout.write('Unmatched\t{0:.3f}\t{1:.3f}\n'.format(basic_aln[0][1]*100.0/tot_weight, 0.0))
        fout.write('Uncertain_match\t{0:.3f}\t{1:.3f}\n'.format(basic_aln[2][1]*100.0/tot_weight, basic_aln[2][1]*100.0/basic_aln[1][1]))
        for g, (r, w) in sorted(match_group.iteritems(), key=lambda x:(x[1][1], x[0]), reverse=True) :
            if w/tot_weight >= 0.000005 :
                fout.write('{0}\t{1:.3f}\t{2:.3f}\t{3} ({4})\n'.format(g, w*100.0/tot_weight, w*100.0/basic_aln[1][1], '|'.join(groups[g]), ','.join(sorted(set(r), key=lambda x:int(x)))))

if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    assert 'workspace' in params and 'r1' in params
    existing_data = os.path.join(params['dbname'], 'db_metadata.msg')
    assert existing_data, 'no data in the database.'
    data = pd.read_msgpack(existing_data)
    
    if params.get('stage', '0') in '0' :
        align(**params)
    if params.get('stage', '0') in '01' :
        bam2matrix(data, **params)
    if params.get('stage', '0') in '012' :    
        qvector = ipopt(**params)
    if params.get('stage', '0') in '0123' :
        qvector = os.path.join(params['workspace'], 'ipopt.gvector.txt')
        assign_reads(data, qvector, **params)
    if params.get('stage', '0') in '01234' :
        assign = os.path.join(params['workspace'], 'read_assignment.gz')
        profiling(data, assign, **params)
    
    import glob
    for fname in glob.glob(os.path.join(params['workspace'], 'r?.fastq')) :
        subprocess.Popen(['gzip', fname]).communicate()
        #os.unlink(fname)
