# sigma integrated
import sys, subprocess, os, gzip, re, shutil, numpy as np, pandas as pd
import utils, time, glob
from multiprocessing import Pool
from scipy.stats import norm
from scipy.special import factorial

def parse_mapping(database, seqinfo, maps, mismatch=0.05, HGT_prior=[[0.05, 0.99, 0.1], [0.02, 0.99, 0.2], [0.01, 0.99, 0.5]], UCE_prior=[487, 2000], **args) :
    ln_fact = np.log(factorial(np.arange(1707, dtype=float)/10.))
    def dpois(freq, n) :
        res = np.zeros(n.shape)
        f1, f2, n = np.array(freq)/np.sqrt(2), np.array(freq)*np.sqrt(2), np.array(n)
        
        x = ((n - np.ceil(f2)) >= 20*np.sqrt(np.ceil(f2)))
        res[x], f2[x] = -230.2585, -1
        f2[n<=np.ceil(f2)], f1[n>=np.floor(f1)] = -1, -1

        x = (f2 > 0) & (f2 <= 20)
        res[x], f2[x] = ( np.log(f2[x])*n[x] - ln_fact[(n[x]*10).astype(int)] ) - \
                                ( np.log(f2[x])*np.ceil(f2[x]) - ln_fact[np.ceil(f2[x]).astype(int)*10] ), -1
        res[f2 > 0] = dnorm(f2[f2 > 0], n[f2 > 0]) - dnorm(f2[f2>0], np.ceil(f2[f2>0]))
        
        x = (f1 > 0) & (f1 <= 170)
        res[x], f1[x] = ( np.log(f1[x])*n[x] - ln_fact[(n[x]*10).astype(int)] ) - \
                                ( np.log(f1[x])*np.floor(f1[x]) - ln_fact[np.floor(f1[x]).astype(int)*10] ), -1
        res[f1 > 0] = dnorm(f1[f1 > 0], n[f1 > 0]) - dnorm(f1[f1>0], np.floor(f1[f1>0]))
        res[res>0] = 0.0
        return res
    def dnorm(freq, n) :
        if len(freq) == 0 :
            return np.array([])
        res = norm.cdf((freq-n+0.5)/np.sqrt(freq)) - norm.cdf((freq-n-0.5)/np.sqrt(freq))
        res[res < 1e-100] = -230.2585
        res[res > 0] = np.log(res[res>0])
        return res

    # 0. read ID
    # 1. taxa ID
    # 2. seq ID
    # 3. align Site
    # 4. No. mutation
    # 5. No. conserved
    # 6. similarity score
    # 7.  summarized score **
    # 8.  sparsity score
    # 9.  conservation score
    # 10. normalized align score
    
    # merge PE alignments
    n_tax = np.max(maps.T[1])+1
    tag =maps.T[0]*n_tax + maps.T[1]
    map_idx, tag_idx = np.unique(tag, return_index=True, return_inverse=True)[1:]
    m2 = np.zeros(shape=[map_idx.shape[0], 11])
    m2[:, :6] = maps[map_idx]
    m2[:, 4], m2[:, 5] = np.bincount(tag_idx, maps.T[4]), np.bincount(tag_idx, maps.T[5])
    
    # penalty for SE alignments
    m2[:, 6] = m2[:, 5] + m2[:, 4]
    m2 = m2[np.argsort(-m2.T[6])]
    map_idx, tag_idx = np.unique(m2.T[0], return_index=True, return_inverse=True)[1:]
    m2[:, 6] = m2[map_idx[tag_idx], 6] - m2[:, 6]
    m2[:, 4] += 0.1*m2[:, 6]
    m2[:, 5] += 0.9*m2[:, 6]
    
    maps = m2
    
    maps[:, 6] = ( maps[:, 4]*np.log(mismatch) + maps[:, 5]*np.log(1-mismatch) ) - \
                   np.sum(maps[:, 4:6], 1) * (mismatch*np.log(mismatch)+(1-mismatch)*np.log(1-mismatch))
    
    print time.time()-o_t, "Got mappings."
    sys.stdout.flush()

    # # similarity based weighting
    maps = maps[np.argsort(maps.T[6])[::-1]]
    read_ids = np.unique(maps.T[0], return_index=True, return_inverse=True)
    maps[:, 10] = np.exp(maps[:, 6] - maps[read_ids[1][read_ids[2]], 6])

    # # prepare cluster information from source DB
    b2 =       np.zeros(shape=[np.max(database['index'].astype(int))+1, len(next(database.itertuples()).barcode.split('.'))], dtype=int)
    barcodes = np.zeros(shape=[np.max(database['index'].astype(int))+1, len(next(database.itertuples()).barcode.split('.'))], dtype=int)
    barcodes.fill(-1)
    b2[database['index'].as_matrix().astype(int)] = database['barcode'].apply(lambda barcode:[int(b[1:]) for b in barcode.split('.')]).tolist()
    tax_ids = np.unique(maps.T[1]).astype(int)
    barcodes[tax_ids] = b2[tax_ids]
    
    # # PER match conservation weighting
    rcnt = np.max(maps.T[0]).astype(int)+1
    conservation = np.ones(shape=rcnt, dtype=float)
    for id, (prior, c_p, a_p) in enumerate(HGT_prior) :
        taxa_cnt = barcodes[np.unique(barcodes.T[id+3], return_index=True)[1], id]
        taxa_cnt = np.bincount(*np.unique(taxa_cnt[taxa_cnt >= 0], return_counts=True))

        if np.max(taxa_cnt) == 1 : break
        
        # (readID + taxID) united key
        maps[:, 7] = barcodes[maps.T[1].astype(int), id+3] + maps.T[0].astype(int)*(taxa_cnt.size + 1)
        maps[:, 8] = barcodes[maps.T[1].astype(int), id] + maps.T[0].astype(int)*(taxa_cnt.size + 1)
        c1 = np.unique(maps.T[7], return_index=True)[1]
        c2 = np.unique(maps[c1, 8], return_index=True, return_counts=True)
        c = [[], c1[c2[1]], c2[2]]
        
        maps[c[1], 8] = taxa_cnt[barcodes[maps[c[1], 1].astype(int), id]]
        maps[c[1], 9] = c[2]
        r_weight = np.bincount( maps[c[1], 0].astype(int), weights=maps[c[1], 10], minlength= rcnt)
        r_weight[r_weight == 0] = 1
        r_present = np.bincount( maps[c[1], 0].astype(int), weights=maps[c[1], 9]*maps[c[1], 10], minlength= rcnt)/r_weight
        r_total = np.bincount( maps[c[1], 0].astype(int), weights=maps[c[1], 8]*maps[c[1], 10], minlength= rcnt )/r_weight
        #p = (r_total > 0)
        lk1 = np.log(prior) + r_present*np.log(c_p) + (r_total-r_present)*np.log(1-c_p)
        lk1[lk1<-700] = -700
        lk2 = np.log(1-prior) + r_present*np.log(a_p) + (r_total-r_present)*np.log(1-a_p)
        lk2[lk2<-700] = -700
        lk1, lk2 = np.exp(lk1), np.exp(lk2)
        r_total = lk1/(lk1 + lk2)
        conservation = conservation*(1-r_total)
        
    maps[:, 9] = 1-conservation[maps[:,0].astype(int)]
    print time.time()-o_t, "conservation score."
    sys.stdout.flush()
    
    # remove a non-specific match if it is 5 SNP more divergent than the best match
    maps = maps[ maps.T[10] >= np.power(mismatch/(1-mismatch), 5) ]
    print time.time()-o_t, "remove highly divergent hits."
    sys.stdout.flush()
    
    # # PER match depth weighting
    # encode every 0.5 / 2 kb 
    seqinfo = seqinfo[np.in1d(seqinfo.T[2], maps.T[1])]
    seqlen = np.bincount(seqinfo.T[0], seqinfo.T[1])
    maps.T[8] = 1
    maps[seqlen[maps.T[2].astype(int)] < 500, 8] = 0
    for fid, fragment_len in enumerate(UCE_prior) :
        cov_score = np.zeros(shape=maps.shape[0], dtype=float)
        seq_info = np.zeros(shape=[np.max(seqinfo.T[0])+1, 5], dtype=int)
        seq_info[seqinfo.T[0], 0] = seqinfo.T[1]
        seq_info[seqinfo.T[0], 4] = seqinfo.T[2]
        seq_info.T[1] = 1+(seq_info.T[0]/fragment_len)
        seq_info.T[2, 0], seq_info.T[2, 1:] = 1, np.cumsum(seq_info.T[1,:-1]+1) + 1
        seq_info.T[3] = np.sum(seq_info[:, 1:3], 1)
        #
        # get depth for every region
        regions = np.zeros(shape=seq_info[-1, 3]+1)
        cov_score = seq_info[maps.T[2].astype(int), 2] + (maps.T[3]/fragment_len).astype(int)
        regions[:(np.max(cov_score).astype(int)+1)] = np.bincount(cov_score.astype(int) )#, weights=maps.T[10])
        # -- correct for continuity
        
        cov_score = regions[cov_score.astype(int)]
        print time.time()-o_t, "Start to run depth weighting."
        sys.stdout.flush()
        
        # tax_list = [[tax_id, region_start, region_end]]
        tax_ids = np.unique(seq_info.T[4], return_index=True)
        tax_list = np.vstack([tax_ids[0],tax_ids[1]])[:, np.argsort(tax_ids[1])].astype(float)
        tax_list = np.vstack([tax_list, np.concatenate([tax_list[1][1:], [seq_info.shape[0]]])  ]).T
        
        # get mean depth
        for tid, (tax_id, s, e) in enumerate(tax_list.astype(int)) :
            reg = regions[seq_info[s][2]:seq_info[e-1][3]] # np.concatenate([ np.arange(seq_info[ss][2], seq_info[ss][3]) for ss in np.arange(s, e) ])
            mean_cov = .5
            if np.max(reg) > 1 :
                mean_cov = np.mean(reg[reg>=0])
            tax_list[tid][1] = mean_cov
            
        mean_covs = np.bincount(tax_list.T[0].astype(int), tax_list.T[1])
        #default_lks = np.bincount(tax_list.T[0].astype(int), tax_list.T[2])
        cov_score = dpois(mean_covs[maps.T[1].astype(int)], cov_score)
        maps[maps.T[8] > 0, 8] = np.exp(cov_score[maps.T[8] > 0]) * maps[maps.T[8] > 0, 8]
    
        # # generalized weights for mapping depth
    map_dense = np.bincount(maps.T[0].astype(int), weights=maps.T[10]*np.power(maps.T[8], 1.0/len(UCE_prior)))
    #map_dense = np.bincount(maps.T[0].astype(int), weights=maps.T[10]*maps.T[8])
    map_score = np.bincount(maps.T[0].astype(int), weights=maps.T[10])
    map_dense[map_score >0] = map_dense[map_score > 0]/map_score[map_score>0]
    maps.T[8] = map_dense[maps.T[0].astype(int)]
    
    # # summarise maps.T[6] (similarity weighting),    maps.T[8] (depth weighting),    maps.T[9] (conservation weighting)
    maps.T[6] = np.exp(maps.T[6])
    maps.T[7] = maps.T[8]*maps.T[9]
    maps = maps[(maps.T[7] > 0) & (maps.T[6] > 0)]
    print time.time()-o_t, "Get all weighting."
    sys.stdout.flush()

    return maps

    
def summary_matrix(data, MapDB, workspace, bowtie_db, mismatch=0.05, n_thread=10, **params) :
    db_list = searchDBs(bowtie_db, MapDB)

    # # prepare taxonomy info for sequences
    seq_files, map_files = [], []
    for dbname in db_list :
        seq_files.append(os.path.join(workspace, dbname.rsplit('/', 1)[-1] + '.seqinfo.npy'))
        map_files.extend(glob.glob(os.path.join(workspace, dbname.rsplit('/', 1)[-1] + '.map.*.npy')))
    print time.time()-o_t, "info prepared."
    sys.stdout.flush()
    
    maps, seqinfo = np.vstack([np.load(mfile) for mfile in map_files]), np.vstack([np.load(sfile) for sfile in seq_files])

    maps = parse_mapping(data, seqinfo, maps, mismatch)
    
    for mid in xrange(0, maps.shape[0], 20000000) :
        maps[mid:(mid+20000000)].dump( os.path.join(workspace, 'read_scores.{0}.npy'.format(mid/20000000)) )
    
    return maps

def searchDBs(bowtie_db, MapDB) :
    db_list = []
    for mdb in MapDB.split(',') :
        master_file = os.path.join(bowtie_db, '{0}.info'.format(mdb))
        assert os.path.isfile(master_file), 'some database is not present'
        for fname in glob.glob(os.path.join(bowtie_db, '{0}.*.taxa.gz'.format(mdb))) :
            db_list.append(fname[:-8])
    return sorted(db_list)

def bowtie2matrix(bowtie_db, MapDB, workspace, mismatch, r1, r2=None, n_thread=10, **params) :
    def convert_file(source, target) :
        with open(target, 'w') as fout :
            fin = gzip.open(source) if r1.endswith('gz') else open(source)
            line = fin.readline()
            fin.close()
            fin = gzip.open(source) if r1.endswith('gz') else open(source)
            if line.startswith('>') :
                line_id, seq = 0, []
                for line in fin :
                    if line.startswith('>') :
                        if len(seq) > 0 :
                            seq = ''.join(seq)
                            fout.write('{0}\n+\n{1}\n'.format(seq, 'H' * len(seq)))
                        fout.write('@{0}\n'.format(line_id))
                        line_id, seq = line_id + 1, []
                    else :
                        seq.extend(line.strip().split())
                seq = ''.join(seq)
                fout.write('{0}\n+\n{1}\n'.format(seq, 'H' * len(seq)))
            else :
                for id, line in enumerate(fin) :
                    fout.write(line if id%4 else '@{0}\n'.format(id/4)) 
            fin.close()
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

    convert_file(r1, os.path.join(workspace, 'r1.fastq'))
    if r2 is None :
        read_info = '-U {0}'.format(os.path.join(workspace, 'r1.fastq'))
    else :
        read_info = '-1 {0} -2 {1}'.format(os.path.join(workspace, 'r1.fastq'), os.path.join(workspace, 'r2.fastq'))
        convert_file(r2, os.path.join(workspace, 'r2.fastq'))
            
    # run bowtie2 in multiple threads
    for db_prefix in db_list :
        out_prefix = os.path.join(workspace, db_prefix.rsplit('/', 1)[-1])
        if os.path.isfile(db_prefix+'.4.bt2') :
            cmd_bt2 = '{0} -k 300 --no-unal --ignore-quals --score-min L,20,1.1 --mp 4,4 --np 4 --sensitive-local -q -p {2} -I 25 -X 800 -x {3} {1}'.format(
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
                            seq_info.append([seqs[name][1],       # seq ID
                                             int(seq_len),        # seq Length
                                             seqs[name][0]] )     # tax ID
                else :
                    part = line.strip().split()
                    if part[2] in seqs :
                        key = (part[0], int(part[1]) & 128)
                        key2 = (part[0], 128-key[1])
                        if key not in prev:
                            if key2 not in prev :
                                prev = {key:{}}
                            else :
                                prev[key] = {}
                            
                        if seqs[part[2]][0] not in prev[key] :
                            s = [int(l) for l in re.findall('(\d+)S', part[5])]
                            if sum(s) > 100 or sum(s) * 2 >= len(part[9]) or (len(s) > 1 and min(s) > 4) :
                                continue
                            else :
                                prev[key][seqs[part[2]][0]] = 1
                        else :
                            continue
                        score = int(part[11][5:]) + (0 if len(s) == 0 else max(s)*1.5)
                        mut = (len(part[9])*2 - score)/6.0
                        nom = len(part[9]) - mut
                        res.append([ int(part[0]),        # read ID
                                     seqs[part[2]][0],    # taxa ID
                                     seqs[part[2]][1],    # seq ID
                                     int(part[3]),        # align Site
                                     mut,                 # No. mutation
                                     nom,                 # No. conserved
                                    ])
        else :
            for read in read_info.split()[1::2] :
                cmd_malt = '{0} -m BlastN -i {1} -a {4} -t {2} -mq 200 -e 0.00001 -id 70 -mem load -d {3}.malt'.format(
                    params['malt_run'], read, n_thread, db_prefix, out_prefix+'.tmp.gz')
                subprocess.Popen(cmd_malt.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                
                prev, res, seq_info = {}, [], []
                used_names = {}
                with gzip.open(out_prefix+'.tmp.gz') as fin :
                    for line in fin :
                        if line.startswith('@') :
                            if line.startswith('@SQ') :
                                name, seq_len = [ p[3:] for p in line.strip().split()[1:3] ]
                                name2 = name.split('|tax', 1)[0]
                                if name2 in seqs :
                                    seqs[name] = seqs.pop(name2)
                                    used_names[name] = 1
                                    seq_info.append([seqs[name][1],       # seq ID
                                                     int(seq_len),        # seq Length
                                                     seqs[name][0]] )     # tax ID
                        else :
                            part = line.strip().split()
                            if part[2] in seqs :
                                key = (part[0], int(part[1]) & 128)
                                key2 = (part[0], 128-key[1])
                                if key not in prev:
                                    if key2 not in prev :
                                        prev = {key:{}}
                                    else :
                                        prev[key] = {}
                                    
                                if seqs[part[2]][0] not in prev[key] :
                                    s = [int(l) for l in re.findall('(\d+)H', part[5])]
                                    if len(s) > 1 and min(s) > 9 :
                                        continue
                                    m = sum([int(l) for l in re.findall('(\d+)[MI]', part[5])])
                                    if sum(s) > m :
                                        continue
                                    prev[key][seqs[part[2]][0]] = 1
                                else :
                                    continue
                                rlen = sum(s) + m
                                score = int(part[14][5:]) + (0 if len(s) == 0 else max(s)*1.6)
                                mut = (rlen*2 - score)/5.0
                                nom = rlen - mut
                                res.append([ int(part[0]),        # read ID
                                             seqs[part[2]][0],    # taxa ID
                                             seqs[part[2]][1],    # seq ID
                                             int(part[3]),        # align Site
                                             mut,                 # No. mutation
                                             nom,                 # No. conserved
                                            ])
                os.unlink(out_prefix+'.tmp.gz')
        for name in used_names: 
            seqs.pop(name, None)
        sfile = out_prefix + '.seqinfo.npy'
        np.array(seq_info).dump(sfile)
        res = np.array(res)
        if res.size > 0 :
            for mid in xrange(0, res.shape[0], 20000000) :
                mfile = out_prefix + '.map.{0}.npy'.format(mid/20000000)
                res[mid:(mid+20000000)].dump( mfile )
    return True

def ipopt(workspace, least_amount=[0.000001, 5], bootstrap=100, **params) :
    def write_down(fout, reads, w) :
        if len(reads) > 0 :
            fout.write('*\t{0}\t{2}\t{1}\n'.format(int(reads.keys()[0]), '\t'.join(['{0}={1}'.format(t, p) for t, p in sorted(reads.values()[0].items())]), w))
    
    # # get total numbers of reads
    if os.path.isfile(os.path.join(workspace, 'r1.fastq')) :
        with open(os.path.join(workspace, 'r1.fastq')) as fin :
            fin.seek(-10000, 2)
            line = fin.readlines()[-4]
            n_read = int(line.strip()[1:])+1
    else :  
        with gzip.open(os.path.join(workspace, 'r1.fastq.gz')) as fin :
            for lid, line in enumerate(fin) :
                pass
            n_read = lid/4 +1 
    
    for ite in xrange(bootstrap+1) :
        maps = []
        for id in xrange(99999) :
            if os.path.isfile(os.path.join(workspace, 'read_scores.{0}.npy'.format(id))) :
                maps.append(np.load(os.path.join(workspace, 'read_scores.{0}.npy'.format(id))))
            else :
                break
        maps = np.vstack(maps)
        
        maps = maps[(maps.T[9] >= 0.01) & (maps.T[8] >= 0.01) & (maps.T[7]>= 0.001) ]
        if ite > 0 :
            maps = maps[ np.random.randint(0, maps.shape[0], maps.shape[0]) ]
    
        picked = []
        workon = np.copy(maps)
        r_limit = max(n_read*least_amount[0], least_amount[1])
        while True :
            tax_match = np.bincount(workon.T[1].astype(int), weights=workon.T[10]*workon.T[7])
            if tax_match.size == 0 : break
            tax_id = np.argmax(tax_match)
            if tax_match[tax_id]*2.5 < r_limit :
                break
            picked.append([tax_id, tax_match[tax_id]])
            workon = workon[(~np.in1d(workon.T[0], workon[(workon.T[1] == tax_id) & (workon.T[10] >= 0.5), 0])) & (workon.T[1] != tax_id)]
        picked = np.array(picked)
        if picked.shape[0] > 5000 :
            r_limit = max(n_read*least_amount[0]*5, least_amount[1]*2)
            picked = picked[ picked.T[1] >= r_limit ]
        
        print time.time()-o_t, "got candidates."
        sys.stdout.flush()
        
        maps = maps[np.in1d(maps.T[1], picked.T[0])]
        maps = maps[np.lexsort([-maps.T[10], maps.T[0]])]
        fname = os.path.join(workspace, 'ipopt.qmatrix') if ite == 0 else os.path.join(workspace, 'ipopt.bootstrap.{0}'.format(ite))
        with open(fname, 'w') as fout :
            fout.write('''#\t+\tMatrixName\tTotalNumberReads\tMatchedReads\tUnmatchedReads\tReadLimit
#\t@\tGenomeIndex\tGenomeName\tMatchedNumber
#\t*\tReadID\tWeight\tGenomeIndex=QValue
+\t{0}\t{1}\t{2}\t{3}\t{4}
'''.format(fname, n_read, (np.unique(maps.T[0])).size, n_read-(np.unique(maps.T[0])).size, r_limit))
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
        
        subprocess.Popen('{ipopt} -t {n_thread} -i {0}'.format(fname, **params).split()).communicate()
    return os.path.join(workspace, 'ipopt.qmatrix.solution')

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
    maps[maps.T[1]+1 >= taxa.shape[0], 1] = -1
    maps[taxa[maps.T[1].astype(int)+1] == 0, 1] = -1
    taxa[0] = r_limit/n_read

    # justified similarity score
    maps.T[10] = maps.T[6]*taxa[maps.T[1].astype(int)+1]
    # remove low level hits
    maps = maps[np.argsort(-maps.T[10])]
    r_max_score = np.unique(maps.T[0], return_index=True, return_inverse=True)
    
    maps.T[10] = (maps.T[10]+1e-301) / (maps[ r_max_score[1][r_max_score[2]] , 10] + 1e-300)
    maps = maps[maps.T[10] >= 0.1]

    read = np.bincount(maps.T[0].astype(int), maps.T[10])
    maps.T[10] = maps.T[10]/read[maps.T[0].astype(int)]
    
    taxa = np.bincount(maps.T[1].astype(int)+1, maps.T[10])
    maps[taxa[maps.T[1].astype(int)+1] < r_limit, 1] = -1
    maps = maps[np.lexsort([maps.T[4], -maps.T[10], maps.T[0]])]
    top_hits = np.unique(maps.T[0], return_index=True)[1]
    top_ref = maps[top_hits, 1][:]
    maps[maps.T[1] == -1, 1] = -2
    maps[top_hits, 1] = top_ref
    maps = maps[maps.T[1] >= -1]

    read = np.bincount(maps.T[0].astype(int), maps.T[10])
    maps.T[10] = maps.T[10]/read[maps.T[0].astype(int)]
    

    # 4: read distance;   5: mean distance
    #significant_aln = (maps.T[10] >= 0.05)
    aln_len = np.bincount(maps[:, 1].astype(int)+1, (maps.T[7]*maps.T[10]*(maps.T[4] + maps.T[5])))
    aln_mut = np.bincount(maps[:, 1].astype(int)+1, (maps.T[7]*maps.T[10]*maps.T[4]))
    aln_len[aln_len==0] = 1.0
    
    maps.T[4] = 100.0*maps.T[4]/(maps.T[4] + maps.T[5])
    maps.T[5] = 100.0*(aln_mut/aln_len)[maps.T[1].astype(int)+1]
    
    cls_map = np.bincount( data['index'].astype(int)+1, data['barcode'].apply(lambda barcode:int(barcode.split('.')[0][1:])).as_matrix() ).astype(int)
    cls_cnt = maps[:, [1,10,9]]
    cls_cnt = cls_cnt[cls_cnt.T[2]>0.001]
    cls_cnt.T[0] = cls_map[cls_cnt.T[0].astype(int)+1]
    n = np.bincount(cls_cnt.T[0].astype(int), cls_cnt.T[1])
    n[n == 0] = 1
    cls_cnt = np.bincount(cls_cnt.T[0].astype(int), cls_cnt.T[1]*cls_cnt.T[2])/n
    maps.T[9] /= cls_cnt[cls_map[maps.T[1].astype(int)+1]]
    n=np.bincount(maps.T[0].astype(int), maps.T[10])
    n[n==0] = 1
    cls_cnt = np.bincount(maps.T[0].astype(int), maps.T[9]*maps.T[10])/n
    maps.T[9] = cls_cnt[maps.T[0].astype(int)]
    
    read = {}
    for m in maps :
        if m[0] not in read :
            read[m[0]] = [m[8], m[9]]
        read[m[0]].append([m[1].astype(int), m[10], m[4], m[5]])
    with gzip.open(os.path.join(workspace, 'read_assignment.gz'), 'w') as fout :
        fout.write('#\tReadID\tReadName\tDepth_Weight\tCore_Weight\tReferenceID:Proportion:Distance:MeanDistance\n')
        fin = gzip.open(params['r1']) if params['r1'].endswith('.gz') else open(params['r1'])
        line = fin.readline()
        fin.close()
        fin = gzip.open(params['r1']) if params['r1'].endswith('.gz') else open(params['r1'])
        if line.startswith('>') :
            rid = 0
            for line in fin :
                if line.startswith('>') :
                    name = line.strip().split()[0][1:]
                    if rid in read :
                        fout.write('{0}\t{1}\t{2:.5f}\t{3:.5f}\t{4}\n'.format(rid, name, read[rid][0], read[rid][1], '\t'.join(['{0}:{1:.2f}:{2:.2f}:{3:.2f}'.format(*x) for x in read[rid][2:]])))
                    else :
                        fout.write('{0}\t{1}\n'.format(rid, name))
                    rid += 1
        else :
            for id, line in enumerate(fin) :
                if id % 4 == 0 :
                    rid = id / 4
                    if rid in read :
                        fout.write('{0}\t{1}\t{2:.5f}\t{3:.5f}\t{4}\n'.format(rid, line.strip().split()[0][1:], read[rid][0], read[rid][1], '\t'.join(['{0}:{1:.2f}:{2:.2f}:{3:.2f}'.format(*x) for x in read[rid][2:]])))
                    else :
                        fout.write('{0}\t{1}\n'.format(rid, line.strip().split()[0][1:]))
        fin.close()
    print 'Read-level assignments  are in {0}'.format(os.path.join(workspace, 'read_assignment.gz'))

def match_group(data, match_ref, **params) :
    pgroup = np.vstack([ data['index'], data['barcode'].apply(lambda b:b.split('.')[3] ) ]).T
    pop = np.unique(pgroup[ np.in1d(pgroup.T[0], [r for r in match_ref.keys()]) , 1])
    groups = {}
    for g in pop :
        groups[g] = [{'':0} for c in params['taxa_columns'] ] 
        ref = pgroup[pgroup.T[1]==g, 0]
        d = data.loc[data['index'].isin(ref)]
        for part in d[['index', 'barcode'] + list(reversed(params['taxa_columns']))].as_matrix() :
            index, barcode, taxa = part[0], part[1], part[2:]
            for t, gg in zip(taxa, groups[g]) :
                if not (t in ('', '-') and t.startswith('*')):
                    gg[t] = gg.get(t, 0) + 1
    results = {}
    for idx, assign in sorted( match_ref.iteritems(), key=lambda x:np.sum(x[1].T[1]), reverse=True ):
        level_abundance = [np.sum(assign[:, 1]), \
                           np.sum(assign[assign.T[0] <= 10., 1]), \
                           np.sum(assign[assign.T[0] <= 5., 1]), \
                           np.sum(assign[assign.T[0] <= 2., 1]), \
                           np.sum(assign[assign.T[0] <= 1., 1]), \
                        ]
        barcodes = data.loc[data['index'] == idx]['barcode'].as_matrix()[0].split('.')
        for grp, abnd in zip (['~' + barcodes[0][1:]] + barcodes[:4], level_abundance) :
            if abnd <= 0 : continue
            if grp not in results :
                results[grp] = [abnd, [idx], [ { label:cnt*abnd for label, cnt in lvl.iteritems() } for lvl in groups[barcodes[3]] ]]
            else :
                results[grp][0] +=abnd
                results[grp][1].append(idx)
                for res, lvl in zip(results[grp][2], groups[barcodes[3]]) :
                    for label, cnt in lvl.iteritems() :
                        res[label] = res.get(label, 0.) + cnt * abnd
    for label, content in results.iteritems() :
        taxonomy = []
        for level in content[2] :
            cnts = sorted(level.items(), key=lambda x:x[1], reverse=True)
            if len(cnts)> 1 and cnts[1][1] >= 0.1*cnts[0][1] :
                taxonomy.append( '{0} ({1})'.format(cnts[0][0], '/'.join([ c[0] for c in cnts[1:] if c[1] >= cnts[0][1]*0.1 ]) ) )
            else :
                taxonomy.append(cnts[0][0])
        if label.startswith('p') :
            taxonomy.append( '{0[0]}: {0[1]}'.format(data.loc[data['index'] == content[1][0], ['organism_name', 'assembly_accession']].as_matrix().tolist()[0]) )
        elif label.startswith('s') :
            taxonomy = taxonomy[:-1]
        elif label.startswith('u') :
            taxonomy = taxonomy[:-2]
        elif label.startswith('~') :
            taxonomy = taxonomy[:-3]
        content[2] = taxonomy
    return results

def profiling(data, assign, minFreq, **params) :
    basic_aln = [0, 0., 0.]
    match_ref = {}
    with gzip.open(assign) as fin :
        for line in fin :
            if line.startswith('#') : continue
            part = line.strip().split('\t')
            basic_aln[0] += 1
            if len(part) >= 3 :
                w = float(part[2]) * float(part[3])
                basic_aln[1] += 1
                for match in part[4:] :
                    ref, prop, d2, d1 = match.split(':')
                    if ref != '-1' :
                        basic_aln[2] +=w*float(prop)
                        d = max(float(d1), float(d2))
                        if ref not in match_ref :
                            match_ref[ref] = []
                        match_ref[ref].append([d, w*float(prop) ])
                        
    match_ref = { idx:np.array(reads) for idx, reads in match_ref.iteritems() }
    groups = match_group(data, match_ref, **params)    
    
    with open(os.path.join(params['workspace'], 'profile.txt'), 'w') as fout :
        fout.write('Total\t{0}\t{1}\n'.format(basic_aln[0], basic_aln[1]))
        fout.write('Unmatched\t{0:.3f}\t{1:.3f}\n'.format((basic_aln[0]-basic_aln[1])*100.0/basic_aln[0], 0.0))
        fout.write('Uncertain_match\t{0:.3f}\t{1:.3f}\n'.format((basic_aln[1] - basic_aln[2])*100.0/basic_aln[0], (basic_aln[1] - basic_aln[2])*100.0/basic_aln[1]))
        for g, (w, r, t) in sorted(groups.iteritems(), key=lambda x:(x[1][0], x[0]), reverse=True) :
            if w/basic_aln[0] >= minFreq :
                fout.write('{0}\t{1:.4f}\t{2:.4f}\t{3} ({4})\n'.format(g, w*100.0/basic_aln[0], w*100.0/basic_aln[1], '|'.join(t), ','.join(r) ))
    print 'Profilling results are in {0}'.format(os.path.join(params['workspace'], 'profile.txt'))
if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    params['bootstrap'] = int(params['bootstrap']) if 'bootstrap' in params else 0
        
    data = utils.load_database(**params)
    
    o_t = time.time()
    if params.get('stage', '0') in '0' :
        bowtie2matrix(**params)
    if params.get('stage', '0') in '01' :
        summary_matrix(data, **params)
    if params.get('stage', '0') in '012' :    
        qvector = ipopt(least_amount=[params['minFreq'], params['minNum']], **params)
    if params.get('stage', '0') in '0123' :
        qvector = os.path.join(params['workspace'], 'ipopt.qmatrix.solution')
        assign_reads(data, qvector, **params)
    if params.get('stage', '0') in '01234' :
        assign = os.path.join(params['workspace'], 'read_assignment.gz')
        profiling(data, assign, **params)
    
    import glob
    for fname in glob.glob(os.path.join(params['workspace'], 'r?.fastq')) :
        subprocess.Popen(['gzip', '-f', fname]).communicate()
        #os.unlink(fname)
