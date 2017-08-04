import sys, os, utils, numpy as np, pandas as pd, subprocess, shutil, multiprocessing

lastal = '/home/zhemin/bin/lastal'
lastdb = '/home/zhemin/bin/lastdb'


def make_alignment( inputs ) :
    def encode(seq) :
        code = {'A':1, 'C':2, 'G':3, 'T':4, '-':5}
        s = np.array(list(seq))
        bases, indices = np.unique(s, return_inverse=True)
        bases = np.array([code.get(b, 0) for b in bases])
        return bases[indices]

    def call_mutation( comparison, rep_mask = 3 ) :
        seq = np.vstack([np.zeros(comparison[6].size), np.zeros(comparison[6].size), comparison[6], comparison[12], comparison[13]]).astype(int).T
        
        v = np.diff(np.concatenate([[0], seq.T[4], [0]]))
        low_qual = np.vstack([ np.where(v>0), np.where(v<0) ]).T
        high_qual = np.vstack([ np.concatenate([[0], low_qual.T[1]]), np.concatenate([low_qual.T[0], [seq.shape[0]]]) ]).T
        for ss, se in high_qual[high_qual.T[1] - high_qual.T[0] < 8] :
            seq[ss:se, 4] = 1
        seq[low_qual[low_qual.T[0] > 0, 0]-1 ,4] = 1
        seq[low_qual[low_qual.T[1] < seq.shape[0], 1] ,4] = 1
        
        direct = [1 if comparison[4] == '+' else -1, 1 if comparison[10] == '+' else -1]
        seq.T[0] = np.arange(comparison[2] *direct[0], comparison[2] *direct[0] + seq.shape[0], 1)
        for site in np.where(seq.T[2] > 4)[0] :
            seq[site:, 0] -= 1
        
        seq.T[1] = np.arange(comparison[8] *direct[1], comparison[8] *direct[1] + seq.shape[0], 1)
        for site in np.where(seq.T[3] > 4)[0] :
            seq[site:, 1] -= 1
        
        v = np.diff(np.concatenate([[0], seq.T[4], [0]]))
        low_qual = np.vstack([ seq[np.where(v>0), 0], seq[np.where(v<0)[0]-1, 0], seq[np.where(v>0), 1], seq[np.where(v<0)[0]-1, 1] ]).T
        
        v = np.diff(np.concatenate([[0], seq.T[3] > 4, [0]]))
        deletes = np.vstack([seq[np.where(v>0), 0], seq[np.where(v<0)[0]-1, 0], seq[np.where(v>0), 1], seq[np.where(v<0)[0]-1, 1]])
        deletes = np.vstack([deletes, np.zeros(shape=[3, deletes.shape[1]])]).T
        deletes.T[5] = 5
        deletes.T[6] = 1
        v = np.diff(np.concatenate([[0], seq.T[2] > 4, [0]]))
        inserts = np.vstack([seq[np.where(v>0), 0], seq[np.where(v<0)[0]-1, 0], seq[np.where(v>0), 1], seq[np.where(v<0)[0]-1, 1]])
        inserts = np.vstack([inserts, np.zeros(shape=[3, inserts.shape[1]])]).T
        inserts.T[4] = 5
        inserts.T[6] = 1
        
        mutations = seq[(seq.T[2] != seq.T[3]) & (seq.T[4] == 0)]
        comparison[0] = (comparison[3]-comparison[2]+1) - 3 * mutations.shape[0] - 7*(deletes.shape[0]+inserts.shape[0]) - 2*np.sum(deletes.T[1] - deletes.T[0]+1) - np.sum(inserts.T[1] - inserts.T[0]+1)
        
        mutations = np.vstack([np.vstack([mutations.T[0],mutations.T[0],mutations.T[1],mutations.T[1],mutations.T[2],mutations.T[3],mutations.T[4]]).T, deletes, inserts])

        comparison[6] = low_qual.astype(int).tolist()
        comparison[12] = ''
        comparison = comparison[:13] + mutations[np.argsort(mutations.T[0])].astype(int).tolist()

        return comparison

    def sub_comparison(comparison, ref_coords=None, qry_coords=None) :
        if (ref_coords is None) == (qry_coords is None) :
            raise Exception('Exact only one set of coords, either qry or ref')
        direct = [1 if comparison[4] == '+' else -1, 1 if comparison[10] == '+' else -1]
        mutations = []
        if ref_coords is not None :
            if direct[0] > 0 :
                rc = [max(ref_coords[0], comparison[2]), min(ref_coords[1], comparison[3])]
            else :
                rc = [max(-ref_coords[1], -comparison[2]), min(-ref_coords[0], -comparison[3])]
            if rc[0] > rc[1] :
                return []
            qc = [0, 0]
            mutations = [ mut[:] for mut in comparison[13:] if mut[1]>= rc[0] and mut[0] <= rc[1] ]
            if len(mutations) > 0 :
                if mutations[0][0] <= rc[0] and mutations[0][5] > 4 :
                    qc[0] = mutations[0][2] + 1
                    mutations[0][0] = rc[0]
                elif mutations[0][4] > 4 :
                    qc[0] = mutations[0][2] - (mutations[0][0]-rc[0]+1)
                elif mutations[0][5] > 4 :
                    qc[0] = mutations[0][2] - (mutations[0][0]-rc[0]-1)
                else :
                    qc[0] = mutations[0][2] - (mutations[0][0]-rc[0])
                if mutations[-1][0] == rc[1] and mutations[-1][4] > 4 :
                    qc[1] = mutations[-1][2]-1
                    mutations.pop(-1)
                else :
                    qc[1] = mutations[-1][3] + (rc[1] - mutations[-1][1])
                    if mutations[-1][1] >= rc[1] and mutations[-1][5] > 4 :
                        mutations[-1][1] = rc[1]
            else :
                pre_mut = [ mut for mut in [[comparison[2]*direct[0], comparison[2], comparison[8]*direct[1], comparison[8]]] + comparison[13:] if mut[0] <= rc[0] ][-1]
                delta = -pre_mut[0] + pre_mut[2]
                qc = [rc[0]+delta, rc[1] + delta ]
        else :
            if direct[1] > 0 :
                qc = [max(qry_coords[0], comparison[8]), min(qry_coords[1], comparison[9])]
            else :
                qc = [max(-qry_coords[1], -comparison[8]), min(-qry_coords[0], -comparison[9])]
            if qc[0] > qc[1] :
                return []
            rc = [0, 0]
            mutations = [ mut[:] for mut in comparison[13:] if mut[3]>= qc[0] and mut[2] <= qc[1] ]
            if len(mutations) > 0 :
                if mutations[0][2] <= qc[0] and mutations[0][4] > 4 :
                    rc[0] = mutations[0][0] + 1
                    mutations[0][2] = qc[0]
                elif mutations[0][4] > 4 :
                    rc[0] = mutations[0][0] - (mutations[0][2]-qc[0]-1)
                elif mutations[0][5] > 4 :
                    rc[0] = mutations[0][0] - (mutations[0][2]-qc[0]+1)
                else :
                    rc[0] = mutations[0][0] - (mutations[0][2]-qc[0])
                if mutations[-1][2] >= qc[1] and mutations[-1][5] > 4 :
                    rc[1] = mutations[-1][0]-1
                    mutations.pop(-1)
                else :
                    rc[1] = mutations[-1][1] + (qc[1]-mutations[-1][3])
                    if mutations[-1][3] >= qc[1] and mutations[-1][4] > 4 :
                        mutations[-1][3] = qc[1]
            else :
                pre_mut = [ mut for mut in [[comparison[2]*direct[0], comparison[2], comparison[8]*direct[1], comparison[8]]] + comparison[13:] if mut[2] <= qc[0] ][-1]
                delta = -pre_mut[0] + pre_mut[2]
                rc = [qc[0]-delta, qc[1] - delta ]

        ref_ins = [mut[3]-mut[2]+1 for mut in mutations if mut[5] > 4]
        qry_ins = [mut[1]-mut[0]+1 for mut in mutations if mut[4] > 4]

        gap_open = len(ref_ins) + len(qry_ins)
        mut_num = len(mutations) - gap_open
        score = (rc[1]-rc[0]+1) - 3 * mut_num - 7*gap_open - 2*sum(ref_ins) - sum(qry_ins)
        return [score, comparison[1], abs(rc[0]), abs(rc[1]), comparison[4], comparison[5], comparison[6], comparison[7], abs(qc[0]), abs(qc[1]), comparison[10], comparison[11], comparison[12]] + mutations


    comparisons = []
    for line in inputs:
        if line[0] == 'a' :
            comparison = [ int(line.split(' ', 2)[1][6:]) ]
        elif line[0] == 's' :
            part = line.strip().split()[1:]
            part[1:5] = [int(part[1]), int(part[2]), part[3], int(part[4])]
            if part[3] == '+' :
                part[1:3] = [part[1]+1, part[1]+part[2]]
            else :
                part[1:3] = [part[4]-part[1], part[4]-part[1]-part[2]+1]
            comparison.extend(part)
            if len(comparison) >= 13 :
                comparison[6], comparison[12] = encode(comparison[6]), encode(comparison[12])
                comparison.append((comparison[6] == 0) | (comparison[12] == 0) | (comparison[6] == 5) | (comparison[12] == 5))
        elif line[0] in 'pq' :
            part = line.strip().split()
            q = np.array(list(part[-1]))
            comparison[13] = (comparison[13] | (q <= '*'))
        elif len(line.strip()) == 0 :
            if comparison[0] >= 100 and (comparison[1] != comparison[7] or comparison[4] != '+' or comparison[2] != comparison[8] or comparison[3] != comparison[9]) : 
                comparisons.append( call_mutation(comparison) )
    
    # remove significant low identity regions in query
    comparisons.sort(key=lambda x: min(x[8:10]) )
    comparisons.sort(key=lambda x: x[7] )

    low_q = []
    for id, regi in enumerate(comparisons) :
        if len(regi) == 0 : continue
        for jd in xrange(id+1, len(comparisons)) :
            regj = comparisons[jd]
            if len(regj) == 0 : continue
            if regi[7] != regj[7] : break
            si, ei = sorted(regi[8:10])
            sj, ej = sorted(regj[8:10])
            s = max(si, sj)
            e = min(ei, ej)
            if e >= s :
                overlap_i = sub_comparison(regi, qry_coords=[s, e])
                overlap_j = sub_comparison(regj, qry_coords=[s, e])
                
                if overlap_i[0] < 0.95 * overlap_j[0] and ( regi[0] < regj[0] or ei < ej ) : 
                    if s - si >= 30 :
                        comparisons[id] = sub_comparison(regi, qry_coords=[si, s-1])
                        if overlap_i[3] < overlap_i[2] and overlap_i[4] == '+' :
                            print comparisons[id]
                        overlap_i[12] = 'E'
                        low_q.append(overlap_i)
                        if overlap_i[3] >= overlap_i[2]  :
                            regi = comparisons[id]
                        if len(regi) == 0: break
                    else :
                        comparisons[id][12] = 'E'
                        break
                elif overlap_i[0] * 0.95 > overlap_j[0] :
                    if ej - e >= 30 :
                        comparisons[jd] = sub_comparison(regj, qry_coords=[e+1, ej])
                        overlap_j[12] = 'E'
                        if overlap_j[3] >= overlap_j[2]  :
                            low_q.append(overlap_j)
                    else :
                        comparisons[jd][12] = 'E'
                elif s == si and e == ei and regj[0] > regi[0]*3 and overlap_i[0] <= overlap_j[0] :
                    comparisons[id][12] = 'E'
                    break
                elif s == sj and e == ej and regi[0] > regj[0]*3 and overlap_i[0] >= overlap_j[0] :
                    comparisons[jd][12] = 'E'
                else :
                    comparisons[id][12] = 'D'
                    comparisons[jd][12] = 'D'
            else :
                break

    # remove significant low identity regions in reference
    comparisons = sorted([x for x in comparisons if len(x) > 0] + low_q, key=lambda x: x[2] )
    comparisons.sort(key=lambda x: x[1] )

    for id, regi in enumerate(comparisons) :
        if len(regi) == 0 : continue
        for jd in xrange(id+1, len(comparisons)) :
            regj = comparisons[jd]
            if len(regj) == 0 : continue
            if regi[1] != regj[1] : break
            si, ei = regi[2:4]
            sj, ej = regj[2:4]
            s = max(si, sj)
            e = min(ei, ej)
            if e >= s :
                overlap_i = sub_comparison(regi, ref_coords=[s, e])
                overlap_j = sub_comparison(regj, ref_coords=[s, e])
                
                if overlap_i[0] < 0.95 * overlap_j[0] and ( regi[0] < regj[0] or ei < ej ) : 
                    if s - si >= 30 :
                        comparisons[id] = sub_comparison(regi, ref_coords=[si, s-1])
                        regi = comparisons[id]
                        if len(regi) == 0: break
                    else :
                        comparisons[id] = []
                        break
                elif overlap_i[0] * 0.95 > overlap_j[0] :
                    if ej - e >= 30 :
                        comparisons[jd] = sub_comparison(regj, ref_coords=[e+1, ej])
                    else :
                        comparisons[jd] = []
                elif overlap_i[0] == overlap_j[0] and len(overlap_i) == len(overlap_j) :
                    if si == sj and ei == ej:
                        diff = 0
                        for i, i_snp in enumerate(overlap_i[13:]) :
                            j_snp = overlap_j[13+i]
                            if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                diff = 1
                                break
                        if diff == 0 :
                            if comparisons[id][12] in 'DE':
                                comparisons[id] = [] 
                                break
                            else: 
                                comparisons[jd] = []
                    elif si <= sj and ei >= ej :
                        diff = 0
                        for i, i_snp in enumerate(overlap_i[13:]) :
                            j_snp = overlap_j[13+i]
                            if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                diff = 1
                                break
                        if diff == 0 :
                            comparisons[jd] = []
                    elif si >= sj and ei <= ej:
                        diff = 0
                        for i, i_snp in enumerate(overlap_i[13:]) :
                            j_snp = overlap_j[13+i]
                            if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                diff = 1
                                break
                        if diff == 0 :
                            comparisons[id] = []
                            break
            else :
                break
        if len(comparisons[id]) > 0 :
            regi = comparisons[id]
            regi[6] = [lq for lq in regi[6] if lq[1] >= regi[2] and lq[0] <= regi[3]]
    
    # mark repetitive regions in query
    repeats = []
    mutations = {}
    comparisons = sorted([x for x in comparisons if len(x) > 0 and x[12] != 'E'], key=lambda x: min(x[8:10]) )
    comparisons = [c for c in comparisons if float(c[0])/(abs(c[3]-c[2])+1) >= 0.7]
    comparisons.sort(key=lambda x: x[7] )
    
    for id, regi in enumerate(comparisons) :
        for jd in xrange(id+1, len(comparisons)) :
            regj = comparisons[jd]
            if regi[7] != regj[7] : break
            si, ei = sorted(regi[8:10])
            sj, ej = sorted(regj[8:10])
            s = max(si, sj)
            e = min(ei, ej)
            if e >= s :
                for mut in regi[13:] :
                    if abs(min(mut[2:4])) <= e and abs(max(mut[2:4])) >= s :
                        mut[6] = 1
                for mut in regj[13:] :
                    if abs(min(mut[2:4])) <= e and abs(max(mut[2:4])) >= s :
                        mut[6] = 1
                overlap_i = sub_comparison(regi, qry_coords=[s, e])
                overlap_j = sub_comparison(regj, qry_coords=[s, e])
                regi[6] = [lq for lq in regi[6] if lq[0] <overlap_i[2] or lq[1] > overlap_i[3]]
                regj[6] = [lq for lq in regj[6] if lq[0] <overlap_j[2] or lq[1] > overlap_j[3]]
                repeats.append(overlap_i[1:4] + [0])
                repeats.append(overlap_j[1:4] + [0])
        
    # identify repetitive regions in the reference
    comparisons.sort(key=lambda x: x[2] )
    comparisons.sort(key=lambda x: x[1] )

    for id, regi in enumerate(comparisons) :
        if len(regi) == 0 : continue
        for jd in xrange(id+1, len(comparisons)) :
            regj = comparisons[jd]
            if regi[1] != regj[1] : break
            si, ei = sorted(regi[2:4])
            sj, ej = sorted(regj[2:4])
            s = max(si, sj)
            e = min(ei, ej)
            if e >= s :
                overlap_i = sub_comparison(regi, ref_coords=[s, e])
                overlap_j = sub_comparison(regj, ref_coords=[s, e])
                if len(overlap_i) == len(overlap_j) :
                    diff = 0
                    for i, i_snp in enumerate(overlap_i[13:]) :
                        j_snp = overlap_j[13+i]
                        if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                            diff = 1
                            break
                    if diff == 1 :
                        for mut in regi[13:] :
                            if abs(mut[0]) <= e and abs(mut[1]) >= s :
                                mut[6] = 1
                        for mut in regj[13:] :
                            if abs(mut[0]) <= e and abs(mut[1]) >= s :
                                mut[6] = 1
                        regi[6] = [lq for lq in regi[6] if lq[0] <s or lq[1] > e]
                        regj[6] = [lq for lq in regj[6] if lq[0] <s or lq[1] > e]
                        repeats.append([regi[1], s, e, 0])
                else :
                    for mut in regi[13:] :
                        if abs(mut[0]) <= e and abs(mut[1]) >= s :
                            mut[6] = 1
                    for mut in regj[13:] :
                        if abs(mut[0]) <= e and abs(mut[1]) >= s :
                            mut[6] = 1
                    regi[6] = [lq for lq in regi[6] if lq[0] <s or lq[1] > e]
                    regj[6] = [lq for lq in regj[6] if lq[0] <s or lq[1] > e]
                    repeats.append([regi[1], s, e, 0])
        for mut in regi[13:] :
            if mut[6] == 0 :
                if (regi[1], mut[0]) not in mutations : 
                    mutations[ (regi[1], mut[0]) ] = [[regi[1], regi[7]] + mut]
                else :
                    mutations[ (regi[1], mut[0]) ].append([[regi[1], regi[7]] + mut])
        repeats.extend([[regi[1]]+ lq[:2] + [1] for lq in regi[6]])
    
    repeats.sort(key=lambda x:x[1])
    repeats.sort(key=lambda x:x[0])
    repetitive_regions = []
    for rep in repeats:
        if len(repetitive_regions) == 0 or repetitive_regions[-1][0] != rep[0] or repetitive_regions[-1][2]+1 < rep[1] :
            repetitive_regions.append(rep)
        elif rep[2] > repetitive_regions[-1][2] :
            repetitive_regions[-1][2] = rep[2]
            if repetitive_regions[-1][3] > 0 :
                repetitive_regions[-1][3] = rep[3]

    return comparisons, repetitive_regions, mutations


def prepare_lastdb( reference, refdb = None ) :
    if refdb == None :
        refdb = reference
    if os.path.isfile(refdb + '.ssp') :
        return refdb
    cmd = '{0} -cR01 {2} {1}'.format( lastdb, reference, refdb )
    db_run = subprocess.Popen( cmd.split() )
    db_run.communicate()
    if db_run.returncode == 0 :
        return refdb
    else :
        return None

def run_lastal( data ) :
    refdb, (qry, query) = data
    output = query + '.lastal'
    out_seq = make_alignment(open(output).readlines())
    return out_seq
    
    cmd = '{0} -j4 -r1 -q2 -a7 -b1 {1} {2}'.format( lastal, refdb, query )
    lastal_run = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    inputs = lastal_run.communicate()[0]
    if lastal_run.returncode != 0 :
        cmd = '{0} -Q1 -j4 -r1 -q2 -a7 -b1 {1} {2}'.format( lastal, refdb, query )
        lastal_run = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE )
        inputs = lastal_run.communicate()[0]
    out_seq = make_alignment(inputs.split('\n'))
    return out_seq

def ref_align(ref, qry) :
    pool = multiprocessing.Pool(10)
    #refdb = prepare_lastdb(ref[1])
    #aligns = pool.map(run_lastal, [[refdb, q] for q in [ref] + qry])
    
    aligns = pool.map(run_lastal, [['', q] for q in [ref] + qry])
    #import pickle
    #pickle.dump(aligns, open('test', 'w'))
    
    #aligns = pickle.load(open('test'))
    
    # get core genome
    snp_sites = {(name, site):[ref_base[0][6]] for region, repeat, mut in aligns[1:] for (name, site), ref_base in mut.iteritems()}
    core_genome = {reg[1]:reg[5] for region, repeat, mut in aligns[1:] for reg in region}
    core_seq = {name:np.zeros(length, dtype=int) for name, length in core_genome.iteritems()}
    for region, repeat, mut in aligns[1:] :
        aln_seq = {name:np.zeros(length, dtype=int) for name, length in core_genome.iteritems()}
        for reg in region :
            aln_seq[reg[1]][reg[2]-1:reg[3]] = 1
        for rep in repeat :
            aln_seq[rep[0]][rep[1]-1:rep[2]] = 0
        for name, seq in aln_seq.iteritems() :
            core_seq[name] += seq
        for name_site in snp_sites :
            if name_site in mut :
                snp_sites[name_site].append(mut[name_site][0][7])
            elif aln_seq[name_site[0]][name_site[1]-1] == 1 :
                snp_sites[name_site].append(snp_sites[name_site][0])
            else :
                snp_sites[name_site].append('-')
    #for reg in aligns[0][0] :
        #core_seq[reg[1]][reg[2]-1:reg[3]] = 0
    for name, seq in core_seq.iteritems() :
        seq[:] = (seq >= 0.95*(len(aligns)-1))
    
    # generate matrix
    snp_sites = {(name, site):bases for (name, site), bases in snp_sites.iteritems() if core_seq[name][site-1] > 0}
    for (name, site), ref_base in snp_sites.iteritems() :
        core_seq[name][site-1] = 0
    constant_sites = {}
    for name, seq in core_seq.iteritems() :
        v = np.diff(np.concatenate([[0], seq, [0]]))
        constant_sites[name] = np.vstack([np.where(v > 0)[0] + 1, np.where(v < 0)[0]]).T
    return constant_sites, sorted(snp_sites.items())

def get_specific_reads(**params) :
    assert 'ref_id' in params and 'workspace' in params, 'missing ref_id or workspace'
    ref_ids = params['ref_id'].split(',')
    
    assignment, readseqs = os.path.join(workspace, 'read_assignment.gz'), sorted(glob.glob(os.path.join(workspace, 'r?.fastq.gz')))
    outputs = [os.path.join(workspace, '{0}.r{1}.fastq.gz'.format(ref_ids[0], id+1)) for id in xrange(len(readseqs))]
    
    ref_ids = set(ref_ids)
    read_ids = {}
    with gzip.open(assignment) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            if len(part) < 3 :
                continue
            ref_map = part[4].split(':')
            if float(part[2]) > 0 and ref_map[0] in ref_ids and float(ref_map[1]) >= params['ratio'] :
                read_ids[part[0]] = [part[1], part[4]]
    
    for readseq, output in zip(readseqs, outputs) :
        with gzip.open(readseq) as fin, gzip.open(output, 'w') as fout :
            for id, line in enumerate(fin) :
                if id%4 == 0 :
                    rid = line[1:].strip().split()[0]
                    if rid in read_ids :
                        write = 1
                        line = '@{0} {1} {2}\n'.format(read_ids[rid][0], rid, read_ids[rid][1])
                    else :
                        write = 0
                if write :
                    fout.write(line)
    return outputs

if __name__ == '__main__' :
    # get main reference
    params = utils.load_params(sys.argv)
    existing_data = os.path.join(params['dbname'], 'db_metadata.msg')
    data = pd.read_msgpack(existing_data)
    
    out_folder = os.path.join(params['dbname'], 'place_db', params['PlaceDB'])
    #if os.path.isdir(out_folder) :
        #shutil.rmtree(out_folder)
    #os.makedirs(out_folder)
    
    if 'sparse_ref' in params :
        if 'external_ref' in params :
            raise Exception('cannot have both sparse_ref and external_ref')
        ref_genome = utils.get_genome_file(params['sparse_ref'], data, out_folder)[0]
    elif 'external_ref' in params :
        ref_genome = params['external_ref']
    else :
        raise Exception('need to have either sparse_ref or external_ref')
    
    shutil.copy(ref_genome[1], os.path.join(out_folder, 'reference'))
    ref_genome[1] = os.path.join(out_folder, 'reference')
    
    # align other genomes onto reference
    qry_genome = []
    if 'sparse_qry' in params :
        sparse_qrys = params['sparse_qry'].split(',')
        qry_genome.extend(utils.get_genome_file(params['sparse_qry'], data, out_folder))
    if 'external_ref' in params :
        with open(params['external_qry']) as fin :
            for line in fin :
                qry_genome.append(line.strip().split()[:2])
    qry_genome = [q for q in qry_genome if q[0] != ref_genome[0]]
    constant_regions, SNP_matrix = ref_align(ref_genome, qry_genome)
    write_down(constant_regions, SNP_matrix, os.path.join(out_folder, 'snp_matrix'))
    # infer phylogeny
    # assign ancestral
    
    
    
    # get specific reads
    # infer SNPs
    # align reads onto reference
    # get SNP coverage
    # place SNPs onto phylogeny