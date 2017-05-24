import sys, os, utils, numpy as np, pandas as pd

def run_lastal(data) :
    r, q = data
    # run lastal
    # infer difference
    return align

def ref_align(ref, qry) :
    aligns = map(run_lastal, [[ref, q] for q in [ref] + qry])
    # get core genome
    # generate matrix
    return matrix

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
    
    if 'sparse_ref' in params :
        if 'external_ref' in params :
            raise Exception('cannot have both sparse_ref and external_ref')
        ref_genome = utils.get_genome_file(params['sparse_ref'], data)[0]
    elif 'external_ref' in params :
        ref_genome = params['external_ref']
    else :
        raise Exception('need to have either sparse_ref or external_ref')
    
    # align other genomes onto reference
    qry_genome = []
    if 'sparse_qry' in params :
        sparse_qrys = params['sparse_qry'].split(',')
        qry_genome.extend(utils.get_genome_file(params['sparse_qry'], data))
    if 'external_ref' in params :
        with open(params['external_qry']) as fin :
            for line in fin :
                qry_genome.append(line.strip().split()[:2])
    
    # get specific reads
    reads = get_specific_reads(**params)
    # infer SNPs
    # infer phylogeny
    # align reads onto reference
    # get SNP coverage
    # place SNPs onto phylogeny