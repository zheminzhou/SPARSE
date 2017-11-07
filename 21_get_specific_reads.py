import sys, os, gzip, glob

if __name__ == '__main__' :
    params = dict([p.split('=', 1) for p in sys.argv[1:]])
    params['ratio'] = float(params.get('ratio', 0.5))
    
    ref_id, workspace = params['ref_id'], params['workspace']
    ref_ids = ref_id.split(',')
    
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
    print 'Reference specific reads are saved in {0}'.format(output)
