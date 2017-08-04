import sys, os

cutoff = 0.01
mode = 'proportion'
polish = 'raw'

if __name__ == '__main__' :
    inputs = []
    for param in sys.argv[1:] :
        if param.startswith('cutoff') :
            cutoff = float(param.split('=')[1])
        elif param.startswith('mode') :
            mode = param.split('=')[1]
        elif param.startswith('polish') :
            polish = param.split('=')[1]
        else :
            inputs.append(param)
    
    outputs = {}
    min_props = {}
    for id, fname in enumerate(inputs) :
        if os.path.isdir(fname) :
            inputs[id] = os.path.join(fname, 'profile.txt')
            
        outputs[inputs[id]] = {}
        with open(inputs[id]) as fin :
            for line in fin :
                if line.startswith('Total') :
                    tot_reads = int(line.strip().split()[1])
                    effective_reads = float(line.strip().split()[2])
                    min_prop = max(cutoff*10000/tot_reads, cutoff)
                    if mode == 'number' :
                        min_prop = min_prop * tot_reads
                    min_props[inputs[id]] = min_prop
                elif line.startswith('U') :
                    continue
                else :
                    part = line.strip().split('\t')
                    taxonomy = part[3].rsplit(' ', 1)[0]
                    if mode != 'number' :
                        if float(part[1]) >= min_prop/2.0 :
                            outputs[inputs[id]][taxonomy] = outputs[inputs[id]].get(taxonomy, 0.0) + float(part[1])
                    else :
                        if float(part[2])*effective_reads >= min_prop/2.0 :
                            outputs[inputs[id]][taxonomy] = outputs[inputs[id]].get(taxonomy, 0.0) + float(part[2])*effective_reads
    
    taxonomies = {}
    for fname, hits in outputs.iteritems() :
        for taxonomy, prop in hits.iteritems() :
            if prop >= min_props[fname] and prop >= -taxonomies.get(taxonomy, [0.0])[0] :
                taxonomies[taxonomy] = [-prop]
    tax = []
    for t, p in sorted(taxonomies.iteritems(), key=lambda x:(x[1], x[0])) :
        tt = t.rsplit('|', 1)[0]
        if tt in taxonomies :
            taxonomies[t] = taxonomies[tt] + [p]
        tax.append([t, taxonomies[t]])
    taxonomies = sorted(tax, key=lambda x:(x[1], x[0]))
    print '#Level\tTaxonomy\t{0}'.format('\t'.join([f.rsplit('/', 1)[0] for f in inputs]))
    for tax in taxonomies :
        t = tax[0].replace('nan', '')
        if polish in ('very-clean', 'very-very-clean') :
            if t.find(' sp. ') >0 :
                continue
        if not t.endswith('|') :
            level = len(t.split('|'))
            level = {7:'Genus', 8:'Species', 9:'Subspecies', 10:'Reference'}[level]
            if level != 'Species' and polish != 'raw' :
                continue
            if polish == 'very-very-clean' :
                t = t.rsplit('|', 1)[-1]
            if mode != 'number' :
                print '{0}\t{1}\t{2}'.format(level, t, '\t'.join(['{:.3f}'.format(outputs[fname].get(tax[0], 0.0)) for fname in inputs]))
            else :
                print '{0}\t{1}\t{2}'.format(level, t, '\t'.join(['{0}'.format(int(round(outputs[fname].get(tax[0], 0.0)/100, 0))) for fname in inputs]))
