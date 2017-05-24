import sys, os

cutoff = 0.005

if __name__ == '__main__' :
    inputs = sys.argv[1:]
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
                    min_prop = max(cutoff*10000/tot_reads, cutoff)
                    min_props[inputs[id]] = min_prop
                elif line.startswith('U') :
                    continue
                else :
                    part = line.strip().split('\t')
                    taxonomy = part[3].rsplit(' ', 1)[0]
                    if float(part[1]) >= min_prop/5.0 :
                        outputs[inputs[id]][taxonomy] = outputs[inputs[id]].get(taxonomy, 0.0) + float(part[1])
    
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
    print '#Level\tTaxonomy\t{0}\n'.format('\t'.join(inputs))
    for tax in taxonomies :
        t = tax[0].replace('nan', '')
        if not t.endswith('|') :
            level = len(t.split('|'))
            level = {7:'Genus', 8:'Species', 9:'Subspecies', 10:'Reference'}[level]
            print '{0}\t{1}\t{2}'.format(level, t, '\t'.join(['{:.3f}'.format(outputs[fname].get(tax[0], 0.0)) for fname in inputs]))
