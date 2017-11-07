import sys, os, re

cutoff = 0.001
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
    species = {}
    for id, fname in enumerate(inputs) :
        if os.path.isdir(fname) :
            inputs[id] = os.path.join(fname, 'profile.txt')
            
        outputs[inputs[id]] = {}
        with open(inputs[id]) as fin :
            for line in fin :
                if line.startswith('Total') :
                    tot_reads = int(line.strip().split()[1])
                    effective_reads = float(line.strip().split()[2])
                    #min_freq = cutoff  if mode == 'number' else cutoff * tot_reads
                    min_props[inputs[id]] = cutoff
                elif line.startswith('U') :
                    continue
                else :
                    part = line.strip().split('\t')
                    cnt = float(part[1]) if mode != 'number' else float(part[2])*effective_reads
                    if cnt < cutoff : continue
                    taxonomy = part[3].rsplit(' ', 1)[0]
                    labels = taxonomy.split('|')
                    if len(labels) > 7 :
                        try:
                            main, other = labels[7].split(' (')
                            other = other[:-1].split('/')
                            if main not in species :
                                species[main] = {}
                            species[main].update({o:1 for o in other })
                        except :
                            if labels[7] not in species :
                                species[labels[7]] = {}
                    taxonomy = re.sub(r' \([^\)]+\)', r'', taxonomy)
                    outputs[inputs[id]][taxonomy] = outputs[inputs[id]].get(taxonomy, 0.0) + cnt
    
    for main, others in species.items() :
        if polish.find('very-clean') >= 0 :
            if main.startswith('*') :
                species[main] = ''
                continue
            others = [o for o in sorted(others.keys()) if not o.startswith('*')]
        else :
            others = sorted(others.keys())
        if len(others) > 0 :
            species[main] = '{0} ({1})'.format(main, '/'.join(others))
        else :
            species[main] = main
    
    taxonomies = {}
    for fname, hits in outputs.iteritems() :
        for taxonomy, prop in hits.iteritems() :
            if prop >= min_props[fname] and -prop < taxonomies.get(taxonomy, [0.])[0] :
                    taxonomies[taxonomy] = [-prop]
    tax = []
    for t, p in sorted(taxonomies.iteritems(), key=lambda x:(x[1], x[0])) :
        tt = t.rsplit('|', 1)[0]
        if tt in taxonomies :
            taxonomies[t] = taxonomies[tt] + p
        tax.append([t, taxonomies[t]])
    taxonomies = sorted(tax, key=lambda x:(x[1], x[0]))
    print '#Level\tTaxonomy\t{0}'.format('\t'.join([f.rsplit('/', 1)[0] for f in inputs]))
    for tax in taxonomies :
        t = tax[0].replace('-', '')
        if not t.endswith('|') :
            labels = t.split('|')
            level = len(labels)
            if level > 7 :
                if species.get(labels[7], '') == '' :
                    continue
                else :
                    labels[7] = species[labels[7]]
            
            level = {6:'Class', 7:'Genus', 8:'Species', 9:'Subspecies', 10:'Reference'}[level]
            if level != 'Species' and polish != 'raw' :
                continue
            if polish == 'very-very-clean' :
                t = labels[-1]
            else :
                t = '|'.join(labels)
            if mode != 'number' :
                print '{0}\t{1}\t{2}'.format(level, t, '\t'.join(['{:.3f}'.format(outputs[fname].get(tax[0], 0.0)) for fname in inputs]))
            else :
                print '{0}\t{1}\t{2}'.format(level, t, '\t'.join(['{0}'.format(int(round(outputs[fname].get(tax[0], 0.0)/100, 0))) for fname in inputs]))
