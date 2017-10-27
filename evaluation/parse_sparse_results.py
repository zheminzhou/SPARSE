import sys, os, numpy as np, gzip, subprocess, re, pandas as pd

tax_col = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def get_taxonomy(tax_db, tax_col = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) :
    names = {}
    children = {}
    nodes = {'1':{}}
    category = { t:1 for t in tax_col }
    named = {}
    with open(os.path.join(tax_db, 'names.dmp')) as fin :
        for line in fin :
            part = line.strip().split('\t')
            if part[6] in ('type material', 'authority') :
                named[part[0]] = 1
                
    with open(os.path.join(tax_db, 'nodes.dmp')) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            nodes[part[0]] = {part[4]: int(part[0])*named.get(part[0], -1)} if part[4] in category else {}
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
                nodes[c].update(nodes[qq])
        q = nq
    return nodes

def sparse_interpret(refset, tax1, allow_fuzzy=0, tax_col = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) :
    g = {}
    for barcode, taxid in refset[['barcode', 'taxid']].as_matrix() :
        barcodes = barcode.split('.')
        pg, eg = barcodes[3], barcodes[6]
        if pg not in g :
            g[pg] = {}
        
        t1 = tax1.get(taxid, {}).get('species', 0)
        if t1 not in g[pg] :
            g[pg][t1] = {eg:1}
        else :
            g[pg][t1][eg] = 1
            
    for pg in g.keys() :
        g[pg] = np.array([[tid,len(egs)] for tid, egs in g[pg].iteritems()])
        g[pg] = g[pg][(g[pg].T[1] >= np.max(g[pg].T[1])*0.1) | (g[pg].T[1] >= 4)]
        
        if np.sum(g[pg].T[0] > 0) > 0 :
            g[pg] = g[pg][g[pg].T[0] > 0]
        elif np.sum(g[pg].T[0] < 0) > 0 :
            if allow_fuzzy :
                g[pg] = -g[pg][g[pg].T[0] < 0]
            else :
                g[pg] = g[pg][g[pg].T[0] < 0]
        s = float(np.sum(g[pg].T[1]))
        g[pg] = { x:y/s for x, y in g[pg] }

    return {idx:g[barcode.split('.')[3]] for idx, barcode in refset[['index', 'barcode']].as_matrix()}

if __name__ == '__main__' :
    refset = pd.read_msgpack(os.path.join(sys.argv[1], 'db_metadata.msg'))
    tax1 = get_taxonomy(sys.argv[3])
    if len(sys.argv) > 4 :
        allow_fuzzy=int(sys.argv[4])
    else :
        allow_fuzzy=0
    refconv = sparse_interpret(refset, tax1, allow_fuzzy)
    
    print '# {0}\t{1}'.format('Read', '\t'.join(tax_col))
    with gzip.open(sys.argv[2]) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            if len(part) <= 2 or part[4].startswith('-1') :
                result = ['0'] * len(tax_col)
            else :
                res = [{} for tcl in tax_col]
                for p in part[4:] :
                    ref, prop, d1, d2 = p.split(':')
                    prop, d1, d2 = float(prop), float(d1), float(d2)
                    div = max(d1, d2)
                    taxids = refconv[ref]
                    
                    for pid, (taxid, taxprop) in enumerate(sorted(taxids.items(), key=lambda x:(x[1], -x[0]))) :
                        if taxid < 0 :
                            taxid, tozero = -taxid, 1
                        else :
                            tozero = 0
                        taxid = str(taxid)
                        for id, tcl in enumerate(tax_col) :
                            ta1 = abs(tax1.get(taxid, {}).get(tcl, 0))
                            if tcl == 'species' and (div >= 5 or tozero) :
                                ta1 = 0
                            elif tcl == 'genus' and div >= 10 :
                                ta1 = 0
                            
                            res[id][ta1] = res[id].get(ta1, 0.) + prop+(pid/100000.0)
                            
                rr2 = [[] for tcl in tax_col]
                for id, r in enumerate(res) :
                    m = [m[0] for m in sorted(r.items(), key=lambda x:-x[1]) if m[1] >= 0.8]
                    rr2[id].extend(m if len(m) > 0 else [0])
                result = [','.join(np.unique(r).astype(str).tolist()) for r in rr2]
            print '{0}\t{1}'.format(part[1], '\t'.join(result))
