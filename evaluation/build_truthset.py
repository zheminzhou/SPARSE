import sys, os, numpy as np, gzip, subprocess, re, time

tax_col = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def get_taxonomy(tax_db, tax_col = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) :
    names = {}
    children = {}
    nodes = {'1':{}}
    category = { t:1 for t in tax_col }
    with open(os.path.join(tax_db, 'nodes.dmp')) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            nodes[part[0]] = {part[4]: part[0]} if part[4] in category else {}
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

def tax_assign(taxid, tax1, tax_col) :
    return [ str(tax1.get(taxid, {}).get(tcl, 0)) for tcl in tax_col ]

if __name__ == '__main__' :
    tax1 = get_taxonomy(sys.argv[3])
    
    print '# {0}\t{1}'.format('Read', '\t'.join(tax_col))
    data_type = sys.argv[1]
    if data_type == 'CAMI' :
        with gzip.open(sys.argv[2]) as fin :
            for id, line in enumerate(fin) :
                if id > 3 :
                    part = line.strip().split('\t')
                    taxid = part[2]
                    res = tax_assign(taxid, tax1, tax_col)
                    print '{0}\t{1}'.format(part[0], '\t'.join(res))
    elif data_type == 'KEGG' :
        kegg_save = {'dda': '69223'}
        with gzip.open(sys.argv[2]) as fin :
            for id, line in enumerate(fin) :
                if id > 0 :
                    part = line.strip().split('\t')
                    kegg_id = part[1][:3]
                    if kegg_id not in kegg_save :
                        url = 'curl http://www.genome.jp/kegg-bin/show_organism?org=' + kegg_id
                        html = subprocess.Popen(url.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
                        taxid = re.findall(r'wwwtax.cgi\?mode=Info\&id=(\d+)', html)[0]
                        kegg_save[kegg_id] = taxid
                    else :
                        taxid = kegg_save[kegg_id]
                    res = tax_assign(taxid, tax1, tax_col)
                    print '{0}\t{1}'.format(part[1], '\t'.join(res))
    elif data_type == 'NCBI' :
        ncbi_save = {
        'NC_008022' :'370552',
        'NC_011566' :'225849',
        'NC_012225' :'565034',
        'NC_012751' :'572265',
        'NC_013971' :'716540',
        'NC_002180' :'1986029'}
        with gzip.open(sys.argv[2]) as fin :
            fin.readline()
            for line in fin :
                    part = line.strip().split('\t')
                    gid = part[1].split('.', 1)[0]
                    if gid not in ncbi_save :
                        try:
                            url = 'curl https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&id=' + gid
                            html = subprocess.Popen(url.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
                            taxid = re.findall(r'<Link>\s+<Id>(\d+)</Id>', html)[0]
                        except :
                            time.sleep(0.3)
                            url = 'curl https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + gid
                            print >> sys.stderr, url
                            html = subprocess.Popen(url.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
                            taxid = re.findall(r'"taxon" ,\s+tag\s+id (\d+) ', html)[0]
                        ncbi_save[gid] = taxid
                    else :
                        taxid = ncbi_save[gid]
                    res = tax_assign(taxid, tax1, tax_col)
                    print '{0}\t{1}'.format(part[1], '\t'.join(res))
