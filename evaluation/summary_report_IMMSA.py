import os, sys


tax_col = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

if __name__ == '__main__' :
    truth, result = sys.argv[1:]
    
    res = {}
    with open(result) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            res[part[0]] = part[-1]
    
    groups = [0,0,0,0]
    with open(truth) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split('\t')
            if part[0] in res :
                tr, rs = part[-1], res[part[0]]
                tr, rs = tr.split(','), rs.split(',')
                
                category = 4
                for ttr in tr :
                    for rrs in rs :
                        if ttr == rrs :
                            if ttr == '0' :
                                category = min(category, 1)
                            else :
                                category = 0
                                break
                        elif ttr != '0' and rrs != '0' :
                            category = min(category, 2)
                        elif ttr != '0' :
                            category = min(category, 3)
                if category == 0 :
                    groups[0] += 1
                elif category == 1 :
                    groups[3] += 1
                elif category == 2 :
                    groups[1] += 1
                    groups[2] += 1
                elif category == 3 :
                    groups[2] += 1
                else :
                    groups[1] += 1
    
    precision = float(groups[0])/(groups[0] + groups[1])
    recall = float(groups[0])/(groups[0] + groups[2])
    f1 = 2/(1/precision + 1/recall)
    print '{0}\tPrecision:\t{1}\tRecall:\t{2}\tRaw:\t{3}\t{4}\t{5}\t{6}'.format(result, precision, recall, *groups)