import os, sys, pandas as pd, numpy as np, json, msgpack
import utils



if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    assert 'query' in params and os.path.isfile(params['query']), 'no query'
    
    existing_data = os.path.join(params['dbname'], 'db_metadata.msg')
    assert existing_data, 'no data in the database.'
    data = pd.read_msgpack(existing_data)
    
    if params.get('dtype', 'fasta') == 'read' :
        msh_file = utils.get_mash(params['query'], is_read=True, **params)
        result = utils.run_mash([msh_file, None, params['n_thread'], params], is_read=True)
    else :
        msh_file = utils.get_mash(params['query'], is_read=False, **params)
        result = utils.run_mash([msh_file, None, params['n_thread'], params])
    os.unlink(msh_file)
    
    if len(result) > 0 :
        r_id = np.array([r[2] for r in result])
        result = np.array([r[1].split('.') for r in result])
        result, r_id = result[(1-r_id >= 0.98* (1-r_id[0])) & (r_id <= params['barcode_dist'][0])], r_id[(1-r_id >= 0.98* (1-r_id[0])) & (r_id <= params['barcode_dist'][0])]

        groups = {}
        m = {'a'+k:[n, a] for k, n, a in data[['index', 'organism_name', 'assembly_accession']].as_matrix()}
        matches = [ dict(record='.'.join(r), similarity=1-i, organism_name=m[r[-1]][0], assembly_accession=m[r[-1]][1]) for r, i in zip(result, r_id) ]
        for id, (dcut, dgroup) in enumerate(zip(params['barcode_dist'], result.T[:-1])) :
            dgroup[r_id > dcut] = ''
            g = np.unique(dgroup, return_index=True)
            tags = ['.'.join(r) for r in result[g[1], :(id+1)] if r[-1] != '']
            info = [ [i, -id, '.'.join(hit)] for i, hit in zip(r_id[g[1]], result[g[1]]) if hit[id] != '' ]
            for t, i in zip(tags, info) :
                groups[t] = i

        #for qry, hit, dist in result :
            #if (1-dist) < (1 - result[0][2])*0.95 :
                #break
            #code = []
            #for id, g in enumerate(hit.split('.')[:-1]) :
                #if params['barcode_dist'][id] >= dist :
                    #code.append(g[1:])
                    #if tuple(code) not in groups :
                        #groups[tuple(code)] = [dist, -id, hit]
                #else :
                    #break
        #matches.append(dict(record=hit, similarity=1.0-dist, organism_name=data[data['index'] == hit.rsplit('.', 1)[-1][1:]]['organism_name'].tolist()[0]))
        groups = [dict(group=c, similarity=1.0-d[0]) for c, d in sorted(groups.iteritems(), key=lambda x:x[1])]
        for g in groups :
            g.update(utils.retrieve_info(g['group'], data=data, **params))
        #result = groups[0]['group']
    else :
        groups, matches, result = [], [], 'unknown'
    print json.dumps(dict(groups=groups, matches=matches), sort_keys=True, indent=2)