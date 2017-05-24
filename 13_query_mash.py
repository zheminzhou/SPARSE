import os, sys, pandas as pd, numpy as np, ujson, json, msgpack
import utils



if __name__ == '__main__' :
    params = utils.load_params(sys.argv)
    assert 'query' in params and os.path.isfile(params['query']), 'no query'
    
    existing_data = os.path.join(params['dbname'], 'db_metadata.msg')
    assert existing_data, 'no data in the database.'
    data = pd.read_msgpack(existing_data)
    
    if params.get('dtype', 'fasta') == 'read' :
        msh_file = utils.get_mash(params['query'], is_read=True, **params)
    else :
        msh_file = utils.get_mash(params['query'], is_read=False, **params)
    result = utils.run_mash([msh_file, None, params['n_thread'], params])
    os.unlink(msh_file)
    
    if len(result) > 0 :
        groups = {}
        matches = []
        for qry, hit, dist in result :
            if (1-dist) < (1 - result[0][2])*0.9 :
                break
            code = []
            for id, g in enumerate(hit.split('.')[:-1]) :
                if params['barcode_dist'][id] * 1.2 >= dist :
                    code.append(g[1:])
                    if tuple(code) not in groups :
                        groups[tuple(code)] = [dist, -id, hit]
                else :
                    break
            matches.append(dict(record=hit, similarity=1.0-dist, organism_name=data[data['index'] == hit.rsplit('.', 1)[-1][1:]]['organism_name'].tolist()[0]))
        groups = [dict(group='.'.join(['{0}{1}'.format(t,k) for t,k in zip(params['barcode_tag'][:len(c)], c)]), similarity=1.0-d[0]) for c, d in sorted(groups.iteritems(), key=lambda x:x[1])]
        for g in groups :
            g.update(utils.retrieve_info(g['group'], data=data, **params))
        result = groups[0]['group']
    else :
        groups, matches, result = [], [], 'unknown record'
    print json.dumps(dict(groups=groups, matches=matches, result=result), sort_keys=True, indent=2)