# msh_barcoding.py

import sys, os, requests
import shutil, utils
from multiprocessing import Pool
from time import sleep

def worker(lnk, n_thread=16) :
    context = zmq.Context()
    sock = context.socket(zmq.REQ)
    sock.connect(lnk)
    
    result = {"n_thread": n_thread}
    while True:
        sock.send_json(status)
        work = sock.recv_json()
        print json.dumps(work)
        if work['command'] == 'exit' :
            break
        elif work['command'] == 'sleep' :
            sleep(work['argv'])
        else :
            result = Tasks.__dict__[work['command']](work['argv'])
        sock.send_json(result)
        sock.recv_json()
        
class Tasks(object) :
    @staticmethod
    def mash_proc(data) :
        idx, file_link, url_link, params = data
        if not os.path.isfile(file_link) :
            file_link = 'nan'
        if file_link == 'nan' :
            fname = url_link.rsplit('/',1)[-1]
            subprocess.Popen('wget -O {1} {0}'.format(url_link, fname).split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
        else :
            fname = file_link
        try :
            sha = utils.fileSHA(fname)
            if os.path.isfile(fname.rsplit('/',1)[-1] + '.msh') :
                os.unlink(fname.rsplit('/',1)[-1] + '.msh')
            msh_file = utils.get_mash(fname, fname.rsplit('/',1)[-1], is_read=False, **params)
            if file_link == 'nan' :
                os.unlink(fname)
            return idx, sha, msh_file, fname
        except :
            return idx, '', '', ''

    def create_db(data) :
        bt_build, prefix, id, genomes = data
        dbname = '{0}.{1}'.format(prefix, id)
        seq_tax = []
        for g in genomes :
            if ':' in g[1] :
                if g[1].upper().endswith('.GZ') :
                    subprocess.Popen('wget -O {0}.tmp.gz {1}'.format(dbname, g[1]).split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE).communicate()
                    fname = '{0}.tmp.gz'.format(dbname)
                else :
                    subprocess.Popen('wget -O {0}.tmp {1}'.format(dbname, g[1]).split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE).communicate()
                    fname = '{0}.tmp'.format(dbname)
            else :
                fname = g[1]
            show_cmd = 'zcat' if fname.upper().endswith('.GZ') else 'cat'
            show_run = subprocess.Popen("{0} {1}|tee -a '{2}'|grep '^>'".format(show_cmd, fname, dbname), shell=True, stdout=subprocess.PIPE)
            for line in iter(show_run.stdout.readline, r'') :
                seqname = line[1:].strip().split()[0]
                seq_tax.append([seqname, g[0]])
            if g[1] != fname :
                os.unlink(fname)
        with gzip.open('{0}.taxa.gz'.format(dbname), 'w') as fout :
            for s, b in seq_tax :
                fout.write('{0}\t{1}\n'.format(s, b))
        r = subprocess.Popen('{0} -o 3 {1} {1}'.format(bt_build, dbname).split(), stdout=subprocess.PIPE)
        r.communicate()
        if r.returncode == 0 :
            subprocess.Popen('gzip {0}'.format(dbname).split()).communicate()
        else :
            for fname in glob.glob(dbname + '.*') :
                os.unlink(fname)
        return [prefix, id, r.returncode]
