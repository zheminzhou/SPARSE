import argparse, numpy as np, sys

def _create(argv) :
    from modules.A1_db_create import db_create
    parser = argparse.ArgumentParser(description='''Create an empty database. Use "SPARSE.py index" to fill in references later''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the new database. ', required=True)
    args = parser.parse_args(argv)
    return db_create(args.__dict__)

def _index(argv) :
    from modules.A2_db_index import db_index
    parser = argparse.ArgumentParser(description='''Fill reference genomes into a SPARSE database''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the database. ', required=True)
    parser.add_argument('-u', '--update', help='Index the current RefSeq database. Overwrite SEQLIST. ', action='store_true')
    parser.add_argument('-s', '--seqlist', help='Load in a tab-delimited file. ')
    parser.add_argument('-t', '--n_thread', help='Number of threads to use. Default: 20 ', type=int, default=20)
    args = parser.parse_args(argv)
    return db_index(args.__dict__)

def _MapDB(argv) :
    from modules.A3_db_MapDB import db_MapDB
    parser = argparse.ArgumentParser(description='''Generate sub-databases for read alignments.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the database. ', required=True)
    parser.add_argument('-m', '--MapDB', help='Name for sub-database. ', required=True)
    parser.add_argument('-s', '--seqlist', help='A tab-delimited list of reference to include. Required. ', required=True)
    parser.add_argument('-t', '--n_thread', help='Number of threads to use. Default: 20 ', type=int, default=20)
    parser.add_argument('--malt', dest='dbtype', help='Use MALT instead of bowtie2 [default]', default='bowtie2', action='store_const', const='malt')
    parser.add_argument('--append', dest='mode', help='Append to existing database instead of overwrite [default]', default='overwrite', action='store_const', const='append')
    args = parser.parse_args(argv)
    return db_MapDB(args.__dict__)

def _predict(argv) :
    from modules.B2_query_reads import query_read
    parser = argparse.ArgumentParser(description='''Alignment based taxonomy prediction.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the database. ', required=True)
    parser.add_argument('-m', '--MapDB', help='Comma delimited names for sub-database. Default: representative,subpopulation,Virus', default='representative,subpopulation,Virus')
    parser.add_argument('-w', '--workspace', help='Folder name for all outputs and intermediate results.', required=True)
    parser.add_argument('-1', '--r1', help='SE read or first part of PE reads.', required=True)
    parser.add_argument('-2', '--r2', help='Second part of PE reads.')
    parser.add_argument('-t', '--n_thread', help='Number of threads to use. Default: 20 ', type=int, default=20)
    args = parser.parse_args(argv)
    return query_read(args.__dict__)

def _sample(argv) :
    from modules.B3_query_sample import query_sample
    parser = argparse.ArgumentParser(description='''Rapid query of assembly or the whole read file.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the database. ', required=True)
    parser.add_argument('-q', '--query', help='A genome in fasta format, or a set of reads in fastq format. Required. ', required=True)
    parser.add_argument('--read', dest='dtype', help='query is read in fastq format rather than assembly', default='fasta', action='store_const', const='read')
    args = parser.parse_args(argv)
    return query_sample(args.__dict__)

def _query(argv) :
    from modules.B1_query_metadata import query_metadata
    parser = argparse.ArgumentParser(description='''Retrieve metadata for a set of references.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the database. ', required=True)
    parser.add_argument('-s', '--seqlist', help='output result. Default: to screen', default=None)
    parser.add_argument('-d', '--default', help='Default MapDB criteria for updates. Choose from:\nrepresentative, subpopulation, Virus, Eukaryota', default=None)
    parser.add_argument('--min', help='Minimum size of genomes to show', type=int, default=None)
    parser.add_argument('--max', help='Maximum size of genomes to show', type=int, default=None)
    parser.add_argument('--group', help='Filter using the prefix of barcode addresses', default=None)
    parser.add_argument('--tag', help='''Filter by relationships between different level of barcodes. i.e., 
"p!=r;p==a" gets references that have the same numbers in p groups and a groups, but different between p groups and r groups''', default=None)
    parser.add_argument('--index', help='Filter by index.', default=None)
    parser.add_argument('--barcode', help='Filter by barcode.', default=None)
    parser.add_argument('--assembly_accession', help='Filter by assembly_accession.', default=None)
    parser.add_argument('--refseq_category', help='Filter by refseq_category.', default=None)
    parser.add_argument('--assembly_level', help='Filter by assembly_level.', default=None)
    parser.add_argument('--taxid', help='Filter by taxid.', default=None)
    parser.add_argument('--organism_name', help='Filter by organism_name.', default=None)
    parser.add_argument('--species', help='Filter by species.', default=None)
    parser.add_argument('--genus', help='Filter by genus.', default=None)
    parser.add_argument('--family', help='Filter by family.', default=None)
    parser.add_argument('--order', help='Filter by order.', default=None)
    parser.add_argument('--class', help='Filter by class.', default=None)
    parser.add_argument('--phylum', help='Filter by phylum.', default=None)
    parser.add_argument('--kingdom', help='Filter by kingdom.', default=None)
    parser.add_argument('--superkingdom', help='Filter by superkingdom.', default=None)
    
    args = parser.parse_args(argv)
    return query_metadata(args.__dict__)

def _update(argv) :
    from modules.C1_modify_metadata import update
    parser = argparse.ArgumentParser(description='''Update metadata.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-n', '--dbname', help='Name for the database. ', required=True)
    parser.add_argument('-s', '--seqlist', help='A tab-delimited list of reference. Required. ', required=True)
    args = parser.parse_args(argv)
    return update(args.__dict__)

def _SSR(argv) :
    from modules.C2_get_specific_reads import get_SSR
    parser = argparse.ArgumentParser(description='''Extract species specific reads associated with particular references.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-w', '--workspace', help='Folder name for all outputs and intermediate results.', required=True)
    parser.add_argument('-i', '--ref_id', help='Comma delimited reference indexes to extract.', required=True)
    parser.add_argument('-r', '--ratio', help='The minimum confidential ratio to report.', default=0.5, type=float)
    args = parser.parse_args(argv)
    return get_SSR(args.__dict__)

def _report(argv) :
    from modules.D1_sparse_parse import report
    parser = argparse.ArgumentParser(description='''Generate a flat-table report for multiple runs.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--level', help='Level to report, default: s.', default='s')
    parser.add_argument('workspaces', metavar='WS', nargs='+', help='Folder names for all outputs.')
    args = parser.parse_args(argv)
    report(args.level, args.workspaces)




def SPARSE() :
    try :
        parser = argparse.ArgumentParser(description='''SPARSE (Strain Prediction and Analysis with Representative SEquences)''', formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('command', help='task to run in SPARSE.')
        sys.argv[0] = ' '.join(sys.argv[:2])
        exec '_{0}(sys.argv[2:])'.format(sys.argv[1])
    except Exception as e :
        print '''
Program: SPARSE (Strain Prediction and Analysis with Representative SEquences)

Usage:   SPARSE.py <command> [options]

Commands:
  create        Create empty folder structures for a new SPARSE database
  index         Load in a list of assemblies (in RefSeq format) and index them into a SPARSE database
  query         Query metadata info in a SPARSE database
  update        Update metadata info in a SPARSE database
  MapDB         Create bowtie2 or MALT sub-databases for metagenomic reads
  predict       Align reads onto MapDB and do taxonomic predictions
  sample        Compare an assembly with all genomes in a SPARSE database
  SSR           Extract species-specific reads (SSR) from a SPARSE read-mapping result
  report        Reformat and merge multiple SPARSE outputs into a flat table.

Use SPARSE.py <command> -h to get help for each command.
'''
        import traceback
        traceback.print_exception(*sys.exc_info())
        
        sys.exit(0)


if __name__ == '__main__' :
    SPARSE()
