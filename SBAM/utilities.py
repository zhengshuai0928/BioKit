from SBAM.fasta import *
from SBAM.table import *

def extractGff3(gff3_path, fasta_path, feature, out_path):
    out = open(out_path, 'w')
    fa = Fasta(fasta_path)
    gf = Gff3(gff3_path)
    for anno in gf:
        seqid  = '>' + anno[0]['seqid']
        strand = anno[0]['strand']
        out.write(seqid + '_' + '_'.join(anno[0]['start':'score']) + '\n')
        seq = ''
        for line in anno:
            if line['type'] == feature:
                start = int(line['start']) - 1
                end = int(line['end'])
                seq += fa[seqid][start:end]
        if strand == '-':
            seq = seq[::-1]
        out.write(seq + '\n')