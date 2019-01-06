from SBAM.fasta import Fasta
import os
import glob

class Cluster:
    
    def __init__(self, clstr_path):
        self._clstr = {}
        self.fasta = ''
        clstr_file = open(clstr_path)
        last_id = clstr_file.readline().strip('\n')
        seq_ids = []
        for line in clstr_file:
            line = line.strip('\n')
            if line.startswith('>'):
                self._clstr[last_id] = seq_ids
                seq_ids = []
                last_id = line
                continue
            id_start = line.index('>')
            id_end = line.index('.')
            seq_ids.append(line[id_start:id_end])
        self._clstr[last_id] = seq_ids
        clstr_file.close()

    def __getitem__(self, clstr_no):
        pre = '>Cluster '
        return self._clstr[pre + str(clstr_no)]

    def getCluSeq(self, clstr_no, fasta_path, out_path='', mode='w'):
        if not self.fasta:
            self.fasta = Fasta(fasta_path)
        out = open(out_path, mode)
        for seq_id in self[clstr_no]:
            out.write(seq_id + '\n')
            out.write(self.fasta[seq_id] + '\n')
        out.close()

    def muscle(self, clstr_no, fasta_path, out_path='', mode='w', fmt=''):
        tmp_file = '.cluster.seq.tmp'
        self.getCluSeq(clstr_no, fasta_path, tmp_file)
        muscle_str = ' '.join(('muscle', '-in', tmp_file, '-out',
                                out_path, fmt))
        os.system(muscle_str)
        os.remove(tmp_file)

    def muscleShow(self, clstr_no, fasta_path, fmt='-clw'):
        out_path = 'cluster_%s_muscle.aln'% clstr_no
        self.muscle(clstr_no, fasta_path, out_path, fmt=fmt)
        out = open(out_path)
        for line in out:
            if not line: break
            print(line, end='')
        out.close() 

def _blast2seqs():
    seq1 = input('Please input first  seq: ')
    seq2 = input('Please input second seq: ')
    _type = ''
    while _type not in ['nucl', 'prot']:
        _type = input('Type: (nucl/prot): ')
    qpath = '.seq1.fas'
    dbpath = '.seq2.fas'
    fas_file = open(qpath, 'w')
    fas_file.write('>seq1' + '\n')
    fas_file.write(seq1.upper() + '\n')
    fas_file.close()   
    fas_file = open(dbpath, 'w')
    fas_file.write('>seq2' + '\n')
    fas_file.write(seq2.upper() + '\n')
    fas_file.close()
    db_cmd = 'makeblastdb -in %s -dbtype %s \
              -logfile .makeblastdb_log'%(dbpath, _type)
    os.system(db_cmd)
    if _type == 'nucl':
        task = 'blastn'
    if _type == 'prot':
        task = 'blastp'
    blast_cmd = '%s -query %s -db %s -task %s'% \
                (task, qpath, dbpath, task)
    out = os.popen(blast_cmd).read()
    os.wait()
    print('-' * 50)
    out_list = out.split('\n')
    for line in out_list:
        if line == '':
            continue
        print(line)
    print('-' * 50)
    os.remove('.makeblastdb_log.perf')
    os.remove('.makeblastdb_log')
    os.remove(qpath)
    os.remove(dbpath)
    for file_name in glob.glob(dbpath + '.*'):
        os.remove(file_name)
    go_on = ''
    while not go_on in ['Y', 'y', 'N', 'n']:
        go_on = input('Go on? (y/Y/n/N): ')
    return go_on

def blast2seqs():
    go_on = 'y'
    while go_on in ['Y', 'y']:
        go_on = _blast2seqs()