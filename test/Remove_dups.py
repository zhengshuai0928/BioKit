# Simplified Bioinformatics Analysis Module
import os
import sys
import re
import time
import gc

def progressBar(symbol, symbol_num, count, total):
    symbol_count = int(count / (total / symbol_num))
    sys.stdout.write('\r')
    sys.stdout.write(('[%-50s]%.2f%%' % ((symbol 
                      * symbol_count), count / total * 100)))
    sys.stdout.flush()

class Fasta:

    def __init__(self, fasta_path):
        print('Reading in...')        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        count = 0
        _seq = ''
        fasta_file = open(fasta_path)
        _id_last = fasta_file.readline().strip()
        for line in fasta_file:
            count += 1
            progressBar('#', 50, count, line_num)
            line = line.strip()
            if line.startswith('>'):
                self._fasta[_id_last] = _seq
                _seq = ''
                _id_last = line
                continue
            _seq += line.upper()
        self._fasta[_id_last] = _seq
        fasta_file.close()
        sys.stdout.write('\n')
        self.setAttrs()

    def setAttrs(self):
        self.number = len(self._fasta)
        self.total_len = 0
        self._idlist = list(self._fasta.keys())
        for key in self._fasta:
            self.total_len += len(self._fasta[key])
    
    def __getitem__(self, _id_str):  
        try:
            return self._fasta[_id_str]
        except KeyError:
            print('Fuzzy matching...')
            for full_id in self._idlist:
                if re.search(_id_str, full_id):
                    print('Matched:', full_id)
                    return self._fasta[full_id]
    
    def rmDups(self):   
        reverse_dict = {}
        visited_seqs = []
        dup_keys = []
        for _id in self._idlist:
            if self._fasta[_id] in visited_seqs:
                dup_keys.append(_id)
                print(_id) 
                print('  duplicates', 
                       reverse_dict[self._fasta[_id]])
                continue
            visited_seqs.append(self._fasta[_id])
            reverse_dict[self._fasta[_id]] = _id
        if dup_keys:
            for key in dup_keys:
                self._fasta.pop(key)
            print('-' * 20)
            print('%s duplicated seqs are removed.'% len(dup_keys))
        else:
            print('No duplicates!')

        del visited_seqs
        del reverse_dict
        gc.collect()
        self.setAttrs()

if __name__ == '__main__':
    fa = Fasta(sys.argv[1])
    time_begin = time.time()
    fa.rmDups()
    print(time.time() - time_begin)
