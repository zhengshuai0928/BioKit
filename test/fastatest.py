import sys
import os
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

class Fasta2:

    def __init__(self, fasta_path):
        print('Reading in...')        
        self._fasta = {}
        fasta_lines = open(fasta_path).readlines()
        line_num = len(fasta_lines)
        count = 0
        _seq = ''
        _id_last = fasta_lines.pop(0)
        
        for line in fasta_lines:
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
        sys.stdout.write('\n')
        del fasta_lines
        gc.collect()

if __name__ == '__main__':
    time_begin = time.time()
    fa = Fasta2(sys.argv[1])
    print(time.time() -  time_begin)