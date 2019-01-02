import os
import sys

class Fasta:

    def __init__(self, fasta_path):        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        
        # Progress bar setting
        progress_symbol = '#'
        symbol_num = 50
        each_symbol = line_num / symbol_num
        #
        count = 0
        for line in open(fasta_path):
            # Progress bar printing
            count += 1
            symbol_count = int(count / each_symbol) 
            sys.stdout.write('\r')
            sys.stdout.write(('[%-50s]%.2f%%' % ((progress_symbol * symbol_count), 
                               count / line_num * 100)))
            sys.stdout.flush()
            #
            line = line.strip()
            if not line: 
                continue
            if line.startswith('>'):
                continue
        sys.stdout.write('\n')

class Fasta1:

    def __init__(self, fasta_path):        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        
        # Progress bar setting
        progress_symbol = '#'
        symbol_num = 50
        each_symbol = line_num / symbol_num
        #
        count = 0
        for line in open(fasta_path):
            # Progress bar printing
            count += 1
            symbol_count = int(count / each_symbol) 
            sys.stdout.write('\r')
            sys.stdout.write(('[%-50s]%.2f%%' % ((progress_symbol * symbol_count), 
                               count / line_num * 100)))
            sys.stdout.flush()
            #
            line = line.strip()
            if not line: 
                continue
            if line.startswith('>'):
                _id = line
                self._fasta[_id] = ''
                continue
            self._fasta[_id] += line.upper()
        sys.stdout.write('\n')

class Fasta2:

    def __init__(self, fasta_path):        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        
        # Progress bar setting
        progress_symbol = '#'
        symbol_num = 50
        each_symbol = line_num / symbol_num
        #
        count = 0
        _seq = ''
        for line in open(fasta_path):
            # Progress bar printing
            count += 1
            symbol_count = int(count / each_symbol) 
            sys.stdout.write('\r')
            sys.stdout.write(('[%-50s]%.2f%%' % ((progress_symbol * symbol_count), 
                               count / line_num * 100)))
            sys.stdout.flush()
            #
            line = line.strip()
            if not line: 
                continue
            if line.startswith('>'):
                _id = line
                if _seq:
                    self._fasta[last_id] = _seq
                    _seq = ''
                last_id = _id
                continue
            _seq += line.upper()
        self._fasta[last_id] = _seq
        sys.stdout.write('\n')

class Fasta3:

    def __init__(self, fasta_path):        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        
        # Progress bar setting
        progress_symbol = '#'
        symbol_num = 50
        each_symbol = line_num / symbol_num
        #
        tmp_file_path = 'biokit.fasta.tmp'
        _tmp = open(tmp_file_path, 'w')
        count = 0
        _seq = ''
        for line in open(fasta_path):
            # Progress bar printing
            count += 1
            symbol_count = int(count / each_symbol) 
            sys.stdout.write('\r')
            sys.stdout.write(('[%-50s]%.2f%%' % ((progress_symbol * symbol_count), 
                               count / line_num * 100)))
            sys.stdout.flush()
            #
            line = line.strip()
            if not line: 
                continue
            if line.startswith('>'):
                _tmp.write('\n')
                _tmp.write(line + '\n')
                continue
            _tmp.write(line.upper())
        _tmp.write('\n')
        _tmp.close()
        sys.stdout.write('\n')

        _tmp = open(tmp_file_path)
        _tmp.readline()
        for line in _tmp:
            line = line.strip()
            if line.startswith('>'):
                _id = line
                continue
            self._fasta[_id] = line.upper()
        _tmp.close()
        os.remove(tmp_file_path)

class Fasta4:

    def __init__(self, fasta_path):        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        
        # Progress bar setting
        progress_symbol = '#'
        symbol_num = 50
        each_symbol = line_num / symbol_num
        #
        count = 0
        for line in open(fasta_path):
            # Progress bar printing
            count += 1
            symbol_count = int(count / each_symbol) 
            sys.stdout.write('\r')
            sys.stdout.write(('[%-50s]%.2f%%' % ((progress_symbol * symbol_count), 
                               count / line_num * 100)))
            sys.stdout.flush()
            #
            line = line.strip()
            if not line: 
                continue
            if line.startswith('>'):
                _id = line
                self._fasta[_id] = []
                continue
            self._fasta[_id].append(line.upper())
        sys.stdout.write('\n')

class Fasta5:

    def __init__(self, fasta_path):        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        
        # Progress bar setting
        progress_symbol = '#'
        symbol_num = 50
        each_symbol = line_num / symbol_num
        #
        count = 0
        _seq = []
        for line in open(fasta_path):
            # Progress bar printing
            count += 1
            symbol_count = int(count / each_symbol) 
            sys.stdout.write('\r')
            sys.stdout.write(('[%-50s]%.2f%%' % ((progress_symbol * symbol_count), 
                               count / line_num * 100)))
            sys.stdout.flush()
            #
            line = line.strip()
            if not line: 
                continue
            if line.startswith('>'):
                _id = line
                if _seq:
                    self._fasta[last_id] = _seq
                    _seq = []
                last_id = _id
                continue
            _seq.append(line.upper())
        self._fasta[last_id] = _seq
        sys.stdout.write('\n')