def progressBar(symbol, symbol_num, count, total):
    symbol_count = int(count / (total / symbol_num))
    sys.stdout.write('\r')
    sys.stdout.write(('[%-50s]%.2f%%' % ((symbol 
                      * symbol_count), count / total * 100)))
    sys.stdout.flush()

from SBM.fasta import *
from SBM.sequence import *
from SBM.utilities import *
print('Simplified Bioinformatics Module (SBM) 1.0')
print('You gonna love using this module!^_^')
