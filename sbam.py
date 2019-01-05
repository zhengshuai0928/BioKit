def progressBar(symbol, symbol_num, count, total):
    symbol_count = int(count / (total / symbol_num))
    sys.stdout.write('\r')
    sys.stdout.write(('[%-50s]%.2f%%' % ((symbol 
                      * symbol_count), count / total * 100)))
    sys.stdout.flush()

from SBAM.fasta import *
from SBAM.sequence import *
from SBAM.utilities import *