import sys
import getopt

import pandas as pd
import phylopandas as ph
import csv
import re
import os

class Select():
    """A class to select informations about some scaffold
        
    Returns:
        A fasta file and two .gtf file corresponding to scaffold
    """

    def __init__(self):
        self.scaffold_number = None
        self.dir = None

    def fasta_file(self, fasta_path, scaffold_number):
        """This method select scaffold's fasta sequence

        Args:
            fasta_path ([str]): [fasta's path]
            scaffold_number ([int]): [scaffold's number]
        """

        self.scaffold_number = scaffold_number
        self.scaffold_number = str(self.scaffold_number)
        fasta = ph.read_fasta(fasta_path)

        if 'scaffold_' + scaffold_number in fasta.id.values:

            self.dir = 'scaffold_' + scaffold_number + '_info'
            os.mkdir(self.dir)

            fasta.drop(columns={'description', 'name'}, inplace=True)
            to_save = fasta[fasta.id == 'scaffold_' + self.scaffold_number].copy()
            to_save.id = '>' + to_save.id
            to_save.to_csv('./'+ self.dir + '/scaffold_' + self.scaffold_number +'.fasta', sep='\n', index=False, header=False, quoting=csv.QUOTE_NONE)

        else:
            print('Scaffold ' + self.scaffold_number + ' nao existente!')
            exit()

    def TE_file(self, te_path):
        """This method select transposable elements in the scaffold

        Args:
            te_path ([str]): [transposable element's path]
        """
    
        self.scaffold_number = str(self.scaffold_number)
        
        TEs = pd.read_csv(te_path, sep='\t',header=None, comment='#')
        TEs.rename(columns={0:'seqId'}, inplace=True)
        TEs.drop(columns=9, inplace=True)
        
        to_save = TEs[TEs.seqId == 'scaffold_' + self.scaffold_number+'_c1'].copy()
        to_save['seqId'] = to_save['seqId'].str.extract(r'(scaffold_\d+)')[0]
        to_save['type'] = 'repeat_region'
        to_save.to_csv('./'+ self.dir + '/TEs_scaffold_' + self.scaffold_number +'.gff', sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    def genoma_file(self, gen_path):
        """This method select the scaffold's regions in genome

        Args:
            gen_path ([str]): [genome's path]
        """
    
        self.scaffold_number = str(self.scaffold_number)
        genoma = pd.read_csv(gen_path, sep='\t',header=None, comment='#')
        genoma.rename(columns={0:'seqId'}, inplace=True)
        
        to_save = genoma[genoma.seqId == 'scaffold_' + self.scaffold_number].copy()
        to_save.to_csv('./'+ self.dir + '/scaffold_' + self.scaffold_number +'_secreted.gff', sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
    
    try:
        OPTS, ARGS = getopt.getopt(sys.argv[1:], 'f:t:g:n:h', ['fasta_path', 
                                                             'te_path', 
                                                             'gen_path', 
                                                             'scaffold_number'
                                                             'help'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)

    SELECT = Select()
    FASTA_PATH = None
    TE_PATH = None
    GEN_PATH = None
    SCAFFOLD_NUMBER = None

    for opt, arg in OPTS:

        if opt in ('-h', '--help'):
            print('''
            This program select all informations about a scaffold. 
            select_scaffold.py -f --fasta_path -t --te_path -g --gen_path -n --scaffold_number -h --help
           ''')
            sys.exit(2)
        
        elif opt in ('-f', '--fasta_path'):
            if re.match('.+fasta$|.+txt$', arg):
                FASTA_PATH = arg
            else:
                print("-f isn't a fasta or txt file")
                sys.exit(3)

        elif opt in ('-t', '--te_path'):
            if re.match('.+gff$', arg):
                TE_PATH = arg
            else:
                print("-t isn't a gff or txt file")
                sys.exit(4)

        elif opt in ('-g', '--gen_path'):
            if re.match('.+gff$', arg):
                GEN_PATH = arg
            else:
                print("-g isn't a gff or txt file")
                sys.exit(5)

        elif opt in ('-n', '--scaffold_number'):
            if re.match('[0-9]+$', arg):
                SCAFFOLD_NUMBER = arg
            else:
                print("-n isn't a number")
                sys.exit(6)

    if FASTA_PATH and TE_PATH and GEN_PATH and SCAFFOLD_NUMBER:

        SELECT.fasta_file(FASTA_PATH, SCAFFOLD_NUMBER)
        SELECT.TE_file(TE_PATH)
        SELECT.genoma_file(GEN_PATH)

    else:
        print('missing arguments -f, -t, -g or -n')
        sys.exit(7)