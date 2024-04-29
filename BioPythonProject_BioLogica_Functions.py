from Bio.SeqUtils import GC
from Bio import Align
import sys
from sys import argv
import getopt
from itertools import product
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
###############################################
# Function 1
def gc(seq):
    for i in seq:
        if i not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")
    my_seq = Seq(seq)
    gc = GC(my_seq)
    print(f"GC_Percentage = {gc} %")
###############################################
# Function 2
def transcribe(dna):
    for i in dna:
        if i not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")
    basecomplement = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 'u', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(dna)
    letters = [basecomplement[base] for base in letters]
    print(''.join(letters))
###############################################
# Function 3
def reverse_complement(seq):
    for i in seq:
        if i not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")
    reversed = seq[::-1]
    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(reversed)
    letters = [base_complement[base] for base in letters]
    print(''.join(letters))
###############################################
# Function 4
def calc_nbases(dna):
    for i in dna:
        if i not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")
    nbases = dna.count('N') + dna.count('n')
    return nbases
###############################################
# Function 5
def is_valid(seq,type):
  Seq=seq.upper()
  bool=True
  if type=='protein':
    for i in range(len(Seq)):
        if Seq[i] not in 'ABCDEFGHIKLMNPQRSTVWXYZ':
            bool=False
            break
  elif type=='dna':
    for i in range(len(Seq)):
        if Seq[i] not in 'ACGTN':
            bool=False
            break
  elif type=='rna':
    for i in range(len(Seq)):
        if Seq[i] not in 'ACGUN':
            bool=False
            break
  else:
    print("Please enter a valid type")
    exit()

  return  bool
###############################################
# Function 6
def filter_nbases(dna):
    for i in dna:
        if i not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")
    filtered_dna = ""
    for i in dna:
        if (i == 'n' or i == 'N'):
            continue
        else:
            filtered_dna+= i
    return filtered_dna
###############################################
# Function 7
def seq_alginment(str1, str2):
    for i in str1:
        if i not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")
    for j in str2:
        if j not in 'ACGTNacgtn':
            raise Exception("Please enter a dna sequence contains (a,c,g,t,A,C,G,T)")

    The_Alignment = Align.PairwiseAligner()
    Seq_Align = The_Alignment.align(str1, str2)
    Seq_Align_Sorted = sorted(Seq_Align)
    for align in Seq_Align_Sorted :
        The_Seq_Alignment = align
        Score = align.score
    print( f"Sequence Alignment : \n{The_Seq_Alignment}\nScore = {Score}")
###############################################
# Function 8
def seq_alignment_files (file1,file2):
    seqs1 = SeqIO.to_dict(SeqIO.parse(open(file1), 'fasta'))
    seqs2 = SeqIO.to_dict(SeqIO.parse(open(file2), 'fasta'))
    for sr1, sr2 in product(seqs1, seqs2):
        for a in pairwise2.align.localxx(str(seqs1[sr1].seq), str(seqs2[sr2].seq)):
            print(format_alignment(*a))
###############################################
# Function 9
def online_alignment(sequence, rest):
    opts, argv = getopt.getopt(rest, 'o:')
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    print("result from searching online....")
    print(result_handle)
    blast_records = NCBIXML.parse(result_handle)
    blast_records = list(blast_records)
    print(blast_records)
    E_VALUE_THRESH = 0.00000000001
    count = 0
    for opt, arg in opts:
        if opt in ['-o']:
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            count += 1
                            print("*Alignment*")
                            print("sequence:", alignment.title)
                            print("length:", alignment.length)
                            print(hsp.query[0:75] + "...")
                            print(hsp.match[0:75] + "...")
                            print(hsp.sbjct[0:75] + "...")
                            with open(arg, 'a') as save_file:
                                save_file.write(alignment.title)
                                save_file.write(str(alignment.length))
                                save_file.writelines(hsp.query[0:75] + "...")
                                save_file.writelines(hsp.match[0:75] + "...")
                                save_file.writelines(hsp.sbjct[0:75] + "...")
                            print()
    print(f"There are {count} similar sequences in Blast output")
###############################################
# Function 10
def merge_fasta(parameters):
    files = ""
    argv = sys.argv[2:]
    if len(argv)<2:
        print("Please enter at least two fasta files to merge.")
        exit()
    if (argv[-2] == "-o"):
        if len(parameters)-2 < 2:
            print("Please enter at least two fasta files to merge.")
            exit()
        else:
            for i in range(len(parameters) - 2):
                records = list(SeqIO.parse(parameters[i], "fasta"))
                dna = records[0]
                desc = dna.description
                seq = dna.seq
                line1 = ">" + desc
                files += line1
                files += "\n"
                i = 0
                j = 70
                while (True):
                    files += seq[i:j]
                    files += "\n"
                    i += 70
                    j += 70
                    if (j > len(seq)):
                        files += seq[i:]
                        files += "\n"
                        break
            argv = sys.argv[-2:]
            opts, args = getopt.getopt(argv, "o:")
            for opt, arg in opts:
                if opt == '-o':
                    f2 = open(arg, "a")
                    f2.write(str(files))
                    f2.close()
                    print("Files merged")

    else:
        for i in parameters:
            records = list(SeqIO.parse(i, "fasta"))
            dna = records[0]
            desc = dna.description
            seq = dna.seq
            line1 = ">" + desc
            files += line1
            files += "\n"
            i = 0
            j = 70
            while (True):
                files += seq[i:j]
                files += "\n"
                i += 70
                j += 70
                if (j > len(seq)):
                    files += seq[i:]
                    files += "\n"
                    break
        print(files)
###############################################
# Function 11
def convert_to_fasta(file):
    with open(file) as input_handle,open("ls_orchid.fasta","w") as output_handle :
        genbank=SeqIO.parse(input_handle,"genbank")
        fasta = SeqIO.write(genbank, output_handle, "fasta")
        print("Converted %i records" % fasta)
###############################################
#running functions
if (sys.argv[1] == "gc"):
    if (len(sys.argv)<2):
        print("Please enter the sequence")
        exit()
    else:
        gc(sys.argv[2])

elif (sys.argv[1] == "transcribe"):
    if (len(sys.argv)<2):
        print("Please enter the sequence")
        exit()
    else:
        transcribe(sys.argv[2])

elif (sys.argv[1] == "reverse_complement"):
    if (len(sys.argv)<2):
        print("Please enter the sequence")
        exit()
    else:
        reverse_complement(sys.argv[2])

elif (sys.argv[1] == "calc_nbases"):
    if (len(sys.argv)<2):
        print("Please enter the sequence")
        exit()
    else:
        print(calc_nbases(sys.argv[2]))

elif (sys.argv[1] == "is_valid"):
    if (len(sys.argv)<3):
        print("Please enter the sequence and it's type")
        exit()
    else:
        print(is_valid(sys.argv[2], sys.argv[3]))

elif (sys.argv[1] == "filter_nbases"):
    if (len(sys.argv)<2):
        print("Please enter the sequence")
        exit()
    else:
        print(filter_nbases(sys.argv[2]))

elif (sys.argv[1] == "seq_alignment"):
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv[-2:], 'o:')
    if len(opts) != 0:
        if len(argv) == 5:
            if opts[0][0] == '-o':
                sys.stdout = open(opts[0][1], "w")
        else:
            print("Please check your parameters")
            exit()
    opt = argv[0]
    seq_alginment(argv[1], argv[2])

elif (sys.argv[1] == "seq_alignment_files"):
    argv = sys.argv[1:]
    options = 'o:'
    opts, args = getopt.getopt(argv[-2:], options)
    if len(opts) != 0:
        if len(argv) == 5:
            if opts[0][0] == '-o':
                sys.stdout = open(opts[0][1], "w")
        else:
            print("Please check your parameters")
            exit()

    opt = argv[0]
    seq_alignment_files(argv[1], argv[2])

elif(sys.argv[1] == "online_alignment"):
    online_alignment(sys.argv[2], sys.argv[3:])

elif(sys.argv[1] == "merge_fasta"):
    merge_fasta(argv[2:])

elif(sys.argv[1] == "convert_to_fasta"):
    if (len(sys.argv) < 2):
        print("Please enter the genbank file")
        exit()
    else:
        convert_to_fasta(sys.argv[2])

else:
    print("Please enter a valid function name.")