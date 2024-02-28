from Bio import SeqIO
from Bio.Seq import Seq

with open('./dbs/ncbi_amr.faa') as original, open('./dbs/ncbi_amr_corrected.faa', 'w') as corrected:
    for record in SeqIO.parse('./dbs/ncbi_amr.faa', 'fasta'):
        record.seq = Seq(str(record.seq).replace("*", ""))
        SeqIO.write(record, corrected, 'fasta')
