import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

input_table = sys.argv[1]

with open(input_table) as f:
  for line in f:
    file, gene, assembly_name, assembly_path = line.strip().split(',')[0:4]
    
    # read names and postions from bed file
    positions = defaultdict(list)
    with open(file) as f:
        for line in f:
            name, start, stop = line.split()[0:3]
            strand = line.split()[5]
            positions[name].append((int(start), int(stop), strand))
            
    # parse fasta file and turn into dictionary
    records = SeqIO.to_dict(SeqIO.parse(open(assembly_path), 'fasta'))
    
    # search for short sequences
    short_seq_records = []
    for name in positions:
        for (start, stop, strand) in positions[name]:
            long_seq_record = records[name]
            long_seq = long_seq_record.seq
            short_seq = str(long_seq)[start:stop]
            # reverse complement the sequence if gene is on minus strand. 
            if strand=='-': 
              short_seq_record = SeqRecord(Seq(short_seq).reverse_complement(), id=name, description='')
            else: 
              short_seq_record = SeqRecord(Seq(short_seq), id=name, description='')
            short_seq_records.append(short_seq_record)
            
    # write to file
    with open('exons_'+gene+'_'+assembly_name+'.fasta', 'w') as f:
      SeqIO.write(short_seq_records, f, 'fasta')
