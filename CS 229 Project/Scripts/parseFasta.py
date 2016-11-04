def fasta_reader(filename):
  from Bio.SeqIO.FastaIO import FastaIterator
  with open(filename) as handle:
    for record in FastaIterator(handle):
      yield record

f = open('groundTruth','w')
for entry in fasta_reader("Homo_sapiens.GRCh38.cdna.all.fa"):
  reads = round(20 * len(str(entry.seq)) / 100)
  print str(entry.id), reads
  f.write("%s %s\n" % (str(entry.id), reads))
f.close()