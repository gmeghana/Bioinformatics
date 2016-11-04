
# coding: utf-8

# In[30]:

from Bio import SeqIO
input_file = "D:\\3-2016 spring\\M229S\\Project\\data\\Homo_sapiens.GRCh38.cdna.all.fa"
dic_gen = {}
dic_seq = {}
handle = open(input_file, "rU")
i = 0
j = 0
for fasta in SeqIO.parse(handle, "fasta"):
    name, des,sequence = fasta.id,fasta.description,fasta.seq.tostring().strip()
    senquence = sequence.replace(" ", "")
    error = "ACACAGCAGCAAACACAGGAAATGAAAGAGATGTATCAAAATGCAGAAGCTAAAGTGAATAAT"
    
    trans = des.strip().split(' ')
    tid = trans[0][:15]
    gid = trans[3].split(":")[1][:15]
    if error in sequence:
        print "gid:"+gid
        print "tid"+tid
    if gid in dic_gen:
        dic_gen[gid].append(tid)
        dic_seq[gid].append(sequence)
    else:
        duam = {gid:[tid]}
        seq = {gid:[sequence]}
        dic_gen.update(duam)
        dic_seq.update(seq)
        j = j+1
    i = i+1
    if i%10000 == 0:
        print i/10000
print j


# In[32]:

print len(dic_seq["ENSG00000274869"][0])


# In[10]:

out_files = "D:\\3-2016 spring\\M229S\\Project\\data\\cdna84_specialized_format.fa"
genid = dic_gen.keys()
with open(out_files, "w") as out_file:
    for i in range(len(dic_gen)):
        recordID = genid[i]+"|"+"|".join(dic_gen[genid[i]])
        recordSeq = "|".join(dic_seq[genid[i]])
        out_file.write(">"+recordID+"\n")
        out_file.write(recordSeq+"\n")
out_file.close


# In[15]:

print sequence


# In[3]:

from Bio import SeqIO
out_files = "D:\\3-2016 spring\\M229S\\Project\\data\\genes.protein_coding.fa"
handle = open(out_files, "rU")
for fasta in SeqIO.parse(handle, "fasta"):
    name, sequence = fasta.id,fasta.seq.tostring()
    print name
    print sequence
    break


# In[ ]:



