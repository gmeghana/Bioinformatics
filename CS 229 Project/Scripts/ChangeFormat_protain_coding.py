
# coding: utf-8

# In[3]:

from Bio import SeqIO
input_file = "D:\\3-2016 spring\\M229S\\Project\\data\\genes.protein_coding.fa"
out_file = "D:\\3-2016 spring\\M229S\\Project\\data\\genes.protein_coding_sk.fa"
f1 = open(out_file,"w")
handle = open(input_file, "rU")
total_gene = 0
total_trans = 0
for fasta in SeqIO.parse(handle, "fasta"):
    name,sequence = fasta.id,fasta.seq.tostring().strip()
    trans = name.split("|")
    seq = sequence.split("|")
    j = 0
    total_gene +=1
    for i in trans[1:]:
        f1.write(">"+i+" gene:"+trans[0]+"\n")
        f1.write(seq[j]+"\n")
        j = j+1
        total_trans +=1
    if j%10000 == 0:
        print j/10000
handle.close()
f1.close()
        


# In[5]:

total_trans


# In[6]:

total_gene


# In[ ]:



