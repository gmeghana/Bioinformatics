
# coding: utf-8

# In[2]:

import random
total_trans = 143838
sub_num = int(total_trans*0.05)
sample_list = []
for i in range(5):
    samples = random.sample(range(1, total_trans), sub_num)
    sample_list.append(samples)
    


# In[3]:

sample_list


# In[4]:

sum(abs(sample_list[1]-sample_list[2]))


# In[12]:

m = 0
for i in sample_list[1]:
    if i in sample_list[0]:
        m = m+1


# In[13]:

print m


# In[7]:

total_trans*0.05


# In[19]:

from Bio import SeqIO
input_file = "D:\\3-2016 spring\\M229S\\Project\\data\\genes.protein_coding_sk.fa"
out_pre = "D:\\3-2016 spring\\M229S\\Project\\data\\genes.protein_coding_sk"
handle = open(input_file, "rU")
transcripts = list(SeqIO.parse(handle, "fasta"))
num = 0
for l in sample_list:
    out_file = out_pre+"_"+str(num)+".fa"
    f1 = open(out_file,"w")
    j = 0
    for i in l:
        des,sequence = transcripts[i].description,transcripts[i].seq.tostring().strip()
        f1.write(">"+des+"\n")
        f1.write(sequence+"\n")
        j = j+1
    f1.close()
    print j
    num+=1
handle.close()


# In[ ]:




# In[ ]:



