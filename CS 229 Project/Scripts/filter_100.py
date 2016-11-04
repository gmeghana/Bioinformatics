
# coding: utf-8

# In[12]:

from Bio import SeqIO
folder1 = ["single_error0","single_error1","paired_error0","paired_error0","paired_error1","paired_error1"]
folder2 = ["rep0","rep1","rep2","rep3","rep4"]
ends = ["","","_1","_2","_1","_2"]
i = 0
for f1 in folder1:
    for f2 in folder2:
        input_file = f1+"/"+f2+"/sample_01"+ends[i]+".fasta"
        out_file = f1+"/"+f2+"/sample_01_100"+ends[i]+".fasta"
        file1 = open(input_file,"rU")
        file2 = open(out_file,"w")
        dic_trans ={}
        for reads in SeqIO.parse(file1,"fasta"):
            des,sequence = reads.description,str(reads.seq).strip()
            if len(sequence)>=100:
                file2.write(">"+des+"\n")
                file2.write(sequence+"\n")
        file1.close()
        file2.close()
    i = i+1
    print i


# In[11]:




# In[ ]:



