
# coding: utf-8

# In[12]:

from Bio import SeqIO
folder1 = ["single_error0","single_error1","paired_error0","paired_error1"]
folder2 = ["rep0","rep1","rep2","rep3","rep4"]
ends = ["","","_1","_1"]
i = 0
for f1 in folder1:
    for f2 in folder2:
        input_file = f1+"/"+f2+"/sample_01"+ends[i]+".fasta"
        out_file = f1+"/"+f2+"/groundTruth.txt"
        file1 = open(input_file,"rU")
        file2 = open(out_file,"w")
        dic_trans ={}
        for reads in SeqIO.parse(file1,"fasta"):
            name= reads.id
            trans = name.split("/")
            if trans[1] in dic_trans:
                dic_trans[trans[1]]+=1
            else:
                duam = {trans[1]:1}
                dic_trans.update(duam)
        for key,item in dic_trans.iteritems():
            file2.write(key+" "+str(item)+"\n")
        file1.close()
        file2.close()
    i = i+1
    print i


# In[11]:




# In[ ]:



