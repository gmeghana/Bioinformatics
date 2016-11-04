number=1000 
total = 176
fold_changes = matrix(c(rep(1,number),rep(1,number)), nrow=number)#total 176241
head(fold_changes)
library(polyester)
library(Biostrings)
fasta_file='/home/xinxin/M229S/polyester-master/inst/extdata/Homo_sapiens.GRCh38.cdna.all.fa'
fasta = readDNAStringSet(fasta_file)
for(i in 1:176){
  begin = (i-1)*number+1
  end= i*number
  small_fasta = fasta[begin:end]
  writeXStringSet(small_fasta, paste("cdna",i,".fa",sep=""))
  readspertx = round(20 * width(small_fasta) / 100)
  simulate_experiment(paste("cdna",i,".fa",sep=""), reads_per_transcript=readspertx, paired = TRUE,error_rate = 0, 
                      num_reps=c(10,10), fold_changes=fold_changes, outdir=paste('/media/My Passport/M229S/paired_error0_cDNA/simulated_reads_paired_error0_',i,sep=""))
}
  
fold_changes = matrix(c(rep(1,241),rep(1,241)), nrow=241)#total 176241
head(fold_changes)

# subset the FASTA file to first 10000 transcripts
small_fasta = fasta[176001:176241]
writeXStringSet(small_fasta, "cdna177.fa")

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)

# simulation call:
simulate_experiment("cdna177.fa", reads_per_transcript=readspertx, paired = TRUE,error_rate = 0, 
                    num_reps=c(10,10), fold_changes=fold_changes, outdir='/media/My Passport/M229S/paired_error0_cDNA/simulated_reads_paired_error0_177') 
