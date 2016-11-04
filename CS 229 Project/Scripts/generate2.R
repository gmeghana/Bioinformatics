number=7191
fold_changes = matrix(c(rep(1,number)),nrow = number)
head(fold_changes)

library(polyester)
library(Biostrings)

# FASTA annotation
fasta_file='/home/xinxin/M229S/polyester-master/inst/extdata/genes.protein_coding_sk_4.fa'
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first 20 transcripts

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(fasta) / 100)

# simulation call:
simulate_experiment(fasta_file, reads_per_transcript=readspertx, paired = TRUE,error_rate = 0, 
                    num_reps=c(1), fold_changes=fold_changes, outdir='/media/My Passport/M229S/paired_error1/rep4')
