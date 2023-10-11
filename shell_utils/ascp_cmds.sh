aspera_link="era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR525/004/SRR5259954/SRR5259954_1.fastq.gz" # need to copy from ENA
ascp -v -QT -l 400m -P33001 -k1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh $aspera_link

