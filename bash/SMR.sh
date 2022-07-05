#!/bin/bash

#gwas=$1
gwas=F1.gz
dis=${gwas//.gz}
#mkdir /lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/$dis
#cd /lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/$dis

#zcat ../$gwas | csvcut -d " " -c CHR,SNP,A1,A2,p,b,se,p,N | csvformat -T | sed "1d" | awk '{print $2,$3,$4,$5,$6,$7,$8,$9 > ""$1".ma"}'
#zcat $gwas | csvcut  -c CHR,SNP,A1,A2,p,b,se,p,N | csvformat -T | sed "1d" | awk '{print $2,$3,$4,$5,$6,$7,$8,$9 > ""$1".ma"}'
#zcat $gwas > int | sed -i 's/\t/,/g' int | gzip > $gwas
for chr in {1..22}
do
sed -i "1i SNP A1 A2 Freq b se p n" $chr.ma 
done
#cd /lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/

for tissue in /lustre/home/acct-bmelgn/bmelgn-3/SMR_format/eqtl/*/;
#for tissue in /lustre/home/acct-bmelgn/bmelgn-3/SMR_format/eqtl/adultbrain/;
do

for chr in {1..22} ;
do
/lustre/home/acct-bmelgn/bmelgn-3/smr_Linux  \
--bfile /lustre/home/acct-bmelgn/bmelgn-3/ldsc/partition/1000G.EUR.QC.$chr \
--gwas-summary $chr.ma \
--beqtl-summary $tissue$chr.besd \
--out $tissue$chr.$dis \
--smr-multi \
--diff-freq 1 \
--peqtl-smr 1e-5 
done

cat $tissue*$dis.msmr > $tissue$dis
rm $tissue*$dis.*

done
#rm -R ./$dis/
