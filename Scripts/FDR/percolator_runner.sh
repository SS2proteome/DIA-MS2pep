file=$1
fasta=$2
ptm=$3
if ((ptm==1))
then
percolator.exe $file.PTM.tsv --results-peptides $file.PTM.pin.target.pep.tsv -A --decoy-results-peptides $file.PTM.pin.decoy.pep.tsv -P 'REV_' -Y -N 200000 -f $fasta   
else
percolator.exe $file.pin --results-peptides $file.pin.target.pep.tsv --results-proteins $file.pin.protein.tsv -A --decoy-results-peptides $file.pin.decoy.pep.tsv -P 'REV_' -Y -N 200000 -f $fasta   
fi
