file=$1
fasta=$2
level=$3

if ((${#@} < 3 ))
then
	echo "Percolator.run.sh [file] [fasta] [leve: psm/peptide/protein/PTM] ";
	exit
fi


if [ $level = "PTM" ]
then
	percolator.exe $file.PTM.pin --results-peptides $file.PTM.pin.target.pep.tsv --decoy-results-peptides $file.PTM.pin.decoy.pep.tsv  -Y -N 200000 #-f $fasta
elif [ $level = "peptide" ]
then
	percolator.exe $file.pin --results-peptides $file.pin.target.pep.tsv  --decoy-results-peptides $file.pin.decoy.pep.tsv  -Y -N 200000   
elif [ $level = "psm" ]
then
	percolator.exe $file.pin --results-peptides $file.pin.target.pep.tsv  --decoy-results-peptides $file.pin.decoy.pep.tsv  -Y -N 200000 --results-psms $file.pin.target.psm.tsv --decoy-results-psms $file.pin.decoy.psm.tsv   
else
	percolator.exe $file.pin --results-peptides $file.pin.target.pep.tsv --results-proteins $file.pin.protein.tsv -A --decoy-results-peptides $file.pin.decoy.pep.tsv -P 'REV_' -Y -N 200000 -f $fasta  --results-psms $file.pin.target.psm.tsv --decoy-results-psms $file.pin.decoy.psm.tsv
fi

