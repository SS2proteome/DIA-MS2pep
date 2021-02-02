

# step by step running of DIA-MS2pep

# read and parse DIA_MS2pep.params file
paramsfile=$1
while read line
do
	if grep -Pq '^[^#=]+\s*=\s*' <<< $line 
	then
		eval $(sed 's/[[:blank:]]\+\|#.*$//g' <<<$line )
	fi
done < $paramsfile
echo "Params:
dir = $dir
filelist = $filelist
ms1ppm = $ms1ppm
ms2ppm = $ms2ppm 
fasta = $fasta
MSFragger_search_engine = $MSFragger_search_engine
MSFragger_search_params = $MSFragger_search_params
PTM = $PTM 
modification_search = $modification_search
" 

# Generating pesudo-ms2 spectra by reading mgf and mzML files

CodeDir=$PWD
filelist=$(sed 's/.raw\|.RAW//g' <<<$filelist)
for file in $filelist
do
	perl $CodeDir/DIA_acquistion_window_generator.pl $file.mzML 
	echo $file 
	echo "Generating pseudo-ms2 spectra in $dir/MS2pep"
	perl $CodeDir/DIA_pesudo_MS2.pl $file MS2pep $ms2ppm # MS2pep is default name of the folder for storage of pseudo-ms2 spectra
done

test ! -f ${dir}/$MSFragger_search_engine && cp ${CodeDir}/$MSFragger_search_engine ${dir}/MS2pep
test ! -f ${dir}/$MSFragger_search_params && cp ${CodeDir}/$MSFragger_search_params ${dir}/MS2pep
test ! -f ${dir}/$fasta && cp ${CodeDir}/$fasta ${dir}/MS2pep

cd ${dir}/MS2pep
for file in $filelist
do
	echo "Performing data searching with MSFragger "
	sh $CodeDir/MSFragger_runner.sh $MSFragger_search_engine $MSFragger_search_params $file $modification_search
done

cd $dir
for file in $filelist
do
	echo "Performing data refinement "
	perl $CodeDir/DIA_data_refinement.pl $file MS2pep $ms1ppm $ms2ppm $PTM $fasta
	(
	cd MS2pep; 
	echo "Performing FDR estimation with percolater  "
	sh $CodeDir/percolator_runner.sh $file $fasta 0 
	sh $CodeDir/percolator_runner.sh $file $fasta 1 # PTM data
	)
done



