
# step by step running of DIA-MS2pep

# read and parse DIA_MS2pep.config file
paramsfile=$1
while read line
do
	if grep -Pq '^[^#=]+\s*=\s*' <<< $line 
	then
		eval $(sed 's/[[:blank:]]\+\|#.*$//g' <<<$line )
	fi
done < $paramsfile
echo "Params:
MS2pepCodedir = $MS2pepCodedir
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

dir=$PWD
filelist=$(sed 's/.raw\|.RAW//g' <<<$filelist)
for file in $filelist
do
	if [ ! -f ${file}.mzML ]
	then
		echo "DIA-MS2pep read mzML file as input, please use 
		msconvert tool to convert RAW data into mzML!" 
		exit
	fi
	
	if [ ! -f ${file}.mgf ]
	then
		perl $MS2pepCodedir/RAW.mzML.parser.pl mgf ${file}.mzML
	fi
	
	if [ ! -f ${file}.ms1 ]
	then
		perl $MS2pepCodedir/RAW.mzML.parser.pl ms1 ${file}.mzML
	fi
	
	perl $MS2pepCodedir/DIA_acquistion_window_generator.pl $file.mzML 
	echo $file 
	echo "Generating pseudo-ms2 spectra in $dir/MS2pep"
	perl $MS2pepCodedir/DIA_pesudo_MS2.pl $file MS2pep $ms2ppm # MS2pep is default name of the folder for storage of pseudo-ms2 spectra
done

test ! -f ${dir}/$MSFragger_search_engine && cp ${MS2pepCodedir}/$MSFragger_search_engine ${dir}/MS2pep
test ! -f ${dir}/$MSFragger_search_params && cp ${MS2pepCodedir}/$MSFragger_search_params ${dir}/MS2pep
test ! -f ${dir}/$fasta && cp ${MS2pepCodedir}/$fasta ${dir}/MS2pep

cd ${dir}/MS2pep
for file in $filelist
do
	echo "Performing data searching with MSFragger "
	sh $MS2pepCodedir/MSFragger_runner.sh $MSFragger_search_engine $MSFragger_search_params $file $modification_search
done

cd $dir
for file in $filelist
do
	echo "Performing data refinement "
	perl $MS2pepCodedir/DIA_data_refinement.pl $file MS2pep $ms1ppm $ms2ppm $PTM $fasta
	(
	cd MS2pep; 
	echo "Performing FDR estimation with percolater  "
	sh $MS2pepCodedir/percolator_runner.sh $file $fasta 0 
	sh $MS2pepCodedir/percolator_runner.sh $file $fasta 1 # PTM data
	)
done
