
# step by step running of DIA-MS2pep

# read and parse DIA_MS2pep.config file
paramsfile=$1
shift
filelist="${@}"
mkdir MS2pep
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
ptm_search = $ptm_search
max_processes = $max_processes
datatype = $datatype
" 

# Generating pesudo-ms2 spectra by reading mgf and mzML files

dir=$PWD
filelist=$(sed 's/.raw\|.RAW//g' <<<"$filelist")
echo $filelist

#>$dir/ms2pep.time
for file in $filelist
do
	date "+%T $file pseudo_generation_start " >> $dir/ms2pep.time 
	if [ ! -f ${file}.mzML ]
	then
		echo "DIA-MS2pep read mzML file as input, please use 
		msconvert tool to convert RAW data into mzML!" 
		exit
	fi
	
	if [ ! -f ${file}.mgf ]
	then
		perl $MS2pepCodedir/DataConvert/RAW.mzML.parser.pl mgf ${file}.mzML
	fi
	
	if [ ! -f ${file}.ms1 ]
	then
		perl $MS2pepCodedir/DataConvert/RAW.mzML.parser.pl ms1 ${file}.mzML
	fi
	
	perl $MS2pepCodedir/pseudo-MS2spectrum/DIA_acquistion_window_generator.pl $file.mzML 
	echo $file 
	echo "Generating pseudo-ms2 spectra in $dir/MS2pep"
	#nohup perl $MS2pepCodedir/DIA_pesudo_MS2.pl $file MS2pep $ms2ppm & # MS2pep is default name of the folder for storage of pseudo-ms2 spectra
	if [ $max_processes ]
	then
		if [ $datatype = "DIA" ]
		then
			perl $MS2pepCodedir/pseudo-MS2spectrum/DIA_pesudo_MS2_multiforks.pl $file MS2pep $ms2ppm $max_processes
		else
			perl $MS2pepCodedir/pseudo-MS2spectrum/SWATH_pesudo_MS2_multiforks.pl $file MS2pep $ms2ppm $max_processes
		fi
	else
		if [ $datatype = "DIA" ]
		then
			perl $MS2pepCodedir/pseudo-MS2spectrum/DIA_pesudo_MS2.pl $file MS2pep $ms2ppm $max_processes # MS2pep is default name of the folder for storage of pseudo-ms2 spectra
		else
			perl $MS2pepCodedir/pseudo-MS2spectrum/SWATH_pesudo_MS2.pl $file MS2pep $ms2ppm $max_processes # MS2pep is default name of the folder for storage of pseudo-ms2 spectra
		fi
	fi
	date "+%T $file pseudo_generation_end" >> $dir/ms2pep.time 
done

#wait

if [ -f $MSFragger_search_engine ]
then
	cp $MSFragger_search_engine ${dir}/MS2pep 
else
	echo -e "No $MSFragger_search_engine in $dir \n";
	exit
fi

if [ -f $MSFragger_search_params ]
then
	cp $MSFragger_search_params ${dir}/MS2pep 
else
	echo -e "No $MSFragger_search_engine in $dir \n";
	exit
fi

if [ -f $fasta ]
then
	cp $fasta ${dir}/MS2pep 
else
	echo -e "No $MSFragger_search_engine in $dir  \n";
	exit
fi


cd ${dir}/MS2pep
for file in $filelist
do
	date "+%T $file MSFragger_run_start" >> $dir/ms2pep.time 
	echo "Performing data searching with MSFragger "
	sh $MS2pepCodedir/DataSearch/MSFragger/MSFragger_runner.sh $MSFragger_search_engine $MSFragger_search_params $file $ptm_search $datatype
	date "+%T $file MSFragger_run_end" >> $dir/ms2pep.time 
done

cd $dir
for file in $filelist
do
	echo "Performing data refinement "
	date "+%T $file SearchDataRefinement_start" >> $dir/ms2pep.time 
	if [ $max_processes ]
	then
		if [ $datatype = "DIA" ]
		then
			perl $MS2pepCodedir/SearchDataRefinement/MSFragger/DIA_data_refinement_multiforks.pl $file MS2pep $ms1ppm $ms2ppm $PTM $ptm_search $fasta $max_processes
		else
			perl $MS2pepCodedir/SearchDataRefinement/MSFragger/SWATH_data_refinement_multiforks.pl $file MS2pep $ms1ppm $ms2ppm $PTM $ptm_search $fasta $max_processes
		fi
		
	else
		if [ $datatype = "DIA" ]
		then
			perl $MS2pepCodedir/SearchDataRefinement/MSFragger/DIA_data_refinement.pl $file MS2pep $ms1ppm $ms2ppm $PTM $ptm_search $fasta
		else
			perl $MS2pepCodedir/SearchDataRefinement/MSFragger/SWATH_data_refinement.pl $file MS2pep $ms1ppm $ms2ppm $PTM $ptm_search $fasta
		fi
	fi
	date "+%T $file SearchDataRefinement_end" >> $dir/ms2pep.time 
done

#wait

for file in $filelist
do   date "+%T $file percolater_end" >> $dir/ms2pep.time 
	(
		cd MS2pep; 
		echo "Performing FDR estimation with percolater  "
		if [ $PTM = 1 ]
		then
			level="PTM"
			sh $MS2pepCodedir/FDR/percolator_runner.sh $file $fasta $level
		else
			level="protein"
			sh $MS2pepCodedir/FDR/percolator_runner.sh $file $fasta $level
			if [ $ptm_search = 1 ]
			then
				sh $MS2pepCodedir/FDR/percolator_runner.sh $file $fasta "PTM"
			fi
		fi
		
	)
	date "+%T $file percolater_end" >> $dir/ms2pep.time 
done

