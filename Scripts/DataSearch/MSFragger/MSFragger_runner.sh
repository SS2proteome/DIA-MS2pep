# MSFragger runner

MSFragger=$1
params=$2
file=$3
ptm_search=$4
datatype=$5

lower=-150
upper=500





if [ $ptm_search  = 1 ]
then
	if [ $datatype = "DIA" ]
	then
		charge_range=$(echo {1..5})
		for wmgf in w_[0-9]*_${file}.mgf
		do
		echo $wmgf
		filesize=$(stat -c "%s" $wmgf)
		if [[ $filesize -gt 2000000000 ]]  # MSFragger can not handle mgf with the size over 2G 
		then
			echo "The size of $mgf is too large for MSFragger, so split it into several files"
			awk -vx=1 -vRS="BEGIN" 'NR>1{if(NR%100000 == 0)x++;print RS $0 > "temp_"x"_"FILENAME}' $wmgf
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				echo $z
				awk -vw=$window -vz=$z -vlower=$lower -vupper=$upper ' /precursor_mass_lower =/{print "precursor_mass_lower = "lower;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = "upper;next}
				/precursor_charge =/{print "precursor_charge = "z" "z;next}
				/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params temp_[0-9]*$wmgf 2>/dev/null 
				awk -vp=1 '/<spectrum_query/{p=0;t=1}p;/<\/msms_run_summary>/{t=0};t;END{print "</msms_run_summary>\n</msms_pipeline_analysis>"}' temp_[0-9]_${wmgf/.mgf/}.pepXML > ${wmgf/.mgf/}.pepXML  
				mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
				
			done
			rm -rf temp_[0-9]*$mgf
			rm -rf temp_[0-9]*${mgf/.mgf/}.pepXML
		else
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				echo $z
				awk -vw=$window -vz=$z -vlower=$lower -vupper=$upper ' /precursor_mass_lower =/{print "precursor_mass_lower = "lower;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = "upper;next}
				/precursor_charge =/{print "precursor_charge = "z" "z;next}/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params $wmgf 2>/dev/null
				mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
			done
		fi
		done
	else
		charge_range=$(echo {2..5})
		for wmgf in w_[0-9]*_${file}.mgf
		do
		echo $wmgf
		filesize=$(stat -c "%s" $wmgf)
		if [[ $filesize -gt 2000000000 ]]  # MSFragger can not handle mgf with the size over 2G 
		then
			echo "The size of $mgf is too large for MSFragger, so split it into several files"
			awk -vx=1 -vRS="BEGIN" 'NR>1{if(NR%100000 == 0)x++;print RS $0 > "temp_"x"_"FILENAME}' $wmgf
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				echo $z
				awk -vw=$window -vz=$z -vlower=$lower -vupper=$upper ' /precursor_mass_lower =/{print "precursor_mass_lower = "lower;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = "upper;next}
				/precursor_charge =/{print "precursor_charge = "z" "z;next}
				/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				
				awk -vx=$z '/TITLE/{$0=$0 x RS "CHARGE=" x "+"}/^0 0/{next}{print $0 > "z"x"."FILENAME}' temp_[0-9]*$wmgf 
				
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params z$z.temp_[0-9]*$wmgf  
				awk -vp=1 '/<spectrum_query/{p=0;t=1}p;/<\/msms_run_summary>/{t=0};t;END{print "</msms_run_summary>\n</msms_pipeline_analysis>"}'  z$z.temp_[0-9]*_${wmgf/.mgf/}.pepXML > z$z.{wmgf/.mgf/}.pepXML  
				#mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
				
			done
			rm -rf z*.temp_[0-9]*$mgf
			rm -rf z*.temp_[0-9]*${mgf/.mgf/}.pepXML
		else
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				awk -vx=$z '/TITLE/{$0=$0 x RS "CHARGE=" x "+"}/^0 0/{next}1' $wmgf > z$z.$wmgf
				echo $z
				awk -vw=$window -vz=$z -vlower=$lower -vupper=$upper ' /precursor_mass_lower =/{print "precursor_mass_lower = "lower;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = "upper;next}
				/precursor_charge =/{print "precursor_charge = "z" "z;next}/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params z$z.$wmgf  
				#mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
			done
			rm -rf z*.$wmgf
		fi
		done
	fi
	
else
	if [ $datatype = "DIA" ]
	then
		charge_range=$(echo {1..5})
		for wmgf in w_[0-9]*_${file}.mgf
		do
		echo $wmgf
		filesize=$(stat -c "%s" $wmgf)
		if [[ $filesize -gt 1000000000 ]]  # MSFragger can not handle mgf with the size over 2G 
		then
			echo "The size of $mgf is too large for MSFragger, so split it into several files"
			awk -vx=1 -vRS="BEGIN" 'NR>1{if(NR%100000 == 0)x++;print RS $0 > "temp_"x"_"FILENAME}' $wmgf
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				echo $z
				awk -vw=$window -vz=$z '
				/precursor_mass_lower =/{print "precursor_mass_lower = -" z*w/2 + 3;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = " z*w/2 + 3;next} 
				/precursor_charge =/{print "precursor_charge = "z" "z;next}
				/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params

				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params temp_[0-9]*$wmgf 
				
				awk -vp=1 '/<spectrum_query/{p=0;t=1}p;/<\/msms_run_summary>/{t=0};t;END{print "</msms_run_summary>\n</msms_pipeline_analysis>"}' temp_[0-9]_${wmgf/.mgf/}.pepXML > ${wmgf/.mgf/}.pepXML  
				mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
				
			done
			rm -rf temp_[0-9]*$mgf
			rm -rf temp_[0-9]*${mgf/.mgf/}.pepXML
		else
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				echo $z
				awk -vw=$window -vz=$z '/precursor_mass_lower =/{print "precursor_mass_lower = -" z*w/2 + 3;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = " z*w/2 + 3;next} 
				/precursor_charge =/{print "precursor_charge = "z" "z;next}
				/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params $wmgf 
				
				mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
			done
		fi
		done
	else
		charge_range=$(echo {2..5})
		for wmgf in w_[0-9]*_${file}.mgf
		do
		echo $wmgf
		filesize=$(stat -c "%s" $wmgf)
		if [[ $filesize -gt 2000000000 ]]  # MSFragger can not handle mgf with the size over 2G 
		then
			echo "The size of $mgf is too large for MSFragger, so split it into several files"
			awk -vx=1 -vRS="BEGIN" 'NR>1{if(NR%100000 == 0)x++;print RS $0 > "temp_"x"_"FILENAME}' $wmgf
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				echo $z
				awk -vw=$window -vz=$z '
				/precursor_mass_lower =/{print "precursor_mass_lower = -" z*w/2 + 3;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = " z*w/2 + 3;next} 
				/precursor_charge =/{print "precursor_charge = "z" "z;next}
				/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				
				awk -vx=$z '/TITLE/{$0=$0 x RS "CHARGE=" x "+"}/^0 0/{next}{print $0 > "z"x"."FILENAME}' temp_[0-9]*$wmgf 
				
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params z$z.temp_[0-9]*$wmgf  
				awk -vp=1 '/<spectrum_query/{p=0;t=1}p;/<\/msms_run_summary>/{t=0};t;END{print "</msms_run_summary>\n</msms_pipeline_analysis>"}'  z$z.temp_[0-9]*_${wmgf/.mgf/}.pepXML > z$z.{wmgf/.mgf/}.pepXML  
				#mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
				
			done
			rm -rf z*.temp_[0-9]*$mgf
			rm -rf z*.temp_[0-9]*${mgf/.mgf/}.pepXML
		else
			for z in $charge_range
			do
				window=$(grep -Po '(?<=w_)\d+(?=_)' <<<$wmgf)
				awk -vx=$z '/TITLE/{$0=$0 x RS "CHARGE=" x "+"}/^0 0/{next}1' $wmgf > z$z.$wmgf
				echo $z
				awk -vw=$window -vz=$z '
				/precursor_mass_lower =/{print "precursor_mass_lower = -" z*w/2 + 3;next}
				/precursor_mass_upper =/{print "precursor_mass_upper = " z*w/2 + 3;next} 
				/precursor_charge =/{print "precursor_charge = "z" "z;next}
				/override_charge =/{print "override_charge = 1";next}1' $params > z$z.$window.$params
				java -d64 -Xmx30G -jar $MSFragger z$z.$window.$params z$z.$wmgf  
				#mv ${wmgf/mgf/pepXML} z$z.${wmgf/mgf/pepXML}
			done
			rm -rf z*.$wmgf
		fi
		done
	fi
fi
