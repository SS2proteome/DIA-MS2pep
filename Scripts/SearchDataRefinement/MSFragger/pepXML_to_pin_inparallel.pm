package pepXML_to_pin_inparallel;

use strict;
use Data::Dumper;
use File::Basename;
use Math::Combinatorics qw(combine);
use lib dirname(__FILE__);
use POSIX qw(floor ceil round);
use List::Util qw(max min sum);
use peptide;
use IsotopeDistribution;
use MassAccuracy;

$| = 1;


my $filename;
my $dir;
my $outfilename;
my $logfile;
my $ms1_ppm_tolerence; #ppm
my $ms2_ppm_tolerence; #ppm
my $ptm;
my $ptmfind;

my $unimod_info = readmodxml();

my %mono_aa;
$mono_aa{'A'} = 71.037114;
$mono_aa{'R'} = 156.101111;
$mono_aa{'N'} = 114.042927;
$mono_aa{'D'} = 115.026943;
$mono_aa{'C'} = 103.009185;
#$mono_aa{'C_'} = 160.030649; # Carboxyamid
$mono_aa{'C'} = 103.009185;
$mono_aa{'E'} = 129.042593;
$mono_aa{'Q'} = 128.058578;
$mono_aa{'G'} = 57.021464;
$mono_aa{'H'} = 137.058912;
$mono_aa{'I'} = 113.084064;
$mono_aa{'L'} = 113.084064;
$mono_aa{'K'} = 128.094963;
$mono_aa{'M'} = 131.040485;
$mono_aa{'M_'} = 131.040485;
#$mono_aa{'M_'} = 131.040485 + $O; # Oxidation
$mono_aa{'F'} = 147.068414;
$mono_aa{'P'} = 97.052764;
$mono_aa{'S'} = 87.032028;
#$mono_aa{'S_'} = 166.9984; # Phos
$mono_aa{'T'} = 101.047679; 
#$mono_aa{'T_'} = 181.0140; # Phos
$mono_aa{'U'} = 150.95363;
$mono_aa{'W'} = 186.079313;
$mono_aa{'Y'} = 163.06332;
#$mono_aa{'Y_'} = 243.0297; # Phos
$mono_aa{'V'} = 99.068414;


my %aa_short = (
"Ala" => "A",
"Arg" => "R",
"Asn" => "N",
"Asp" => "D",
"Cys" => "C",
"Gln" => "Q",
"Glu" => "E",
"Gly" => "G",
"His" => "H",
"Leu" => "L",
"Lys" => "K",
"Met" => "M",
"Phe" => "F",
"Pro" => "P",
"Ser" => "S",
"Thr" => "T",
"Trp" => "W",
"Tyr" => "Y",
"Val" => "V",
"Xle" => "I|L"
);

my $min_n_iso_required = 4;
my $proton = 1.00727646677;
my $massdiff_C12_C13 = 1.0033548378;
my $max_peakdepth = 10;

my $ms1_file;
my $mzML_file;
my $mgf_file;
my @mgffiles;
my @pepxmlfiles;
my $window_file;

my $mgffile;
my $pepxmlfile;

my @pre_ppms;
my @frag_ppms;
my $scan;
my $ms1scan;
my $ms2scan;
my %ms2toms1;
my $index_ms1;
my %index_ms1scan;
my %index_ms2scan;

my $ms1scan_mz;
my %ms1scan_mz_low_resolution;

my %DIA_window;
my %isolationlist;
my $num_of_windows;


my %ms1scan_start;
my %ms1scan_cnt_index;
my %ms1scan_basepeak;
my %ms2scan_last_index;
my %ms2scan_next_index;
my %premz_index_ms2scan;
my $mslevel;
my $premz_;

my @ms1_charges;
my %charge_assaignment;
my %ms1_peaks;
my %ms1_file_position;

my %ms1scan_scan_index;
my $ms1scan_cnt;

my %ms2_file_position;
my %pseudomgffilehandles;
my %title_to_mgffile;
my %pseudomgffilehandle_position;


my $end_scan;
my $retention_time_sec;
my $precursor_intensity;
my $peptide_prev_aa;
my $peptide_next_aa;
my $num_matched_ions;
my $tot_num_ions;
my $num_missed_cleavages;
my $calc_neutral_pep_mass;
my $protein;
my $start_scan;
my $start_scan_temp;
my $assumed_charge;
my $precursor_neutral_mass;
my $spectrum;
my $peptide;
my $strip_peptide;
my $massdiff;
my $mod_aminoacid_mass;
my $mod_aminoacid_position;
my %mod_aminoacid;
my $modified_peptide;
my $hyperscore;
my $nextscore;
my $expect;
my $search_data_line;
my $pepLen;
my $hit_rank;
my %search_res;
my %multiple_hits_rm;
my %multiple_peptides;
my %ms1scan_clean_data;
my %peptide_charge_cnt;
my %peptide_expect;
my %keep_specs;
my %keep_pep_list;
my %keep_pep_list_expect_strip_pep;
my %iso_distribution_buffer;
my %factorial_buffer;

my %n;
my %uniq_pl;
my %multiple_hits;

my %PTM_Types;
my @PTM_variable;
my %fixed_variable;
my %pseudoms2_data;
my %deisotope_pseudoms2_data;
my %ms2_data;
my ($ms1profile_,$ms2profile_);

my %pepxml_filehandle_position;
my %pepxml_filehandle;
my %pepxml_scan_data;
my @pepxml_start_scan;

my $n;
my $tell_pos_last;
my $fastafile;

sub pepXML_refinement{
	#print Dumper $_[0];
	$filename = $_[0];
	$dir = $_[1];
	#$outfilename = "$filename/${filename}_pseudo.mgf";
	#$logfile = "${filename}.run.log";
	$ms1_ppm_tolerence = $_[2]; #ppm
	$ms2_ppm_tolerence = $_[3]; #ppm
	$ptm = $_[4];
	$ptmfind = $_[5];
	$fastafile = $_[6];
	my $index_ = $_[7];
	
	my $window_info = $_[8];
	%DIA_window = %{$window_info->[0]};
	%isolationlist = %{$window_info->[1]};
	$num_of_windows = $window_info->[2];
	
	my $mzML_info = $_[9];
	%ms1scan_start = %{$mzML_info->[0]};
	%ms1scan_cnt_index = %{$mzML_info->[1]};
	%ms1scan_basepeak = %{$mzML_info->[2]};
	%ms2scan_last_index = %{$mzML_info->[3]};
	%ms2scan_next_index = %{$mzML_info->[4]};
	%premz_index_ms2scan = %{$mzML_info->[5]};
	%ms2toms1 = %{$mzML_info->[6]};
	
	
	$ms1_file = "$dir/$filename/index${index_}_$filename.ms1";
	$mgf_file = "$dir/$filename/index${index_}_$filename.mgf";

	@mgffiles = glob "$dir/$filename/pseudo_index${index_}_w_[0-9]*_".$filename.".mgf";
	@pepxmlfiles = glob "$dir/$filename/index${index_}_[z][0-9].w_[0-9]*_".$filename.".pepXML";
	#print Dumper \@pepxmlfiles,\@mgffiles,$ms1_file,$mgf_file;

	
	open(MS1,$ms1_file) or die "$! $ms1_file\n";
	while(<MS1>){
		chomp;
		if(/^S/){
			($ms1scan) = /S\t(\d+)\t\d+/; 
			#$ms1scan_cnt++;
			my $pos = tell(MS1);
			$ms1_file_position{sprintf("%06.f",$ms1scan)} = $pos;
			#printf "\r$ms1scan";
			}
	}
	
	open(MS2,$mgf_file) or die "$! $mgf_file\n";
	while(<MS2>){
		chomp;
		if(/^TITLE/){
			($ms2scan) = /.+\.(\d+)\.\1/; 
			my $pos = tell(MS2);
			$ms2_file_position{sprintf("%06.f",$ms2scan)} = $pos;
			#printf "\r$ms2scan";
		}
	}

	foreach my $pseudomgffile(@mgffiles){
		#print $pseudomgffile,"\n";
		open($pseudomgffilehandles{$pseudomgffile},"$pseudomgffile") or die "no $pseudomgffile $!\n";
		my $fh = $pseudomgffilehandles{$pseudomgffile};
		while(<$fh>){
			chomp;
			if(/^TITLE/){
				/TITLE=(.+\.(\d+)\.\2)/; 
				my $title = $1;
				my $pos = tell($fh);
				$pseudomgffilehandle_position{$title} = $pos;
				$title_to_mgffile{$title} = $pseudomgffile;
				#printf "\r$title";
			}
		}
	}

		
	foreach my $pepxmlfile(@pepxmlfiles){
		#print $pepxmlfile,"\n";
		$n++;
		open($pepxml_filehandle{$pepxmlfile},"$pepxmlfile") or die "$pepxmlfile $!\n";
		my $fh = $pepxml_filehandle{$pepxmlfile};
		#print "\n$pepxmlfile\n";
		while(<$fh>){
			chomp;
			my $line = $_;
			#<aminoacid_modification aminoacid="C" massdiff="57.0215" mass="160.0307" variable="N"/>
			if(m{<aminoacid_modification aminoacid="([^"]+)" massdiff="([^"]+)" mass="([^"]+)" variable="Y"/>}){
				#$modifications{$1."[".$3."]"} = $2;
				if($n == 1){
					push @PTM_variable, $1."[".$3."]";
				}
			}elsif(m{<aminoacid_modification aminoacid="([^"]+)" massdiff="([^"]+)" mass="([^"]+)" variable="N"/>}){
				#$modifications{$1."[".$3."]"} = $2;
				if($n == 1){
					$fixed_variable{$1} = $2;
				}
			}elsif(m{<terminal_modification massdiff="[^"]+" protein_terminus="([^"]+)" mass="([^"]+)" terminus="([^"]+)" variable="Y"/>}){
				if($n == 1){
					if($1 eq "Y"){
						if($3 eq "N"){
							push @PTM_variable, "n(".$2.")";
						}elsif($3 eq "C"){
							push @PTM_variable, "(".$2.")c";
						}
					}elsif($1 eq "N"){	
						if($3 eq "N"){
							push @PTM_variable, "n(".$2.")";
						}elsif($3 eq "C"){
							push @PTM_variable, "(".$2.")c";
						}
					}	
				}
			}
			
			if(m{<spectrum_query .*spectrum=}){
				($start_scan,$assumed_charge,$spectrum,$end_scan,$precursor_neutral_mass,$retention_time_sec ) = ($line =~ m{<spectrum_query start_scan="(\d+)".*assumed_charge="(\d+)".*spectrum="([^"]+)".*end_scan="(\d+)".*index="\d+".*precursor_neutral_mass="([^"]+)".*retention_time_sec="([^"]+)"});
				$spectrum =~ /([^.]+\.(\d+)\.\2).*/;
				$start_scan = $2;
				printf "\r$start_scan";
				my $pos = $tell_pos_last;
				#push @pepxml_start_scan,$start_scan;
				push @{$pepxml_scan_data{$start_scan}},[$pepxmlfile,$spectrum];
				$pepxml_filehandle_position{$pepxmlfile}{$spectrum} = $pos;
			}	
			$tell_pos_last = tell($fh);
		}
	}


	@pepxml_start_scan = sort{$a <=> $b}(keys %pepxml_scan_data);
	

	my %buffer_ms1;
	open(TMPOUT,">${dir}/$filename/${index_}_masscal_tmp.$filename");
	
	#print TMPOUT Dumper @pepxmlfiles, @pepxml_start_scan;

	foreach my $scan_list (@pepxml_start_scan){
		
		my $prems1scan = $ms2toms1{$scan_list}->[0];
		if(! $buffer_ms1{$prems1scan}++){
			undef %pseudoms2_data;
			undef %deisotope_pseudoms2_data;
			undef %ms2_data;
			undef %ms1scan_clean_data;
		}
		foreach my $pepxml_info (@{$pepxml_scan_data{$scan_list}}){
			my ($pepxmlfile,$spectrumlist) = @{$pepxml_info};
			my $pepxml_pos = $pepxml_filehandle_position{$pepxmlfile}{$spectrumlist};
			my $pepxml_fh = $pepxml_filehandle{$pepxmlfile};
			
			seek($pepxml_fh,$pepxml_pos,0);
			while(<$pepxml_fh>){
				chomp;
				if(m{<spectrum_query .*spectrum=}){
					#($start_scan,$assumed_charge,$spectrum,$precursor_neutral_mass) = ($line =~ m{<spectrum_query start_scan="(\d+)" assumed_charge="(\d+)" spectrum="([^"]+)" end_scan="\d+" index="\d+" precursor_neutral_mass="([^"]+)" retention_time_sec="[^"]+"});
					($start_scan,$assumed_charge,$spectrum,$end_scan,$precursor_neutral_mass,$retention_time_sec ) = (m{<spectrum_query start_scan="(\d+)".*assumed_charge="(\d+)".*spectrum="([^"]+)".*end_scan="(\d+)".*index="\d+".*precursor_neutral_mass="([^"]+)".*retention_time_sec="([^"]+)"});
					$spectrum =~ s/([^.]+\.(\d+)\.\2).*/\1/;
					$start_scan = $2;
					undef $ms1profile_;
					undef $ms2profile_;
					
					
				}elsif(m{<search_hit .* hit_rank=}){
				#($peptide,$massdiff,$calc_neutral_pep_mass,$protein) = ($line =~ m{peptide="([^"]+)" massdiff="([^"]+)" calc_neutral_pep_mass="([^"]+)" peptide_next_aa="[^"]+" num_missed_cleavages="\d+" num_tol_term="\d+" num_tot_proteins="\d+" tot_num_ions="[^"]+" hit_rank="[^"]+" num_matched_ions="[^"]+" protein="([^"]+)" peptide_prev_aa="[^"]+" is_rejected="\d+">});
				($peptide,$massdiff,$calc_neutral_pep_mass,$peptide_next_aa,$num_missed_cleavages,$tot_num_ions,$hit_rank,$num_matched_ions,$protein,$peptide_prev_aa) = (m{peptide="([^"]+)" massdiff="([^"]+)" calc_neutral_pep_mass="([^"]+)" peptide_next_aa="(.)" num_missed_cleavages="(\d+)" num_tol_term="\d+" num_tot_proteins="\d+" tot_num_ions="([^"]+)" hit_rank="([^"]+)" num_matched_ions="([^"]+)" protein="([^"]+)" peptide_prev_aa="(.)" is_rejected="\d+">});
				$pepLen = length($peptide);
				$strip_peptide = $peptide;
				}elsif(m{<mod_aminoacid_mass mass="([^"]+)" position="(\d+)"/>}){
					$mod_aminoacid_mass = $1;
					$mod_aminoacid_position = $2;
					$mod_aminoacid{$mod_aminoacid_position} = $mod_aminoacid_mass;
					
				}elsif(m{<modification_info (mod_nterm_mass|mod_cterm_mass)="([^"]+)">}){
					$mod_aminoacid_mass = $2;
					$mod_aminoacid_position = $1;
					$mod_aminoacid{$mod_aminoacid_position} = $mod_aminoacid_mass;
				}elsif(m{</modification_info>}){
					
					my @pep = ($peptide =~ /./g);
					$peptide = join "",map{$pep[$_].($mod_aminoacid{$_+1}? "[" :"" ).$mod_aminoacid{$_+1}.($mod_aminoacid{$_+1}? "]" : "")}0..$#pep;
					
					if(exists $mod_aminoacid{mod_nterm_mass}){
						$peptide = "n(".$mod_aminoacid{mod_nterm_mass}.")".$peptide;
					}
					if(exists $mod_aminoacid{mod_cterm_mass}){
						$peptide = $peptide . "(".$mod_aminoacid{mod_cterm_mass}.")c";
					}
					
					while($peptide =~ /.\[[^][]+\]/g){
						my $ptm_ = $&;
						$PTM_Types{$ptm_}{$peptide_prev_aa.".".$peptide.".".$peptide_next_aa} = 1;
					}
					
					if(exists $mod_aminoacid{mod_nterm_mass}){
						my $ptm_ = "n(".$mod_aminoacid{mod_nterm_mass}.")";
						$PTM_Types{$ptm_}{$peptide_prev_aa.".".$peptide.".".$peptide_next_aa} = 1;
					}
					if(exists $mod_aminoacid{mod_cterm_mass}){
						my $ptm_ = "(".$mod_aminoacid{mod_cterm_mass}.")c";
						$PTM_Types{$ptm_}{$peptide_prev_aa.".".$peptide.".".$peptide_next_aa} = 1;
					}

					$strip_peptide = $peptide;
					$strip_peptide =~ s/\[[^][]+\]|c\([^()]+\)|n\([^()]+\)//g;
					undef %mod_aminoacid;
				}
				elsif(m{<search_score name="hyperscore" value="([^"]+)"/>}){
					$hyperscore = $1;
				}
				if(m{<search_score name="nextscore" value="([^"]+)"/>}){
					$nextscore = $1;
				}
				elsif(m{<search_score name="expect" value="([^"]+)"/>}){
					$expect = $1;
					if(length($peptide) and $strip_peptide =~ /^[A-Z_]+$/ ){# and $start_scan =~ /^108741$/ ){#){ 
						printf "\r$start_scan";
						#next if $start_scan < 20000;
						#exit if $start_scan > 25000;
						my $exp_premz;
						my $exp_int;
						my $exp_SN;
						my $curr_ppm_obs;
						my $back_ppm_obs;
						my $cal_premz;
						$cal_premz = ($calc_neutral_pep_mass + $assumed_charge * $proton )/$assumed_charge;

						my $Current_Iso_tag;
						my $Current_exp_premz;
						my $Current_exp_int;
						my $pass = 0;
						my $ios_topic_peaks;
						my $ios_topic_mz;
						my $ios_topic_int;
						my $ios_topic_ppm;
						my $ios_corr;
						my $curr_SN;
						my $back_corr;
						
						my $prems1scan = $ms2toms1{$start_scan}->[0];
						($Current_Iso_tag,$Current_exp_premz,$Current_exp_int,$curr_ppm_obs,$curr_SN) = iso_peak_check($prems1scan,$cal_premz,$assumed_charge,$min_n_iso_required,$start_scan,$peptide);
						$ios_topic_peaks = $Current_Iso_tag;
						$ios_topic_mz = $Current_exp_premz;
						$ios_topic_int = $Current_exp_int;
						$ios_topic_ppm = $curr_ppm_obs;
						$exp_SN = $curr_SN;


						my $next_tag = 0;
						if($ios_topic_peaks == 0){
							my $next_prems1scan = $index_ms1scan{$index_ms2scan{$start_scan} + 1};
							my ($next1scan_Iso_tag,$next_exp_premz,$next_exp_int,$next_ppm_obs,$next_SN);
							($next1scan_Iso_tag,$next_exp_premz,$next_exp_int,$next_ppm_obs,$next_SN) = iso_peak_check($next_prems1scan,$cal_premz,$assumed_charge,$min_n_iso_required,$start_scan,$peptide);
							if($ios_topic_peaks < $next1scan_Iso_tag){
								$ios_topic_peaks = $next1scan_Iso_tag;
								$ios_topic_mz = $next_exp_premz;
								$ios_topic_int = $next_exp_int;
								$ios_topic_ppm = $next_ppm_obs;
								$exp_SN = $next_SN;
								$next_tag = 1;
							}elsif($ios_topic_peaks == $next1scan_Iso_tag){
								if($ios_topic_int < $next_exp_int){
									$ios_topic_peaks = $next1scan_Iso_tag;
									$ios_topic_mz = $next_exp_premz;
									$ios_topic_int = $next_exp_int;					
									$ios_topic_ppm = $next_ppm_obs;
									$exp_SN = $next_SN;
									$next_tag = 1;
								}
							}
						}
						$exp_premz = $ios_topic_mz;
						$exp_int = $ios_topic_int;
						
						my ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$matched_y,$matched_b,$matched_int_sd,$ppm_m,$ppm_s,$matched_int,$matched_coverage,$matched_mod_type,$matched_int_with_mod,$mod_SN,$matched_mzs_type,$match_type,$position,$mod_name,$unimod_mass,$mod_peptide,$PTMSiteScore,$PTMSiteProb,$new_calmz);

						my $tag_ = $ios_topic_peaks ? "Iso" : "Mod"; 


						if( $ptmfind && $ios_topic_peaks == 0 && $expect < 1 && $next_tag == 0 and !$multiple_hits{$tag_}{$start_scan}{$protein}++){

							$ms1profile_ = defined $ms1profile_ ? $ms1profile_ : MS1_profile($prems1scan,$start_scan);
							$ms2profile_ = defined $ms2profile_ ? $ms2profile_ : MS2_profile($start_scan);

							($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$matched_y,$matched_b,$matched_int_sd,$ppm_m,$matched_int,$matched_coverage,$matched_mod_type,$matched_int_with_mod,$mod_SN,$position,$mod_name,$unimod_mass,$mod_peptide,$PTMSiteScore,$PTMSiteProb,$new_calmz,$ppm_s) = find_mod($prems1scan,$cal_premz,$assumed_charge,$start_scan,$peptide,$pepLen,$strip_peptide,$spectrum,$ms1profile_,$ms2profile_,$peptide_prev_aa,$peptide_next_aa);

							if($mod_Iso_tag){
								$exp_premz = $mod_exp_premz;
								$exp_int = $mod_exp_int;
								$cal_premz = $new_calmz;
								$exp_SN = $mod_SN;
								$calc_neutral_pep_mass = $cal_premz * $assumed_charge  - $assumed_charge * $proton;
								$peptide = $mod_peptide;
								$ios_topic_ppm = ppm($exp_premz,$cal_premz);

							}
						}

						if ($ios_topic_peaks and !$multiple_hits{$tag_}{$start_scan}{$protein}++){
							my $peaks ;
							
							if(exists $pseudoms2_data{$spectrum}){
								$peaks = $pseudoms2_data{$spectrum};
							}else{
								$peaks = readpseudoMS2($spectrum);
								$pseudoms2_data{$spectrum} = $peaks;
							}
							my $theoretical_mz_table = peptide::fragmentation($peptide,$assumed_charge);
							($matched_mzs_type,$matched_int,$matched_int_sd,$ppm_m,$ppm_s) = ms2spectra_match($theoretical_mz_table,$peaks);
							($match_type,$matched_coverage) = tag($matched_mzs_type);
							($matched_y,$matched_b) = ($match_type->{"Y"},$match_type->{"B"});
						}
						
						# print  join "\t",
						# $expect,
								# $spectrum,
								# $assumed_charge,
								# $cal_premz,
								# $strip_peptide,
								# $protein,
								# $ios_topic_ppm,
								# $mod_Iso_tag,
								# $cal_premz,$exp_int,$ios_topic_peaks,$peptide,$expect,$protein,$strip_peptide,$peptide_prev_aa.".".$peptide.".".$peptide_next_aa,#$centermz,
								
									# $spectrum.".".$assumed_charge.".".$hit_rank,($protein =~ /REV_/ ? -1 : 1),$start_scan,$ios_topic_peaks ? $precursor_neutral_mass : $calc_neutral_pep_mass,$calc_neutral_pep_mass,$hyperscore,($hyperscore- $nextscore),log($expect),$num_matched_ions/$tot_num_ions,$precursor_neutral_mass,$pepLen,$assumed_charge == 1 ? 1 :0, $assumed_charge == 2 ? 1 :0, $assumed_charge == 3 ? 1 :0, $assumed_charge == 4 ? 1 :0, $assumed_charge == 5 ? 1 :0, $assumed_charge == 6 ? 1 :0,$num_missed_cleavages,($precursor_neutral_mass - $calc_neutral_pep_mass),abs($precursor_neutral_mass - $calc_neutral_pep_mass), $ios_topic_peaks == 1 ? 1 : 0,$ios_topic_peaks == 2 ? 1 : 0, $ios_topic_peaks == 3 ? 1 : 0, $ios_topic_peaks >= 4 ? 1 : 0, log($exp_int+1),$retention_time_sec,$exp_SN
								# ,
								# $matched_y,$matched_b,$matched_int_sd,$ppm_m,$matched_int,$matched_coverage,$mod_name,$PTMSiteScore,$PTMSiteProb,$unimod_mass,$ppm_s
						# ,"\n";

						if ($expect <= 0.01 && !$uniq_pl{PreCheckPPM}{$peptide}++){
							$n{$filename}{PreCheckPPM}{$protein=~/REV_/?"Decoy":"Target"}++ ;
							$n{$filename}{PreCheck_Isotopic}{$protein=~/REV_/?"Decoy":"Target"}{$ios_topic_peaks}++ ;
						}
						if( ($matched_coverage )&& ($mod_Iso_tag || $ios_topic_peaks)
							){
							
							if ($expect <= 0.01 && !$uniq_pl{CheckPPM}{$peptide}++){
								$n{$filename}{CheckPPM}{$protein=~/REV_/?"Decoy":"Target"}++ ;
								$n{$filename}{CheckPPM_Isotopic}{$protein=~/REV_/?"Decoy":"Target"}{$ios_topic_peaks}++; 
							}
							my $ppm = $ios_topic_ppm;
							if ($expect <= 0.01 && $protein !~ /REV_/ && $ios_topic_peaks ){
								push @pre_ppms, $ppm ;
								push @frag_ppms, $ppm_m ;
							}
							my $centermz = $ms2toms1{$start_scan}->[1];
							
							if($mod_Iso_tag > 0){
								push @{$multiple_hits_rm{$exp_premz}{2}},[
								$peptide,
								log($matched_int_with_mod + 1) + $matched_coverage
								];
							}
							if($ios_topic_peaks > 0){
								push @{$multiple_hits_rm{$exp_premz}{1}},[
								$peptide,
								$expect
								];
							}
							$precursor_neutral_mass = $assumed_charge * $exp_premz - $assumed_charge * $proton;
							$ios_topic_peaks = $ios_topic_peaks ? $ios_topic_peaks : $mod_Iso_tag;
							$multiple_peptides{$mod_Iso_tag > 0 ? 2 : 1}{$peptide} = [
								$expect,
								$spectrum,
								$assumed_charge,
								$exp_premz,
								$strip_peptide,
								$protein,
								$ppm,
								
								[$cal_premz,$exp_int,$ios_topic_peaks,$peptide,$expect,$protein,$strip_peptide,$peptide_prev_aa.".".$peptide.".".$peptide_next_aa,$centermz],
								[
									$spectrum.".".$assumed_charge.".".$hit_rank,($protein =~ /REV_/ ? -1 : 1),$start_scan,$ios_topic_peaks ? $precursor_neutral_mass : $calc_neutral_pep_mass,$calc_neutral_pep_mass,$hyperscore,($hyperscore- $nextscore),log($expect),$num_matched_ions/$tot_num_ions,$precursor_neutral_mass,$pepLen,$assumed_charge == 1 ? 1 :0, $assumed_charge == 2 ? 1 :0, $assumed_charge == 3 ? 1 :0, $assumed_charge == 4 ? 1 :0, $assumed_charge == 5 ? 1 :0, $assumed_charge == 6 ? 1 :0,$num_missed_cleavages,($precursor_neutral_mass - $calc_neutral_pep_mass),abs($precursor_neutral_mass - $calc_neutral_pep_mass), $ios_topic_peaks == 1 ? 1 : 0,$ios_topic_peaks == 2 ? 1 : 0, $ios_topic_peaks == 3 ? 1 : 0, $ios_topic_peaks >= 4 ? 1 : 0, log($exp_int+1),$retention_time_sec,$exp_SN
								],
								[$matched_y,$matched_b,$matched_int_sd,$ppm_m,$matched_int,$matched_coverage,$mod_name,$PTMSiteScore,$PTMSiteProb,$unimod_mass,$ppm_s,$ppm]
							];
							
						
						}
					}
				}elsif(m{/spectrum_query}){
					
					foreach my $exp_mz (keys %multiple_hits_rm){
						if(defined $multiple_hits_rm{$exp_mz}{1} && ! defined $multiple_hits_rm{$exp_mz}{2}){
							my @sort = map{$_->[0]}(sort{$a->[1] <=> $b->[1]}@{$multiple_hits_rm{$exp_mz}{1}});
							my $keep_peptide = $sort[0];
							$search_res{$start_scan}{$keep_peptide} = $multiple_peptides{1}{$keep_peptide};
						}elsif(defined $multiple_hits_rm{$exp_mz}{1} && defined $multiple_hits_rm{$exp_mz}{2}){
							my @sort1 = map{$_->[0]}(sort{$a->[1] <=> $b->[1]}@{$multiple_hits_rm{$exp_mz}{1}});
							my $keep_peptide1 = $sort1[0];
							my $keep_peptide1_matched_int = $multiple_peptides{1}{$keep_peptide1}->[9]->[4];
							my $keep_peptide1_matched_coverage = $multiple_peptides{1}{$keep_peptide1}->[9]->[5];
							my $keep_strip_pep1 = $multiple_peptides{1}{$keep_peptide1}->[4];
							
							my @sort2 = map{$_->[0]}(sort{$a->[1] <=> $b->[1]}@{$multiple_hits_rm{$exp_mz}{2}});
							my $keep_peptide2 = $sort2[0];
							my $keep_peptide2_matched_int = $multiple_peptides{2}{$keep_peptide2}->[9]->[4];
							my $keep_peptide2_matched_coverage = $multiple_peptides{2}{$keep_peptide2}->[9]->[5];
							my $keep_strip_pep2 = $multiple_peptides{2}{$keep_peptide2}->[4];
							if($keep_peptide1_matched_int < $keep_peptide2_matched_int && $keep_peptide1_matched_coverage < $keep_peptide2_matched_coverage && 
								! ($keep_strip_pep1 =~ $keep_strip_pep2 || $keep_strip_pep2 =~ $keep_strip_pep1)
							){
								$search_res{$start_scan}{$keep_peptide2} = $multiple_peptides{2}{$keep_peptide2};
							}else{
								$search_res{$start_scan}{$keep_peptide1} = $multiple_peptides{1}{$keep_peptide1};
							}
							
						}else{
							my @sort = map{$_->[0]}(sort{$b->[1] <=> $a->[1]}@{$multiple_hits_rm{$exp_mz}{2}});
							my $keep_peptide = $sort[0];
							$search_res{$start_scan}{$keep_peptide} = $multiple_peptides{2}{$keep_peptide};
						}
					}
					
					undef %multiple_hits_rm;
					undef %multiple_peptides;
					undef %multiple_hits;
					
					last;
				}
			}	
		}		
		foreach my $pep (keys %{$search_res{$start_scan}}){
			push @{$keep_pep_list{$start_scan}},$pep;
			my $d = $search_res{$start_scan}{$pep};
			$keep_pep_list_expect_strip_pep{$pep} = [$d->[0],$d->[4]];
			print TMPOUT join "\t", $start_scan,$pep,map{my $x = $_; ref($x) eq "ARRAY" ? (join "#",map{$_}@{$x}) : $x}@{$d};
			print TMPOUT "\n";
		}
		
		undef %search_res;
	}			
				
	foreach my $pepxmlfile(@pepxmlfiles){
		my $fh = $pepxml_filehandle{$pepxmlfile};			
		close($fh);	
	}		

	close(TMPOUT);
	#print Dumper \%n;

}


sub median{
	my $ref = shift;
	my @sorted = sort{$a <=> $b} @{$ref};
	my $len = scalar @sorted;
	my $median;
	if ($len % 2){
		$median = $sorted[($len-1)/2];
	}else{
		$median = ($sorted[($len-2)/2] + $sorted[$len/2])/2;
	}
	return $median;
}

sub ave{ # Average calculation
	my $data = shift;
	if($#$data>=0){
		my ($sum);
		map {$sum += $_}@{$data};
		my $m = $sum/(1 + $#$data);
		return $m;
	}else{
		return "NA";
	}
}

sub stdev{ # Standard deviation calculation
	my ($data,$ave_) = (@_);
	if($#$data>=0){
		my $delta_sum;
		map {$delta_sum += ($_ - $ave_)**2}@{$data};
		return sqrt($delta_sum/(1 + $#$data));
	}else{
		return "NA";
	}
}

sub seq_similar{
	my ($seq1,$seq2) = (@_);
	my $len1 = length($seq1);
	my $len2 = length($seq2);
	my $diff;
	if(abs($len1 - $len2) > 1){
		return 0;
	}elsif($len1 == $len2){
		for(1..$len1){
			if(substr($seq1,$_-1,1) ne substr($seq2,$_-1,1)){
				$diff++;
			}
		}
		if($diff <= 2){
			return 1;
		}else{
			return 0;
		}
	}else{
		my $commCnt;
		for(1..$len1){
			if(substr($seq1,$_-1,1) eq substr($seq2,$_-1,1)){
				$commCnt++;
			}else{
				last;
			}
		}
		for($len1..1){
			if(substr($seq1,$_-1,1) eq substr($seq2,$_-1,1)){
				$commCnt++;
			}else{
				last;
			}
		}
		if(($len1 - $commCnt <= 2) && ($len2 - $commCnt <= 2)){
			return 1;
		}else{
			return 0;
		}
	}
}


sub ppm{
	my ($observed_mz,$theoretical_mz) = (@_);
	return 1000000*($observed_mz-$theoretical_mz)/$theoretical_mz;
}


sub readmodxml(){
	my %mods;
	my @mod_mass;
	my $mod_title;
	my $mod_delta;
	my %neutral_loss;
	my $site;
	my $position;
	my %uniq;
	my @tmp;
	my $dir = dirname(__FILE__);
	open(MOD,"${dir}/unimod.xml") or die "No ${dir}/unimod.xml, Please download the file of \"umimod.xml\" from http://www.unimod.org/downloads.html $!\n";
	while(<MOD>){
		if(m{<umod:mod title="([^"]+)"}){
			$mod_title = $1;
			$mod_title =~ s/-\&gt;/->/g;
		}
		if(m{<umod:specificity hidden="[^"]+" site="([^"]+)" position="([^"]+)"}){
			$site = $1;
			$position = $2;
			push @tmp,[$site,$position];
		}
		if(m{<umod:NeutralLoss mono_mass="([^"]+)"}){
			$neutral_loss{$site}{$position} = $1 if $1 != 0;
			
		}
		if(m{<umod:delta mono_mass="([^"]+)"}){
			$mod_delta = $1;
			push @mod_mass,$mod_delta if !$uniq{$mod_delta}++;
			foreach my $i (@tmp){
				($site,$position) = @{$i};
				push @{$mods{$mod_delta}{$site}{$position}},[$neutral_loss{$site}{$position},$mod_title];
			}
			undef @tmp;
			undef %neutral_loss;

		}
	}
	close(MOD);
	return [\%mods,\@mod_mass];
}


sub MS1_profile{
	my ($curr_ms1scan,$curr_ms2scan) = (@_);
	my $curve_dot_num = 5;
	my $com_mz;
	my ($prems1scan,$ms2scan);
	my ($curr_ms1scan_mz_list,$curr_lowest_sig,$curr_highest_sig);
	
	if(exists $ms1scan_clean_data{$curr_ms1scan}{$curr_ms2scan}){
		($curr_ms1scan_mz_list,$curr_lowest_sig,$curr_highest_sig) = @{$ms1scan_clean_data{$curr_ms1scan}{$curr_ms2scan}};
	}else{
		($curr_ms1scan_mz_list,$curr_lowest_sig,$curr_highest_sig) = MS1file($curr_ms1scan,$curr_ms2scan);
		$ms1scan_clean_data{$curr_ms1scan}{$curr_ms2scan} =[$curr_ms1scan_mz_list,$curr_lowest_sig,$curr_highest_sig];
	}
	
	my @ms1_data_list;
	my ($ms1scan_mz_list,$lowest_sig,$highest_sig);
	# look forward
	for(1..2){
		if($_ == 1){
			$ms2scan = $ms2scan_last_index{$curr_ms2scan};
			$prems1scan = $ms2toms1{$ms2scan}->[0];
		}else{
			$ms2scan = $ms2scan_last_index{$ms2scan};
			$prems1scan = $ms2toms1{$ms2scan}->[0];
		}
		
		if(exists $ms1scan_clean_data{$prems1scan}{$ms2scan}){
			($ms1scan_mz_list,$lowest_sig,$highest_sig) = @{$ms1scan_clean_data{$prems1scan}{$ms2scan}};
		}else{
			($ms1scan_mz_list,$lowest_sig,$highest_sig) = MS1file($prems1scan,$ms2scan);
			$ms1scan_clean_data{$prems1scan}{$ms2scan} =[$ms1scan_mz_list,$lowest_sig,$highest_sig];
		}
		push @ms1_data_list,$ms1scan_mz_list;
	}
	@ms1_data_list = reverse(@ms1_data_list);
	push @ms1_data_list,$curr_ms1scan_mz_list;
	# look afterward
	for(1..2){
		if($_ == 1){
			$ms2scan = $ms2scan_next_index{$curr_ms2scan};
		}else{
			$ms2scan = $ms2scan_next_index{$ms2scan};
		}
		$prems1scan = $ms2toms1{$ms2scan}->[0];
		if(exists $ms1scan_clean_data{$prems1scan}{$ms2scan}){
			($ms1scan_mz_list,$lowest_sig,$highest_sig) = @{$ms1scan_clean_data{$prems1scan}{$ms2scan}};
		}else{
			($ms1scan_mz_list,$lowest_sig,$highest_sig) = MS1file($prems1scan,$ms2scan);
			$ms1scan_clean_data{$prems1scan}{$ms2scan} =[$ms1scan_mz_list,$lowest_sig,$highest_sig];
		}
		push @ms1_data_list,$ms1scan_mz_list;
	}
	
	foreach my $cmp_spec (@ms1_data_list){
		$com_mz = Spectrum_cmp($curr_ms1scan_mz_list,$cmp_spec,$com_mz);
	}
	return $com_mz;
	
}

sub MS2_profile{
	my ($currms2scan) = (@_);
	my $ms2scan;
	my $ms2peaks;
	my $com_mz;
	my @ms2_data_list;
	
	my $currpeaks;
	#my $peaks;
	if(exists $ms2_data{sprintf("%06.f",$currms2scan)}){
		$currpeaks = $ms2_data{sprintf("%06.f",$currms2scan)};
	}else{
		$currpeaks = readMS2(sprintf("%06.f",$currms2scan));
		$ms2_data{sprintf("%06.f",$currms2scan)} = $currpeaks;
	}
	
	for(1..2){
		if($_ == 1){
			$ms2scan = $ms2scan_last_index{$currms2scan};
		}else{
			$ms2scan = $ms2scan_last_index{$ms2scan};
		}
		$ms2peaks = readMS2(sprintf("%06.f",$ms2scan));
		push @ms2_data_list,$ms2peaks;
	}
	@ms2_data_list = reverse(@ms2_data_list);
	push @ms2_data_list, $currpeaks;
	for(1..2){
		if($_ == 1){
			$ms2scan = $ms2scan_next_index{$currms2scan};
		}else{
			$ms2scan = $ms2scan_next_index{$ms2scan};
		}
		$ms2peaks = readMS2(sprintf("%06.f",$ms2scan));
		push @ms2_data_list,$ms2peaks;
	}
	foreach my $cmp_spec (@ms2_data_list){
		$com_mz = Spectrum_cmp($currpeaks,$cmp_spec,$com_mz);
	}
	
	return $com_mz;
}

sub boundary_index{
    my $data = shift;
	my $bp = shift;
	my $m = $#$data/2;
	
	my ($start,$end);
	
	for($m+1..$#$data){
		if($data->[$_] < 0.05){
			$end = $_;
			last;
		}
	}
	
	for(0..$m-1){
		if($data->[$_] < 0.05){
			$start = $_;
		}
	}
	$start = defined $start ? $start  : 0;
    $end = defined $end ? $end : $#$data;
    return ($start,$end);
}

sub find_mod(){
	my ($prems1scan,$cal_premz,$assumed_charge,$ms2scan,$pep,$pepLen,$strip_pep,$spectrum,$ms1profile,$ms2profile,$peptide_prev_aa,$peptide_next_aa) = (@_);
	my $peaks ;
	my $deisotope_peaks; 
	
	if(exists $pseudoms2_data{$spectrum}){
		$peaks = $pseudoms2_data{$spectrum};
	}else{
		$peaks = readpseudoMS2($spectrum);
		$pseudoms2_data{$spectrum} = $peaks;
	}
	if(exists $deisotope_pseudoms2_data{$spectrum}){
		$deisotope_peaks = $deisotope_pseudoms2_data{$spectrum};
	}else{
		$deisotope_peaks = deisotope($peaks);
		$deisotope_pseudoms2_data{$spectrum} = $deisotope_peaks;
	}
	
	my ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$mod_pos,$position,$NL,$mod_name,$unimod_mass);
	my ($PTMSiteScore,$PTMSiteProb);
	# print STDERR join "\t","\ntest:","prems1scan:$prems1scan,cal_premz:$cal_premz,assumed_charge:$assumed_charge,ms2scan:$ms2scan,pep:$pep,pepLen:$pepLen,strip_pep:$strip_pep","\n";
	
	my ($ms1scan_mz_list,$lowest_sig,$highest_sig);
	if(exists $ms1scan_clean_data{$prems1scan}{$ms2scan}){
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = @{$ms1scan_clean_data{$prems1scan}{$ms2scan}};
	}else{
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = MS1file($prems1scan,$ms2scan);
		$ms1scan_clean_data{$prems1scan}{$ms2scan} =[$ms1scan_mz_list,$lowest_sig,$highest_sig];
	}
	
	
	my $len = $#$ms1scan_mz_list;
	my $obs_mono_mz;
	my $obs_mono_int;
	my $tag_iso;
	my $measured_ppm_mono;
	my $corr_distribution;
	
	my $ms1_start_mz = $ms1scan_mz_list->[0]->[0];
	my $ms1_end_mz = $ms1scan_mz_list->[$len]->[0];
	if(! defined $ms1_start_mz){
		return (0);
		next;
	}
	
	# print STDERR "MS1 start:$ms1_start_mz  MS1 end:$ms1_end_mz \n";
	my $theoretical_mz_table = peptide::fragmentation($pep,$assumed_charge);
	my ($matched_mzs_type,$matched_int_without_mod,$matched_int_sd,$m,$match_mzs_topN5) = ms2spectra_match($theoretical_mz_table,$peaks);
	my ($match_tag,$match_coverage_without_mod) = tag($matched_mzs_type);
	
	my @theoretical_iso_distribution; 
	if(exists $iso_distribution_buffer{$strip_pep}){
		@theoretical_iso_distribution = @{$iso_distribution_buffer{$strip_pep}};
	}else{
		@theoretical_iso_distribution = @{IsotopeDistribution::CAPTT($strip_pep)};
		$iso_distribution_buffer{$strip_pep} =[@theoretical_iso_distribution];
	}
	my $len_of_matched_y_without_mod = $match_tag->{"Y"};
	my $len_of_matched_b_without_mod = $match_tag->{"B"};
	
	my $index = 0;
	my %trace_mz;
	my $iso_require = 10;
	my $match_coverage_with_mod;# = $match_coverage;
	my $mod_matched_coverage;# = $mod_match_coverage;
	my $len_of_matched_y_with_mod = $len_of_matched_y_without_mod;
	my $len_of_matched_b_with_mod = $len_of_matched_b_without_mod;
	my $matched_int_with_mod = 0;
	my $matched_mod_type;
	my $matched_mod_type_num;
	my $matched_isotops;
	my $matched_isotops_int;
	my $corr_mod;
	my $median_corr_keep;
	my ($matched_int_sd_mod,$ppm_ave_mod,$ppm_sd_mod,$matched_int_reported,$match_type_report);
	my %last_candidate_trace_mz;
	my $pre_SN;
	
	while($index <= $len){

		my $curr_mz = $ms1scan_mz_list->[$index]->[0];
		my $curr_int = $ms1scan_mz_list->[$index]->[1];
		$pre_SN = $curr_int / (1 + $lowest_sig );
		#print STDERR "--check1 mz: $curr_mz\n";
		#print STDERR join "##", exists $trace_mz{$curr_mz}, exists $last_candidate_trace_mz{$curr_mz},$pre_SN, "\n";
		if (exists $trace_mz{$curr_mz} || exists $last_candidate_trace_mz{$curr_mz} || $pre_SN < 10 ){
			$index++;
			next;
		}
		
		#print STDERR "--check2 mz: $curr_mz\n";
		#print STDERR Dumper "MS1Profile: ",$ms1profile->{$curr_mz};
		if (! defined $ms1profile->{$curr_mz} || ($ms1profile->{$curr_mz}->[1] == 0 && $ms1profile->{$curr_mz}->[3] == 0) ){
			$index++;
			next ;
		}
		
		my $mod_mass = ($curr_mz - $cal_premz) * $assumed_charge;
		#print STDERR "unimod_site_test\n";
		#print STDERR join ":",($mod_mass,$strip_pep,$ms1_ppm_tolerence,$peptide_prev_aa,$peptide_next_aa),"\n";
		
		my $site = unimod_site_test($mod_mass,$assumed_charge,$cal_premz,$strip_pep,$ms1_ppm_tolerence,$peptide_prev_aa,$peptide_next_aa);
		
		#print STDERR "--unimod_site_test: $curr_mz\n";
		#print STDERR Dumper $site;
		if( ! defined $site){
			$index++;
			next;
		}

		my $non_zero_in_MS1;
		
		my @ms1arr = @{$ms1profile->{$curr_mz}};
		
		my @obs_iso_distribution;
		push @obs_iso_distribution,$curr_int;
		my $index_ = $index;
		my @predict_ios_peaks;
		my $upper_mz;
		for(1..$iso_require){
			$upper_mz = ($curr_mz + $_ * $massdiff_C12_C13 / $assumed_charge);
			push @predict_ios_peaks, $upper_mz;
		}

		my @predict_ios_peaks_high;
		if($assumed_charge){
			my $iso_i = 1;
			my $c = 2 * $assumed_charge;
			while (1){
				my $i = $curr_mz + $iso_i * $massdiff_C12_C13 / $c;
				last if $upper_mz < $i;
				push @predict_ios_peaks_high, $i;
				$iso_i++;
			}
		}
		
		my $n_of_ios_matched;
		my %tmp_trace_mz;
		
		my $index_predict_ion_peaks = 0;
		my $isotopic_ints;
		my $iso_top_intensity = $curr_int;
		my $control = 0;
		my $last_int;
		my $second_last_int;
		my $iso_top_mz = $curr_mz;
		
		
		my $n_of_ios_matched_high;
		my $index_predict_ion_peaks_high = 0;
		my $isotopic_ints_high;
		my $iso_top_intensity_high = $curr_int;
		my $control_high = 0;
		my %tmp_trace_mz_high;
		
		my @obs_mz;
		while($index_++ < $len){
			my $curr_mz_ = $ms1scan_mz_list->[$index_]->[0];
			my $curr_int_ = $ms1scan_mz_list->[$index_]->[1];
			
			last if $curr_mz_ > $predict_ios_peaks[$#predict_ios_peaks] + 0.1;
			if($curr_int_ < 3 * $lowest_sig ){
				next;
			}

			my $pred_iso = $predict_ios_peaks[$index_predict_ion_peaks];

			my $measured_ppm_iso = ppm($curr_mz_,$pred_iso);

			if(abs( $measured_ppm_iso ) <= $ms1_ppm_tolerence ){
				
				if($second_last_int > $last_int && $last_int < $curr_int_){
					$index_--;
					last;
				}
				push @obs_mz,[$curr_int_,$curr_mz_];
				if($iso_top_intensity < $curr_int_){
					$iso_top_intensity = $curr_int_;
					$iso_top_mz = $curr_mz_;
				}
				
				push @obs_iso_distribution,$curr_int_;
				$n_of_ios_matched++;
				$index_predict_ion_peaks++;
				$isotopic_ints += $curr_int_;
				$tmp_trace_mz{$curr_mz_} = 1;
				$second_last_int = $last_int;
				$last_int = $curr_int_;
			}
			
			my $pred_iso_high = $predict_ios_peaks_high[$index_predict_ion_peaks_high];
			my $measured_ppm_iso_high = ppm($curr_mz_,$pred_iso_high);
			if(abs( $measured_ppm_iso_high ) <= $ms1_ppm_tolerence ){
				$n_of_ios_matched_high++;
				$index_predict_ion_peaks_high++;
				$isotopic_ints_high += $curr_int_;
				$tmp_trace_mz_high{$curr_mz_} = 1;
			}
			
			last if $index_predict_ion_peaks == $iso_require;
				
		}
		# print STDERR "----\t check mz:$curr_mz: mod_mass: $mod_mass n_of_ios_matched_high:$n_of_ios_matched_high\tn_of_ios_matched:$n_of_ios_matched\n";
		# print STDERR Dumper \%tmp_trace_mz,\%tmp_trace_mz_high;
		# print STDERR "--check4 mz: $curr_mz\n";
		if($n_of_ios_matched < 3){
			$index++;
			next;
		}

		my $overlap;
		foreach my $i (keys %tmp_trace_mz_high){
			if(exists $tmp_trace_mz{$i}){
				$overlap++;
			}
		}
		# print STDERR Dumper \%tmp_trace_mz,\%tmp_trace_mz_high;
		# print STDERR "--check5 mz: $curr_mz\n";
		if(  $n_of_ios_matched >= 3 && $overlap == $n_of_ios_matched && $n_of_ios_matched < $n_of_ios_matched_high && $isotopic_ints_high > $isotopic_ints){
			$index++;
			foreach my $i (keys %tmp_trace_mz_high){
				$trace_mz{$i} = 1;
			}
			next;
		}
		
		my $ind_corr_arr = $#obs_iso_distribution > 4 ? 4 : $#obs_iso_distribution;
		my $corr = pearson_corr([@theoretical_iso_distribution[0..$ind_corr_arr]],[@obs_iso_distribution[0..$ind_corr_arr]]);
		
		# print STDERR Dumper (@obs_mz,[@theoretical_iso_distribution[0..$ind_corr_arr]],[@obs_iso_distribution[0..$ind_corr_arr]]);
		# print STDERR "corr check : $corr\n";
		if($corr < 0.8){
			$index++;
			next;
		}
		  
		my ($matched_mzs_type,$matched_int,$matched_ppms,$match_mzs_topN5);
		my ($match_type,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$ppm_sd_mod_,$matched_int_reported_,$matched_coverage_,$match_mzs_topN5,$matched_mod_type_,$mod_matched_coverage_);
		my ($PTMsiteScore_modification_Res) =  PTMsiteScore_modification($pep,$assumed_charge,$site,$deisotope_peaks,$strip_pep,$cal_premz,$ms1_ppm_tolerence);
		
		
		if(! defined $PTMsiteScore_modification_Res){
			# print STDERR "no PTMSiteScore\n";
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}
			$index++;
			next;
		}
		
		my ($PTMSiteScore_,$PTMSiteProb_,$position_,$unimod_mass_,$NL_,$mod_name_,$mod_AA) = @{$PTMsiteScore_modification_Res};
		
		my $theoretical_mz_table_mod = peptide::fragmentation_mod($pep,$assumed_charge,$unimod_mass_,$position_,$mod_name_,$NL_,$pepLen);

		($matched_mzs_type,$matched_int,$matched_ppms,$match_mzs_topN5) = ms2spectra_match_mod($theoretical_mz_table_mod,$peaks,$ms1profile->{$iso_top_mz},$ms2profile);
		
		($match_type,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$matched_int_reported_,$matched_coverage_,$match_mzs_topN5,$matched_mod_type_,$mod_matched_coverage_,$ppm_sd_mod_)= tag_mod($matched_mzs_type,$matched_int,$matched_ppms);
		# print STDERR "--check5 mz: $curr_mz\n";
		my ($y,$b) =  ($match_type->{"Y"}, $match_type->{"B"}); 
		# print STDERR Dumper "match_mzs_topN5",$match_mzs_topN5;
		
		my @corrlist;
		my $ave_corr;
		my $median_corr;
		my ($s,$e) = boundary_index(normalization_int($ms1profile->{$iso_top_mz})); 
		foreach my $ion_mz (@{$match_mzs_topN5}){
			next if ! defined $ms2profile->{$ion_mz};
			my $n_non_zero;
			map{if($_){$n_non_zero++}}@{$ms2profile->{$ion_mz}};
			if ($n_non_zero <= 2 || !( $ms2profile->{$ion_mz}->[1] && $ms2profile->{$ion_mz}->[2] && $ms2profile->{$ion_mz}->[3])){
				next;
			}
			# print STDERR Dumper (normalization_int($ms2profile->{$ion_mz}),normalization_int($ms1profile->{$iso_top_mz}),$s,$e);
			push @corrlist, pearson_corr([@{$ms2profile->{$ion_mz}}[$s..$e]],[@{$ms1profile->{$iso_top_mz}}[$s..$e]]);

		}

		$median_corr = median([@corrlist]);
		# print STDERR "median_corr: $median_corr; median_corr_keep: $median_corr_keep  @corrlist\n";
		# print STDERR "--check7 mz: $curr_mz\n";
		if($median_corr < 0.9  || scalar @corrlist < 4 ){
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}
			$index++;
			next;
		}
	
		# print STDERR "$PTMSiteScore_ > $PTMSiteScore   
			# $matched_coverage_ >= $match_coverage_without_mod + 2  
			# $y >= $len_of_matched_y_without_mod + 2 
			# $b >= $len_of_matched_b_without_mod + 2
			# ($median_corr > $median_corr_keep && $isotopic_ints > $matched_isotops_int)
			# ";
		# print STDERR "--check8 mz: $curr_mz\n";
		%last_candidate_trace_mz = %tmp_trace_mz;
		if(
			$PTMSiteScore_ > $PTMSiteScore && 
			(($y >= $len_of_matched_y_without_mod || $b >= $len_of_matched_b_without_mod) && $matched_coverage_ >= $match_coverage_without_mod )&&
			($median_corr > $median_corr_keep && $isotopic_ints > $matched_isotops_int)
		){
			($PTMSiteScore,$PTMSiteProb,$position,$unimod_mass,$NL,$mod_name) = ($PTMSiteScore_,$PTMSiteProb_,$position_,$unimod_mass_,$NL_,$mod_name_."[".$mod_AA."]");
			$match_coverage_with_mod = $matched_coverage_;
			$mod_matched_coverage = $mod_matched_coverage_;
			$len_of_matched_y_with_mod = $y;
			$len_of_matched_b_with_mod = $b;
			$matched_int_with_mod = $matched_int_with_mod_;
			$matched_isotops_int = $isotopic_ints;
			$matched_mod_type = $matched_mod_type_;
			$matched_mod_type_num = $#$matched_mod_type + 1;
			$median_corr_keep = $median_corr;
			($matched_int_sd_mod,$ppm_ave_mod,$matched_int_reported,$match_type_report,$ppm_sd_mod) = ($matched_int_sd_mod_,$ppm_ave_mod_,$matched_int_reported_,$match_type,$ppm_sd_mod_);
			
			# print STDERR "keep: $curr_mz iso_top_mz: $iso_top_mz   $y,$b,$matched_int_with_mod_ $mod_mass $mod_name  $len_of_matched_y_with_mod, $len_of_matched_b_with_mod, $matched_int_with_mod\n";
			
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}

			$mod_Iso_tag = $n_of_ios_matched +1;
			$mod_exp_premz = $curr_mz;
			$mod_exp_int = $curr_int;
			$mod_ppm_obs = 0;
			$massdiff = $mod_mass;
			$index = $index_;
		}else{
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}
		}
		$index++;
	} 

	my $mod_peptide;
	if($mod_Iso_tag >= 2 ){
		my $mod_AA;
		if($position eq "N-" or $position eq "ProtN-"){
			$mod_peptide = "n[".$unimod_mass."]".$pep;
		}elsif($position eq "C-" or $position eq "ProtC-"){
			$mod_peptide = $pep."[".$unimod_mass."]c";
		}else{
			
			my ($n,$t,$a);
			while($pep =~ m{([A-Z])(\[([^][]+)\])?}g){
				$a = $1;
				$n += length($2);
				$t = $3;
				my $pos_ = pos($pep)-$n;
				#print join "-",$a,$n,$t,$mod_peptide,"\n";
				if($pos_ == $position){
					if($a !~ /\[/){
						$mod_peptide .= $a."[".($mono_aa{$a} + $unimod_mass)."]";
					}else{
						$mod_peptide .= $a."[".($t + $unimod_mass)."]";
					}
				}else{
					$mod_peptide .= $&;
				}
			}
		}
		# print STDERR join "\t","@@@",$prems1scan,$cal_premz,$assumed_charge,$ms2scan,$pep,$pepLen,$strip_pep,$mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$mod_pos,$position,$mod_name,$unimod_mass,$mod_peptide,$PTMSiteScore,$PTMSiteProb,"\n";
		
		my $new_calmz = peptide::calmz($mod_peptide,$assumed_charge);
		return ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$match_type_report->{'Y'},$match_type_report->{'B'},$matched_int_sd_mod,$ppm_ave_mod,$matched_int_reported,$match_coverage_with_mod,$matched_mod_type,$matched_int_with_mod,$pre_SN,$position,$mod_name,$unimod_mass,$mod_peptide,$PTMSiteScore,$PTMSiteProb,$new_calmz,$ppm_sd_mod);	
	}else{
		# print STDERR "\n---FAILED\n";
		return (0,$cal_premz);
	}
}

sub iso_peak_check{
	#print Dumper \@_;
	my ($prems1scan,$cal_premz,$assumed_charge,$iso_require,$ms2scan,$peptide) = @_;
	
	my ($ms1scan_mz_list,$lowest_sig,$highest_sig);
	if(exists $ms1scan_clean_data{$prems1scan}{$ms2scan}){
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = @{$ms1scan_clean_data{$prems1scan}{$ms2scan}};
	}else{
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = MS1file($prems1scan,$ms2scan);
		$ms1scan_clean_data{$prems1scan}{$ms2scan} =[$ms1scan_mz_list,$lowest_sig,$highest_sig];
	}
	my $len = $#$ms1scan_mz_list; #print Dumper $ms1scan_mz_list;exit;
	my $obs_mono_mz;
	my $obs_mono_int;
	my $tag_iso;
	my $measured_ppm_mono;
	my $measured_pre_SN;
	my $isotopic_highest_peak;
	my $iso_require = 10;
	
	my $index = 0;
	if(! defined $ms1scan_mz_list->[$index]->[0]){
		return (0);
		next;
	}
	
	
	while($index <= $len){
		my $curr_mz = $ms1scan_mz_list->[$index]->[0];
		my $curr_int = $ms1scan_mz_list->[$index]->[1];
		if($curr_mz > $cal_premz + 0.1){
			last;
		}
		if($curr_mz <= $cal_premz - 0.1){
			$index++;
			next;
		}else{
			my $measured_ppm = ppm($curr_mz,$cal_premz);
			$measured_pre_SN = $curr_int / (1 + $lowest_sig);
			# if($measured_pre_SN < 3){
				# $index++;
				# next;
			# }
			my @obs_iso_distribution;
			my @obs_ppm;
			my @obs_mz;
			if(abs($measured_ppm) <= $ms1_ppm_tolerence ){
				$measured_ppm_mono = $measured_ppm;
				$obs_mono_mz = $curr_mz;
				$obs_mono_int = $curr_int;
				push @obs_iso_distribution,$obs_mono_int;
				push @obs_ppm,$measured_ppm;
				push @obs_mz,$obs_mono_mz;
				$tag_iso = $tag_iso <= 1 ? 1 : $tag_iso;
				$isotopic_highest_peak = $obs_mono_int;
				if($iso_require == 1){
					$index = $len+1;
				}else{

					my $index_ = $index;
					my @predict_ios_peaks;
					my $upper_mz;
					for(1..$iso_require-1){
						$upper_mz = ($curr_mz + $_ * $massdiff_C12_C13 / $assumed_charge);
						push @predict_ios_peaks, $upper_mz;
					}
					my $n_of_ios_matched;
					my %tmp_trace_mz;
					my $index_predict_ion_peaks = 0;
					my $isotopic_ints;
					my $control = 0;
					my $last_int;
					my $second_last_int;
					
					
					my $n_of_ios_matched_high;
					my $index_predict_ion_peaks_high = 0;
					my $isotopic_ints_high;
					my %tmp_trace_mz_high;

					while($index_++ < $len){
						my $curr_mz_ = $ms1scan_mz_list->[$index_]->[0];
						my $curr_int_ = $ms1scan_mz_list->[$index_]->[1];
						
						last if $curr_mz_ > $predict_ios_peaks[$#predict_ios_peaks] + 0.1;
						my $pred_iso = $predict_ios_peaks[$index_predict_ion_peaks];
						my $measured_ppm_iso = ppm($curr_mz_,$pred_iso);
						if(abs( $measured_ppm_iso ) <= $ms1_ppm_tolerence ){
							if($second_last_int > $last_int && $last_int < $curr_int_){
								$index_--;
								last;
							}
							push @obs_ppm,$measured_ppm_iso;
							push @obs_mz,$curr_mz_;
							$n_of_ios_matched++;
							$index_predict_ion_peaks++;
							push @obs_iso_distribution,$curr_int_;
							$tmp_trace_mz{$curr_mz_} = 1;
							$second_last_int = $last_int;
							$last_int = $curr_int_;
						}

						last if $index_predict_ion_peaks == $iso_require-1;
					}
							
					$n_of_ios_matched = $n_of_ios_matched > $iso_require ? $iso_require : $n_of_ios_matched; # in case two more peak within 0.01Da, use '>=', rather '=='
					$tag_iso = $tag_iso < $n_of_ios_matched +1 ? $n_of_ios_matched + 1 : $tag_iso;
					
				}
			}
			$index++;
		}
	}
	#print join "\t",($tag_iso, $obs_mono_mz, $obs_mono_int,$measured_ppm_mono,$measured_pre_SN,"\n");#exit;
	return ($tag_iso, $obs_mono_mz, $obs_mono_int,$measured_ppm_mono,$measured_pre_SN) ;
}


sub unimod_site_test{
	my $mass = shift;
	my $assumed_charge = shift;
	my $cal_premz = shift;
	my $peptide = shift;
	my $ppm_tol = shift;
	my ($peptide_prev_aa,$peptide_next_aa) = @_;
	my $mod_name;
	my $term_tag;
=head site
 site="A"
 site="C"
 site="C-term"
 site="D"
 site="E"
 site="F"
 site="G"
 site="H"
 site="I"
 site="K"
 site="L"
 site="M"
 site="N"
 site="N-term"
 site="P"
 site="Q"
 site="R"
 site="S"
 site="T"
 site="U"
 site="V"
 site="W"
 site="Y"

 position="Any C-term"
 position="Any N-term"
 position="Anywhere"
 position="Protein C-term"
 position="Protein N-term"

'0.984016' => {
                          'R' => {
                                   'Anywhere' => [
                                                   'Deamidated'
                                                 ]
                                 },
                          'Q' => {
                                   'Anywhere' => [
                                                   'Deamidated',
                                                   'Gln->Glu'
                                                 ]
                                 },
                          'N' => {
                                   'Anywhere' => [
                                                   'Deamidated',
                                                   'Asn->Asp'
                                                 ]
                                 }
                        }

=cut	
	my @unimod_mass_list = sort{$a <=> $b}@{$unimod_info->[1]};
	my $unimods = $unimod_info->[0];
	foreach my $unimod_mono_mass (@unimod_mass_list){
		if($unimod_mono_mass < $mass - 0.1){
			next;
		}elsif($unimod_mono_mass > $mass + 0.1){
			last;
		}else{
			my $ppm_ = ppm($cal_premz+$unimod_mono_mass/$assumed_charge,$cal_premz+$mass/$assumed_charge);
			# print STDERR join "\t","unimod_site_test:",$unimod_mono_mass,$mass,$ppm_,$ppm_tol,"\n";
			my $fix_tag = 0;
			foreach my $fixed_mod_aa (keys %fixed_variable){
				if(abs($fixed_variable{$fixed_mod_aa} + $unimod_mono_mass ) < 0.01){
					$fix_tag = 1;
					while($peptide =~ /$fixed_mod_aa/g){
						my $pos_fixed_aa = pos($peptide);
						push @{$mod_name},[0,$pos_fixed_aa,$fixed_mod_aa,"no-fixed_mod",undef];
					}
				}
			}
			if($fix_tag){next;}
			
			if(abs($ppm_) < $ppm_tol){
				foreach my $site_ (keys %{$unimods->{$unimod_mono_mass}}){
					foreach my $specific (keys %{$unimods->{$unimod_mono_mass}->{$site_}}){
						foreach my $item (@{$unimods->{$unimod_mono_mass}->{$site_}->{$specific}}){
							my $mod_name_ = $item->[1];
							my $NL = $item->[0];
							if($site_ eq "N-term"){
								if($specific eq "Any N-term"){
									push @{$mod_name},[$unimod_mono_mass,"N-",$site_,$mod_name_,$NL];
								}
								if($specific eq "Protein N-term" and $peptide_prev_aa eq "-"){
									push @{$mod_name},[$unimod_mono_mass,"ProtN-",$site_,"ProtN-".$mod_name_,$NL];
									$term_tag->{'ProtN-'} = 1;
								}
							}elsif($site_ eq "C-term"){
								if($specific eq "Any C-term"){
									push @{$mod_name},[$unimod_mono_mass,"C-",$site_,$mod_name_,$NL];
								}
								if($specific eq "Protein C-term" and $peptide_next_aa eq "-"){
									push @{$mod_name},[$unimod_mono_mass,"ProtC-",$site_,"ProtC-".$mod_name_,$NL];
									$term_tag->{'ProtC-'} = 1;
								}
							}else{
								if($specific eq "Anywhere"){
									while($peptide =~ m{$site_}g){
										my $position = pos($peptide);
										push @{$mod_name},[$unimod_mono_mass,$position,$site_,$mod_name_,$NL];
									}
								}elsif($specific eq "Any N-term"){
									if($peptide =~ m{^$site_}){
										push @{$mod_name},[$unimod_mono_mass,"N-",$site_,$mod_name_,$NL];
									}
								}elsif($specific eq "Protein N-term" and $peptide_prev_aa eq "-"){
									push @{$mod_name},[$unimod_mono_mass,"ProtN-",$site_,$mod_name_,$NL];
								}elsif($specific eq "Protein C-term" and $peptide_next_aa eq "-"){
									push @{$mod_name},[$unimod_mono_mass,"ProtC-",$site_,$mod_name_,$NL];
								}elsif($specific eq "Any C-term"){
									if($peptide =~ m{${site_}$}){
										push @{$mod_name},[$unimod_mono_mass,"C-",$site_,$mod_name_,$NL];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	# slim $mod_name, specific for N-term, ProN-term, C-term, ProC-term
	
		# [
            # '42.04695',
            # 20,
            # 'G',
            # 'Gly->Val',
            # undef
          # ],
		  # [
            # '42.04695',
            # 'N-',
            # 'N-term',
            # 'Propyl',
            # undef
          # ],
		  # [
            # '42.010565',
            # 'ProN-',
            # 'N-term',
            # 'Acetyl',
            # undef
          # ],
		   # [
            # '42.04695',
            # 'N-',
            # 'N-term',
            # 'Propyl',
            # undef
          # ],


	if($term_tag->{'ProtC-'} || $term_tag->{'ProtN-'} ){
		my $slim_mod_name;
		foreach my $i (@{$mod_name}){
			my ($massdiff,$pos,$short_name,$mod,$NL) = @{$i};
			if($short_name eq "N-term" and $pos ne "ProtN-"){
				next;
			}
			if($short_name eq "C-term" and $pos ne "ProtC-"){
				next;
			}
			push @{$slim_mod_name},$i;
		}
		return $slim_mod_name;
	}else{
		# print STDERR Dumper $mod_name;
		return $mod_name;
	}
	#}
}


sub normalization_int{
	my $data = shift;
	my $max = max(@{$data});
	my @n = map{ $max ? $_ / $max : 0}@{$data};
	return ([@n]);
}

sub MS1file{
	#** 需要input：指针位置，scan起始和终止，scan mz window **
	my ($ms1scan_,$ms2scan ) = (@_);
	#print Dumper "MS1file", \@_;

	my $premz = $ms2toms1{$ms2scan}->[1];

	my ($mz_start,$mz_end) = ($DIA_window{$premz}{lower} - 5 , $DIA_window{$premz}{higher} + 3);
	my $curr_ms1scan_pos = $ms1_file_position{sprintf("%06.f",$ms1scan_)};

	my $cur_scan = $ms1scan_ + 0;

	my ($curr_ms1scan_data,$ms1scan_lowsig_data,$ms1scan_highest_sig) = readMS1($curr_ms1scan_pos,$mz_start,$mz_end,$ms1scan_basepeak{$cur_scan});

	return ($curr_ms1scan_data,$ms1scan_lowsig_data,$ms1scan_highest_sig);
}

sub readMS1{
	my ($pos,$mz_start,$mz_end,$bp) = (@_);
	#print Dumper "readMS1", \@_;
	my @ms1_peaks_;
	my $lowest_sig = 100000000;
	my $highest_sig = 0;
	seek(MS1,$pos,0);
	my $exit;
	while(<MS1>){
		#print $_;
		chomp;
		if(/^\d+/){
			my @line = split / /,$_;
			if ($line[0] >= $mz_start && $line[0] <= $mz_end){
				push @ms1_peaks_,[@line[0,1]];
				if($lowest_sig > $line[1]){
					$lowest_sig = $line[1];
				}
				if($highest_sig < $line[1]){
					$highest_sig = $line[1];
				}
			}
			last if $line[0] > $mz_end;
		}
		if(/^S/){
			last if !$exit++;
		}
	}
	return ([@ms1_peaks_],$lowest_sig, $highest_sig);	
	
}

sub readMS2{
	my ($ms2scan) = (@_);
	my @ms2_peaks_;
	my $pos = $ms2_file_position{$ms2scan};
	
	seek(MS2,$pos,0);
	my $exit;
	while(<MS2>){
		#print $_;
		chomp;
		if(/^\d+/){
			my @line = split / /,$_;
			push @ms2_peaks_,[@line];
		}
		if(/END/){
			last if !$exit++;
		}
	}
	return ([@ms2_peaks_]);
}

sub readpseudoMS2{
	my ($title) = (@_);
	my @ms2_peaks_;
	my $pos = $pseudomgffilehandle_position{$title};
	my $fh = $pseudomgffilehandles{$title_to_mgffile{$title}};
	seek($fh,$pos,0);
	my $exit;
	while(<$fh>){
		#print $_;
		chomp;
		if(/^\d+/){
			my @line = split / /,$_;
			push @ms2_peaks_,[@line];
		}
		if(/END/){
			last if !$exit++;
		}
	}
	return ([@ms2_peaks_]);
}


sub Spectrum_cmp{
	my ($sp1,$sp2,$ref) = (@_); # sp1 is the control spectrum, sp2 is the test spectrum
	
	my ($sp1_len,$sp2_len) = ($#$sp1,$#$sp2);
	my ($sp1_index,$sp2_index) = (0,0);
	while(1){
		if($sp1_index > $sp1_len){
			last;
		}elsif($sp2_index > $sp2_len){
			while(1){
				my ($mz1,$int1) = (@{$sp1->[$sp1_index]});
				push @{$ref->{$mz1}},0;
				$sp1_index++;
				last if $sp1_index > $sp1_len;
			}
			last;
		}else{
			my ($mz1,$int1) = (@{$sp1->[$sp1_index]});
			my ($mz2,$int2) = (@{$sp2->[$sp2_index]}); 
			if($mz1 && abs(1000000*($mz1-$mz2)/$mz1) <= $ms1_ppm_tolerence){
				push @{$ref->{$mz1}},$int2;
				$sp1_index++;
				$sp2_index++;
			}else{
				if($mz1 > $mz2){
					$sp2_index++;
				}else{
					push @{$ref->{$mz1}},0;
					$sp1_index++;
				}
			}
		}
	}
	return $ref;
}



sub ms2spectra_match_mod{
	my ($pred,$obs,$ms1profile_i,$ms2profile_i) = (@_); # predicted, observed mz
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	my $matched_mzs_type;
	my $matched_mzs;
	my $match_mzs_topN5;
	my $ppms_frag;
	my $matched_int;
	my @matched_ints;
	my $matched_type;
	my %checked_obs_mz;
	
	my ($s,$e) = boundary_index(normalization_int($ms1profile_i)); 
		
	while($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
		
		if(defined $checked_obs_mz{$obs_index2}){
			$obs_index2++;
		}
		if($pred_index1 > $pred_len1){
			while($obs_index2 <= $obs_len2){
				$obs_index2++;
			}
			last;
		}
		if($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
			my ($pred_mz1,$obs_mz2,$obs_int2) = ($pred_list[$pred_index1],@{$obs->[$obs_index2]});
			my $measured_ppm = ppm($obs_mz2,$pred_mz1);
			if(abs($measured_ppm) <= $ms2_ppm_tolerence ){
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				my $corr_i = pearson_corr([@{$ms2profile_i->{$obs_mz2}}[$s..$e]],[@{$ms1profile_i}[$s..$e]]);
				if($corr_i > 0.8){
					push @{$matched_mzs},[$obs_mz2,$obs_int2];
					push @{$matched_type},$predicted_mz_ion_type;
				
					push @{$matched_mzs_type},$predicted_mz_ion_type;
					$ppms_frag->{$predicted_mz_ion_type} = $measured_ppm;
					$matched_int->{$predicted_mz_ion_type} = [$obs_mz2,$obs_int2];
				}
				$checked_obs_mz{$obs_index2} = 1;
				$pred_index1++;
				$obs_index2++;
			}else{
				if($pred_mz1 < $obs_mz2){
					$pred_index1++;
				}else{
					$obs_index2++;
				}
			}
		}
	}
	
	my @sorted = map{$_->[0]}sort{$b->[1] <=> $a->[1]}@{$matched_mzs};
	$match_mzs_topN5 = [@sorted[0..4]];
	return ($matched_mzs_type,$matched_int,$ppms_frag,$match_mzs_topN5);
}



sub ms2spectra_match{
	my ($pred,$obs) = (@_); # predicted, observed mz
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	my $matched_mzs_type;
	my $matched_mzs;
	my @ppms_frag;
	my $matched_int;
	my @matched_ints;
	my %checked_obs_mz;
	while($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
		
		if(defined $checked_obs_mz{$obs_index2}){
			$obs_index2++;
		}
		if($pred_index1 > $pred_len1){
			while($obs_index2 <= $obs_len2){
				$obs_index2++;
			}
			last;
		}
		if($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
			my ($pred_mz1,$obs_mz2,$obs_int2) = ($pred_list[$pred_index1],@{$obs->[$obs_index2]});
			my $measured_ppm = ppm($obs_mz2,$pred_mz1);
			if(abs($measured_ppm) <= $ms2_ppm_tolerence ){
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				push @{$matched_mzs},$obs_mz2;
				push @{$matched_mzs_type},$predicted_mz_ion_type;# if $predicted_mz_ion_type =~ /[by]\+/;;
				push @ppms_frag,$measured_ppm;# if $predicted_mz_ion_type =~ /[by]\+/;
				$matched_int += $obs_int2;# if $predicted_mz_ion_type =~ /[by]\+/;
				push @matched_ints, log($obs_int2);# if $predicted_mz_ion_type =~ /[by]\+/;
				$checked_obs_mz{$obs_index2} = 1;
				$pred_index1++;
				$obs_index2++;
			}else{
				if($pred_mz1 < $obs_mz2){
					$pred_index1++;
				}else{
					$obs_index2++;
				}
			}
		}
	}
	my $ppm_md; 
	my $ppm_sd;
	
	
	my $m_int_ave;
	my $m_int_sd;
	if(scalar @matched_ints > 1){
		$m_int_ave = ave(\@matched_ints);
		$m_int_sd = stdev(\@matched_ints,$m_int_ave);
		$ppm_md = median(\@ppms_frag);
		$ppm_sd = stdev(\@ppms_frag,ave(\@ppms_frag));
	}
	#print join "\t",($matched_mzs_type,$matched_int,$m_int_sd || "NA",$ppm_ave,$matched_mzs,$ppm_sd),"\n";
	return ($matched_mzs_type,$matched_int,$m_int_sd || "NA",$ppm_md,$ppm_sd);
}


sub tag{
	my $ref = shift;
	my $return_type = shift; # type or coverage
	my %ion_list;
	my %ion_for_coverage;
	my %uniq;
	
		foreach my $ion (@{$ref}){    
			my ($index,$fragment_type) = split /\./,$ion;
			next if $fragment_type =~ /_NH3|_H2O/;
			my $type;
			if($fragment_type =~ /y/){
				$type = "Y";
			}
			if($fragment_type =~ /b/){
				$type = "B";
			}
			push @{$ion_list{$type}},$index if !$uniq{$type}{$index}++;
			$ion_for_coverage{$index} = 1;
		}
		my %tag_len;
		foreach my $type (keys %ion_list){
			my @n = sort{$a <=> $b}@{$ion_list{$type}};
			my $tag=1;
			my $max_len;
			foreach my $i(0..$#n-1){
				if($n[$i]+1 == $n[$i+1]){
					$tag++;
					
				}elsif($n[$i]+2 == $n[$i+1]){
					$tag++;
					$i++;
				}else{
					$max_len = $tag if $max_len < $tag;
					$tag=1;
				}
			}
			$max_len = $tag if $max_len < $tag;
			
			$tag_len{$type} = $max_len; 
		}
		
	my @n = sort{$a <=> $b}(keys %ion_for_coverage);
	my $tag=1;
	my $max_len;
	foreach my $i(0..$#n-1){
		if($n[$i]+1 == $n[$i+1]){
			$tag++;
		}elsif($n[$i]+2 == $n[$i+1]){
			$tag++;
			$i++;
		}else{
			$max_len = $tag if $max_len < $tag;
			$tag=1;
		}
	}
	$max_len = $tag if $max_len < $tag;
	my $coverage_n = $max_len;
		
	return (\%tag_len,$coverage_n);	
}


sub tag_mod{
	my $matched_ions = shift;
	my $matched_ints = shift;
	my $matched_ppms = shift;
	my $matched_ints_mod;
	my $matched_mod_type;
	my @matched_mods;
	my $match_mods_topN5;
	my $matched_ints_mod_reported;
	my @ppms;
	my @ints;
	my %mod_types;
	my %ion_list;
	my %ion_for_coverage;
	my %uniq;
	my ($type,$index,$fragment_type);
	foreach my $ion (@{$matched_ions}){    
		($index,$fragment_type) = split /\./,$ion;
		if($fragment_type =~ /y/){
			$type = "Y";
		}
		if($fragment_type =~ /b/){
			$type = "B";
		}
		
		if($uniq{$index.$type}->[0] < $matched_ints->{$ion}->[1]){
			$uniq{$index.$type} = [$matched_ints->{$ion}->[1],$index,$index.".".$fragment_type,$matched_ints->{$ion}->[0]];
		}	
	}
	foreach my $i (keys %uniq){
		$fragment_type = $uniq{$i}->[2];
		if($fragment_type =~ /y/){
			$type = "Y";
		}
		if($fragment_type =~ /b/){
			$type = "B";
		}
		if ($fragment_type =~ /mod/){
			$matched_ints_mod += $uniq{$i}->[0];
			push @{$matched_mod_type},$fragment_type;
			$mod_types{$uniq{$i}->[1]} = $fragment_type ;
		}
		push @{$ion_list{$type}},$uniq{$i}->[1];
		push @matched_mods,[$uniq{$i}->[3],$uniq{$i}->[0],$uniq{$i}->[2]];
		$ion_for_coverage{$uniq{$i}->[1]} = 1;
		$matched_ints_mod_reported += $uniq{$i}->[0];
		push @ints,$uniq{$i}->[0];
		push @ppms,$matched_ppms->{$uniq{$i}->[2]};
		
	}
	my %tag_len;
	foreach my $type (keys %ion_list){
		my @n = sort{$a <=> $b}@{$ion_list{$type}};
		my $tag=1;
		my $max_len;
		foreach my $i(0..$#n-1){
			if($n[$i]+1 == $n[$i+1]){
				$tag++;
			
			}elsif($n[$i]+2 == $n[$i+1]){
				$tag++;
				$i++;
			}else{
				$max_len += $tag if $tag >= 3;
				$tag=1;
			}
		}
		$max_len = $tag if $max_len < $tag;
		$tag_len{$type} = $max_len; 
	}
	my @n = sort{$a <=> $b}(keys %mod_types);
	my $tag=1;
	my $max_len;
	foreach my $i(0..$#n-1){
		if($n[$i]+1 == $n[$i+1]){
			$tag++;

		}elsif($n[$i]+2 == $n[$i+1]){
			$tag++;
			$i++;
		}else{
			$max_len += $tag if $tag >= 3;
			$tag=1;
		}
	}
	$max_len = $tag if $max_len < $tag;
	my $mod_coverage_n = $max_len;
	
	
	my @n = sort{$a <=> $b}(keys %ion_for_coverage);
	my $tag=1;
	my $max_len;
	foreach my $i(0..$#n-1){
		if($n[$i]+1 == $n[$i+1]){
			$tag++;
		}elsif($n[$i]+2 == $n[$i+1]){
			$tag++;
			$i++;
		}else{
			$max_len += $tag if $tag >= 3;
			$tag=1;
		}
	}
	$max_len = $tag if $max_len < $tag;
	my $coverage_n = $max_len;
	
	my $ppm_md = median(\@ppms);
	my $ppm_sd;
	if(scalar @ppms){
		$ppm_sd = stdev(\@ppms, ave(\@ppms));
	}
	my $m_int_ave;
	my $m_int_sd;
	if(scalar @ints > 1){
		$m_int_ave = ave(\@ints);
		$m_int_sd = stdev(\@ints,$m_int_ave);
	}
	
	my @sorted = map{$_->[0]}sort{$b->[1] <=> $a->[1]}@matched_mods;
	$match_mods_topN5 = [@sorted[0..4]];
	# print STDERR "match_mods_type: ",join "\t",@{$matched_mod_type},"\n" if defined $matched_mod_type;
	# print STDERR "match_topN5_type: ",join "\t",(map{if($_){$_->[2]}}(sort{$b->[1] <=> $a->[1]}@matched_mods)[0..4]),"\n" if scalar @matched_mods;
	# print STDERR "mod_coverage_n: $mod_coverage_n\n";
	return (\%tag_len,$matched_ints_mod,$m_int_sd || "NA",$ppm_md,$matched_ints_mod_reported,$coverage_n,$match_mods_topN5,$matched_mod_type,$mod_coverage_n,$ppm_sd);
}


sub pearson_corr{
    my ($ref_a, $ref_b) = @_;

    my @x = @{$ref_a};
    my @y = @{$ref_b};
	my $mean_x = sum(@x)/($#x+1);
	my $mean_y = sum(@y)/($#y+1);
	my $correlation;
		 my $N = $#x;
		 my $sum_sq_x = 0;
		 my $sum_sq_y = 0;
		 my $sum_coproduct = 0;
		 for(my $i=0;$i<=$N;$i++){
				  my $delta_x = $x[$i] - $mean_x;
				  my $delta_y = $y[$i] - $mean_y;
				  $sum_sq_x += $delta_x * $delta_x ;
				  $sum_sq_y += $delta_y * $delta_y ;
				  $sum_coproduct += $delta_x * $delta_y;
		 }
		 my $pop_sd_x = sqrt( $sum_sq_x );
		 my $pop_sd_y = sqrt( $sum_sq_y );
		 my $cov_x_y = $sum_coproduct;
		 $correlation = $cov_x_y / ($pop_sd_x * $pop_sd_y) if $pop_sd_x * $pop_sd_y;
	return $correlation;
}

sub deisotope{
	my $peaks = shift;
	
	my @deisotope_peaks;
	my $peakslen = $#$peaks;
	my $ind = 0;
	my %rm_inds;
	my ($mz_last,$int_last) = @{$peaks->[$ind]};
	push @deisotope_peaks, $peaks->[$ind];
	my ($mz_next,$ind_next);
	my @tmp_ind;
	while(++$ind <= $peakslen){
		if($rm_inds{$ind} == 1){
			next;
		}
		($mz_next,$ind_next) = @{$peaks->[$ind]};
		($mz_last,$int_last) = ($mz_next,$ind_next);
		my $ppm_next;

			# charge 2
			my $ind_ = $ind;
			while(++$ind_ < $peakslen){
				if($rm_inds{$ind_} == 1){
					next;
				}			
				($mz_next,$ind_next) = @{$peaks->[$ind_]};
				if($mz_next > $mz_last + 0.5+0.1){
					last;
				}
				my $mz_z2 = $mz_last + $massdiff_C12_C13/2;
				$ppm_next = ppm($mz_next,$mz_z2);

				if(abs($ppm_next) <= $ms2_ppm_tolerence && $ind_next < $int_last * 1.5){
					($mz_last,$int_last) = ($mz_next,$ind_next);
					$rm_inds{$ind_} = 1;
				}
			}
				
			# charge 1
			my $ind_ = $ind;
			while(++$ind_ < $peakslen){
				if($rm_inds{$ind_} == 1){
					next;
				}
				($mz_next,$ind_next) = @{$peaks->[$ind_]};
				if($mz_next > $mz_last + 1.1){				
					last;
				}
				my $mz_z1 = $mz_last + $massdiff_C12_C13;
				$ppm_next = ppm($mz_next,$mz_z1);
				if(abs($ppm_next) <= $ms2_ppm_tolerence && $ind_next < $int_last * 1.5){
					($mz_last,$int_last) = ($mz_next,$ind_next);
					$rm_inds{$ind_} = 1;
				}
			}

		push @deisotope_peaks, $peaks->[$ind];

	}
	return [@deisotope_peaks];
}


sub peakpick{
	my $peaks = shift;
	my %window;
	my %sorted_peaks;
	my %extracted_ions_num;
	foreach my $i (@{$peaks}){
		my ($mz,$int) = @{$i};
		my $range = ceil($mz/100);
		push @{$window{$range}},$i;
	}
	foreach my $w (grep{$_}keys %window){
		my @sorted_int = sort{$b->[1] <=> $a->[1]}grep{$_->[0]}@{$window{$w}};
		foreach my $i (0..$max_peakdepth -1 ){
			if($sorted_int[$i]->[0]){
				$sorted_peaks{$sorted_int[$i]->[0]} = [$i+1,$w];
				map{$extracted_ions_num{$_}++}$i+1..$max_peakdepth;
			}
		}
	}
	return (\%sorted_peaks,\%extracted_ions_num);
}


sub ms2spectra_Site{
	my ($pred,$obs,$top8_mz) = (@_); # predicted, observed mz
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	my $matched_mzs_type;
	my $matched_mzs;
	my @ppms_frag;
	my $matched_int;
	my @matched_ints;
	my $therotical_peak_n;

	my %checked_obs_mz;
	while($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){

		if(defined $checked_obs_mz{$obs_index2}){
			$obs_index2++;
		}
		if($pred_index1 > $pred_len1){
			while($obs_index2 <= $obs_len2){
				$obs_index2++;
			}
			last;
		}
		if($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
			my ($pred_mz1,$obs_mz2,$obs_int2) = ($pred_list[$pred_index1],@{$obs->[$obs_index2]});
			$therotical_peak_n->{$pred_mz1} = 1 if $pred_mz1 < $obs->[$#$obs]->[0];
			my $measured_ppm = ppm($obs_mz2,$pred_mz1);
			if(abs($measured_ppm) <= $ms2_ppm_tolerence ){
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				my $top_n = $top8_mz->{$obs_mz2}->[0];
				foreach my $i (1..$max_peakdepth){
					push @{$matched_mzs_type->{$i}},$predicted_mz_ion_type if $i >= $top_n and $top_n;
				}
				$checked_obs_mz{$obs_index2} = 1;
				$pred_index1++;
				$obs_index2++;
			}else{
				if($pred_mz1 < $obs_mz2){
					$pred_index1++;
				}else{

					$obs_index2++;
				}
			}
		}
	}
	my $therotical_n = scalar (keys %{$therotical_peak_n});
	return ($matched_mzs_type,$therotical_n);
}

sub PeakdepthOptimizationScore{

	my $matched_info = shift;
	my $w = shift;
	my $n = shift;
	my $site_pos = shift;
	my $Score = shift;
	my $Npeak = shift;
	my $_score = shift;

	my $P;

	# p = Npeak * d / w;
	# d = ms2_tolerence 0.02;
	# w = full range of ms2
	
	# P = n!/k!(n-k)! * p^k * (1-p)^(n-k)
	# Score = -10 * log(P)
	my $peakdepth;
	for($peakdepth=1;$peakdepth<=$max_peakdepth;$peakdepth++){
		my $s = 0;;
		if(defined $matched_info->{$peakdepth}){
			my $mod_frag_cnt = 0;
			map{if($_ =~ /mod/){
				$mod_frag_cnt++;
				
				}
			}@{$matched_info->{$peakdepth}};
			if($mod_frag_cnt == 0){
				$s = 0;
				$Score->{$peakdepth}->{$site_pos} = $s;
				next;
			}
			my $k = scalar @{$matched_info->{$peakdepth}};
			my $p = $Npeak->{$peakdepth} * 0.02 / $w;
			if($p == 0){
				$s = 0;
			}else{	
				$P = Cumulative_binominal($n,$k,$p);
				if($P == 0){
					$s = 0;
				}else{
					$s = -10 * log($P)/log(10);
				}
			}
		}else{
			$s = 0;
		}
		$Score->{$peakdepth}->{$site_pos} = $s;
		$_score->{$site_pos}+= $s;
	}
	return ($Score,$_score);
}

sub factorial{
	my $x = shift;
	my $i;
	my $f = 1;
	for($i = 2; $i <= $x; $i++){
		$f *= $i;
	}
	return $f;
}

sub comb_p{
	my ($n,$k) =(@_);
	my $comb = 1;
	my $product_k;
	if(exists $factorial_buffer{$k}){
		$product_k = $factorial_buffer{$k};
	}else{
		$product_k = factorial($k);
		$factorial_buffer{$k} = $product_k;
	}
	$comb /= $product_k;
	my $i;
	for($i=$n;$i>$n-$k;$i--){
		$comb *= $i;
	}
	return $comb;
}

sub binominal{
	my ($n, $k, $p) = (@_);
	my $biprob = $n == $k ?  ( $p ** $k ) * ((1 - $p) ** ($n - $k)) : comb_p($n,$k) * ( $p ** $k ) * ((1 - $p) ** ($n - $k));
	if($biprob eq "Inf"){
		$biprob = 0;
	}
	return $biprob;
}
sub Cumulative_binominal(){
	my ($n, $k, $p) = (@_);
	my $cum_prob;
	foreach my $i($k..$n){
		$cum_prob += binominal($n,$i,$p);
	}
	return $cum_prob;
}
sub sitecombination{
	my $sites = shift;
	my $n = shift;
	my @comb = combine($n,@{$sites});
	return [@comb];
}

sub PhosphoSiteScore{
	my ($peptide_mod,$charge,$Site,$deisotope_peaks,$site_num,$strip_pep) = (@_);
	my ($TOP6_peaks,$extracted_ions_n) = peakpick($deisotope_peaks);
	my $score;
	my $score_for_peakdepth;
	my @SitePos;
	my %modfied_sites;
	while($strip_pep =~ m{[$Site]}g){
		my $pos = pos($strip_pep);
		push @SitePos,$pos;
		$modfied_sites{$pos} = $&;
	}
	
	my @CombSitePos = @{sitecombination([@SitePos],$site_num)};
	my %start_points_adj;
	foreach my $c(@CombSitePos){
		my $putative_frags = peptide::fragmentation_phosSite($peptide_mod,$charge,$c,length($strip_pep));
		my ($matched_frags,$therotical_peak_n) = ms2spectra_Site($putative_frags,$deisotope_peaks,$TOP6_peaks);
		($score,$score_for_peakdepth) = PeakdepthOptimizationScore($matched_frags,$deisotope_peaks->[$#$deisotope_peaks]->[0],$therotical_peak_n,(join "#",@$c),$score,$extracted_ions_n,$score_for_peakdepth);
	}
	if(! defined $score){
		return ($peptide_mod,"-","-","-",0);
	}else{
		my ($peakdepth_of_max_score,$pos_of_peakdepth) = max_diff_score($score);
		
		my $modification_sitepos = $pos_of_peakdepth;
		my %PTMSiteScore_res;
		my $sumScore;
		
		foreach my $c(@CombSitePos){
			my $pos = join "#",@{$c};
			# print STDERR Dumper "pos",$pos;
			foreach my $site (@{$c}){
				$PTMSiteScore_res{$site} += 1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;
				# print STDERR "$score->{$peakdepth_of_max_score}->{$pos} > 0;\n"
				
			}
			$sumScore +=  1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;
		}

		my $PTM_peptide = $peptide_mod;
		my $PTM_score_loc;
		my $ptm;
		my $max_prob;
		my $localized_site;
		my $localized_or_not;
		if($sumScore){
			# print STDERR "PTMsiteScore\t",$peptide_mod,":\t", $peakdepth_of_max_score, "\t",$score->{$peakdepth_of_max_score}->{$modification_sitepos},"\t";
			# print STDERR join ";",map{sprintf "%d(%.2f%%)",$_,100*$PTMSiteScore_res{$_}/$sumScore}@SitePos;
			
			my @siteprob;
			map{
				push @siteprob,[$_,100*$PTMSiteScore_res{$_}/$sumScore];
			}@SitePos;
			my @sorted_siteprob = sort{$b->[1] <=> $a->[1]}@siteprob;
			$localized_site = join ";",map{$_->[0] . "(" . $_->[1] . ")"}@sorted_siteprob[0..$site_num-1];
			$ptm = join ";",map{sprintf "%d(%.2f%%)",$_,100*$PTMSiteScore_res{$_}/$sumScore}@SitePos;
			map{$PTM_peptide = peptide_phos_add($PTM_peptide, $_->[0])}@sorted_siteprob[0..$site_num-1];	
			$localized_or_not = $sorted_siteprob[$site_num-1]->[1] >= 75 ? 1 : 0;
			$PTM_score_loc = $score->{$peakdepth_of_max_score}->{$modification_sitepos};
			# print STDERR "\n";
		}
		return ($PTM_peptide,$PTM_score_loc,$ptm,($localized_site ? $localized_site  : "-"),$localized_or_not);
	}
}

sub max_diff_score{
	my $score_ = shift;
	my %score_diff;
	my %diff_depth;
	my %depth_pos;
	my %max_score_depth;
	my $max_score;
	# print STDERR Dumper $score_;
	foreach my $depth(keys %{$score_}){
		my @sorted_pos = (sort{$score_->{$depth}->{$b} <=> $score_->{$depth}->{$a}}keys %{$score_->{$depth}});
		push @{$max_score_depth{$score_->{$depth}->{$sorted_pos[0]}}},$depth;
		if($max_score < $score_->{$depth}->{$sorted_pos[0]}){
			$max_score = $score_->{$depth}->{$sorted_pos[0]};
		}
		my $diff; 
		if((scalar @sorted_pos) == 1){
			$diff = $score_->{$depth}->{$sorted_pos[0]};
		}else{
			$diff = $score_->{$depth}->{$sorted_pos[0]} - $score_->{$depth}->{$sorted_pos[1]};
		
			if($diff == 0){
				$diff = $score_->{$depth}->{$sorted_pos[0]} - $score_->{$depth}->{$sorted_pos[2]};
			}
			if($diff == 0){
				$diff = $score_->{$depth}->{$sorted_pos[0]} - $score_->{$depth}->{$sorted_pos[3]};
			}
		}
		$depth_pos{$depth} = $sorted_pos[0];
		$score_diff{$depth} = $diff;
		push @{$diff_depth{$diff}},$depth;
	}
	# print STDERR Dumper \%score_diff,\%diff_depth,\%depth_pos;
	my $max_diff = max(map{$score_diff{$_}}@{$max_score_depth{$max_score}});
	if(defined $diff_depth{$max_diff}){
		if(scalar @{$diff_depth{$max_diff}} < 10 ){
			my $peakdepth_of_max_score = min(@{$diff_depth{$max_diff}});
			my $pos_with_max_score = $depth_pos{$peakdepth_of_max_score};

			return ($peakdepth_of_max_score,$pos_with_max_score);
		}
	}
# $VAR1 = {
          # '3' => {
                   # '3' => '332.656787280441',
                   # '4' => '426.172892165666',
                   # '6' => '332.656787280441'
                 # },
          # '1' => {
                   # '3' => '130.812210619349',
                   # '6' => '130.812210619349',
                   # '4' => '192.011234096827'
                 # },
          # '4' => {
                   # '4' => '406.667522431782',
                   # '6' => '316.813202271232',
                   # '3' => '346.266201499673'
                 # },
          # '6' => {
                   # '3' => '469.18119275818',
                   # '6' => '409.212079955519',
                   # '4' => '499.893259964674'
                 # },
          # '5' => {
                   # '6' => '391.459850995456',
                   # '4' => '482.850999033465',
                   # '3' => '451.904279789357'
                 # },
          # '2' => {
                   # '3' => '230.618822167801',
                   # '6' => '172.150267469091',
                   # '4' => '291.583582954626'
                 # }
        # };	
	
}

sub peptide_phos_add{
	my ($pep,$site) = (@_);
	my ($n,$t,$a);
	my $mod_peptide;
	my %phos_mass = (
		'S' => 166.9984,
		'T' => 181.0140,
		'Y' => 243.0297
	);
	while($pep =~ m{([A-Z])(\[([^][]+)\])?}g){
		$a = $1;
		$n += length($2);
		$t = $3;
		my $pos_ = pos($pep)-$n;
		if($pos_ == $site){
			$mod_peptide .= $a."[".$phos_mass{$a}."]";
			
		}else{
			$mod_peptide .= $&;
		}
	}
	return $mod_peptide;
}


sub PTMsiteScore_modification{
	my ($peptide_mod,$charge,$Site,$deisotope_peaks,$strip_pep,$cal_premz,$ms1_ppm_tolerence) = (@_);
	my ($TOP6_peaks,$extracted_ions_n) = peakpick($deisotope_peaks);
	
	
	my %mod_aa;
	while($peptide_mod =~ /(.)\[([^][]+)\]/g){
		$mod_aa{$1} = $2;
	}
	my $CombSitePos;
	my $Site_uniq;
	my %modfied_sites;
	
=head Dumper $Site
$Site = [
          [
            '79.956815',
            1,
            'S',
            'Sulfo',
            '79.956815'
          ],
          [
            '79.966331',
            2,
            'H',
            'Phospho',
            undef
          ]
]
=cut	
	
	foreach my $i (@{$Site}){
		my $pos = $i->[1];
		my $unimod_mass = $i->[0];
		my $NL = $i->[-1];
		my $unimod_name = $i->[-2];
		my $mod_AA = $i->[2];
		if(exists $mod_aa{$mod_AA}){
			
			my $new_mod_mass = $mod_aa{$mod_AA} - $mono_aa{$mod_AA} + $unimod_mass;
			my $new_site = unimod_site_test($new_mod_mass,$charge,$cal_premz,$strip_pep,$ms1_ppm_tolerence,$peptide_prev_aa,$peptide_next_aa);
			# print STDERR Dumper $new_site;
			foreach my $new_i (@{$new_site}){
				my $new_pos = $new_i->[1];
				my $new_unimod_mass = $new_i->[0];
				my $new_NL = $new_i->[-1];
				my $new_unimod_name = $new_i->[-2];
				my $new_mod_AA = $new_i->[2];
				if($new_mod_AA ne $mod_AA){next;}
				my $new_uniq_key = join ":",$new_pos,$new_unimod_mass,$new_NL,$new_mod_AA,"re_loc";
				$Site_uniq->{$new_uniq_key}{$new_unimod_name} = 1;
			}
			next;
		}
		
		my $uniq_key = join ":",$pos,$unimod_mass,$NL,$mod_AA,"loc";
		$Site_uniq->{$uniq_key}{$unimod_name} = 1;
		my $aa = $i->[2];
		$modfied_sites{$pos} = $aa;
	}
	# print STDERR Dumper $Site_uniq;
	foreach my $key (keys %{$Site_uniq}){
		my ($pos,$unimod_mass,$NL,$mod_AA,$loc) = split /:/,$key;
		my $unimod_name = join ";",sort{$a cmp $b}(keys %{$Site_uniq->{$key}});
		push @{$CombSitePos->{$unimod_mass}}, [$pos,$unimod_mass,$NL,$unimod_name,$mod_AA,$loc];
	}
	
	my $ptmscore_res;
	foreach my $unimod_mass(keys %{$CombSitePos}){
		my $score;
		my $score_for_peakdepth;
		my %start_points_adj;
		my $pos_cnt  = scalar @{$CombSitePos->{$unimod_mass}};
		foreach my $i (@{$CombSitePos->{$unimod_mass}}){
			
			my ($site_position,$unimod_mono_mass,$NL,$mod_name_,$aa,$loc) = @{$i};
			# print STDERR join "\t",($site_position,$unimod_mono_mass,$NL,$mod_name_,$aa) ,"\n";
			my $new_peptide_mod;
			if($loc eq "re_loc"){
				
				my $old_peptide_mod = $peptide_mod;
				my ($n,$t,$a);
				while($old_peptide_mod =~ m{([A-Z])(\[([^][]+)\])?}g){
					$a = $1;
					$n += length($2);
					$t = $3;
					my $pos_ = pos($old_peptide_mod)-$n;
					if($pos_ == $site_position){
						$new_peptide_mod .= $a;
						
					}else{
						$new_peptide_mod .= $&;
					}
				}
				# print STDERR Dumper $peptide_mod,$new_peptide_mod;
			}else{
				$new_peptide_mod = $peptide_mod;
			}
			my $putative_frags = peptide::fragmentation_mod($new_peptide_mod,$charge,$unimod_mono_mass,$site_position,$mod_name_,$NL,length($strip_pep));
			my ($matched_frags,$therotical_peak_n) = ms2spectra_Site($putative_frags,$deisotope_peaks,$TOP6_peaks);
			($score,$score_for_peakdepth) = PeakdepthOptimizationScore($matched_frags,$deisotope_peaks->[$#$deisotope_peaks]->[0],$therotical_peak_n,$site_position,$score,$extracted_ions_n,$score_for_peakdepth);
		}
		# print STDERR Dumper "score",$score;
		if(! defined $score){
			$ptmscore_res = undef;
		}else{
			my ($peakdepth_of_max_score,$pos_of_peakdepth) = max_diff_score($score);
			my $modification_sitepos;
			my $modification_unimod_mono_mass;
			my $modification_NL;
			my $modification_mod_name;
			my $modification_mod_AA;
			my $max_diff_score;
			
			my %PTMSiteScore_res;
			my $sumScore;
			foreach my $c(@{$CombSitePos->{$unimod_mass}}){
				my $pos = $c->[0];
				$PTMSiteScore_res{$pos} = 1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;

				$sumScore +=  1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;
				if($pos eq $pos_of_peakdepth){
					($modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA) = (@{$c});
				}	
			}
			my $PTM_peptide;
			my $ptmscore;
			my $PTMSiteScore;
			my $new_calmz;
			if($sumScore){
				# print STDERR $peptide_mod,":\t unimod_mass: $unimod_mass\t","peakdepth_of_max_score: ", $peakdepth_of_max_score, "\t", $score->{$peakdepth_of_max_score}->{$modification_sitepos},"\t";
				# print STDERR join ";",map{sprintf "%d(%.2f%%)",$_->[0],100*$PTMSiteScore_res{$_->[0]}/$sumScore}@{$CombSitePos->{$unimod_mass}};
				# print STDERR "\n";
				my $ind;

				$ptmscore = sprintf ("%.2f",$PTMSiteScore_res{$modification_sitepos}/$sumScore);
				if( $ptmscore >= 0.75){
					$PTM_peptide = $peptide_mod;
					$PTMSiteScore = $score->{$peakdepth_of_max_score}->{$modification_sitepos};		
				}
			}
			if($ptmscore_res->[0] < $PTMSiteScore){
				$ptmscore_res = [$PTMSiteScore,$ptmscore,$modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA];
			}elsif($ptmscore_res->[0] == $PTMSiteScore && $ptmscore_res->[2] == $modification_sitepos){
				$ptmscore_res = [$PTMSiteScore,$ptmscore,$modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name.";".$ptmscore_res->[5],$modification_mod_AA];
			}
		}
	}
	return ($ptmscore_res);
}

1;
