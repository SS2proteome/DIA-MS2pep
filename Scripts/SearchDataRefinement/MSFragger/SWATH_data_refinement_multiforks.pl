#!/usr/bin/env perl
#
# pepxml_post_script [$file] [$charge_dir] [${file}.mgf] [ms1 ppm]

use strict;
use Data::Dumper;
use File::Basename;
use lib dirname(__FILE__);
use POSIX qw(floor ceil round);
use List::Util qw(max min sum);
use peptide;

use IsotopeDistribution;
use MassAccuracy;
use Parallel::ForkManager;
use SWATH_pepXML_to_pin_inparallel;

$| = 1;

if($#ARGV < 7){
	print "Usage: script [filename] [dir] [ms1ppm] [ms2ppm] [PTM_or_not] [ptmfind] [fasta] [processes]
	";
	exit;
}

my $filename = $ARGV[0];
my $dir = $ARGV[1];
my $outfilename = "${filename}.pesudo.mgf";
my $ms1_ppm_tolerence = $ARGV[2]; #ppm
my $ms2_ppm_tolerence = $ARGV[3]; #ppm
my $ptm = $ARGV[4];
my $ptmfind = $ARGV[5];
my $fastafile = $ARGV[6];
my $max_processes = $ARGV[7];

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

#my $ppm_tolerence = 15;
#my $iso_low_precison = 25;
#my $iso_high_precison = 0.001;
my $min_n_iso_required = 4;
my $massdiff_C12_C13 = 1.0033548378;
my $max_peakdepth = 10;
my $proton = 1.00727646677;

my $ms1_file = "$filename.ms1";
my $mzML_file = "$filename.mzML";
my $mgf_file = "$filename.mgf";
my $ms1_min_intensity ;


my @pseudomgffiles = glob "$dir/w_[0-9]*_".$filename.".mgf"; #glob "tmp_[0-9]*_$file.mgf";
#my $outfilename = $ARGV[2];
#my @pepxmlfiles = glob  "$dir/[z][1234]*sub*".$filename.".pepXML";
my @pepxmlfiles = glob  "$dir/[z][12345].w_[0-9]*_".$filename.".pepXML";
#my @pepxmlfiles = "$dir/test.xml";
@pepxmlfiles = glob  "$dir/[z][1234]*".$filename.".pepXML" if ! (scalar @pepxmlfiles);
my $window_file = "$filename.mzML.SWATH_acquisition_window.txt";
print Dumper \@pseudomgffiles,\@pepxmlfiles;

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

my %SWATH_window;
my %isolationlist;
my $num_of_windows;
my %SWATH_window_count;
open(FH,$window_file) or die "$!\n $window_file";
while(<FH>){
	chomp;
	next if $. == 1;
	my @line = split /\t/,$_;
	my $size= $line[-1];
	my $premz = $line[-2];
	my $start = $line[1];
	my $end = $line[2];
	if($size >100){$size = 100;}
	$num_of_windows = $line[0];
	$SWATH_window_count{$num_of_windows}{lower} = $start;
	$SWATH_window_count{$num_of_windows}{higher} = $end;
	$SWATH_window_count{$num_of_windows}{center} = $premz;
	$SWATH_window_count{$num_of_windows}{size} = $size;
	$SWATH_window{$premz}{window} = $num_of_windows;
	
	$SWATH_window{$premz}{size} = $size;
	$SWATH_window{$premz}{lower} = $start;
	$SWATH_window{$premz}{higher} = $end;
	$isolationlist{$premz} = $line[0];

}
close(FH);
#print Dumper \%SWATH_window_count; exit;

#my %ms1toms2;
my $cycle;
my $index;
my $experiment;
my $isolationwindow;
# READ mzML file
print "*****Read mzML file: $mzML_file!\n";
my %ms1scan_start;
my %ms1scan_cnt_index;
my %ms1scan_basepeak;
my %ms2scan_last_index;
my %ms2scan_next_index;
my %premz_index_ms2scan;
my $mslevel;
open(FH,$mzML_file) or die "$! $mzML_file\n";
while(<FH>){
	if(m{<spectrum index="(\d+)" id="sample=[^ ]+ period=[^ ]+ cycle=(\d+) experiment=(\d+)"}){
		($index, $cycle,$experiment) = ($1,$2,$3-1);
	}elsif(/name="ms level" value="(\d+)"/){
		$mslevel = $1;
		if($mslevel == 2){
			$ms2scan = $index;
		}else{
			$ms1scan++;
		}
	}elsif(m{accession="MS:1000827" name="isolation window target m/z" value="([\d.]+)" }){
		my $premz_ = $1;
		$ms2toms1{$ms2scan} = [sprintf("%06.f",$ms1scan),$premz_];
		printf "\r$ms1scan\t$ms2scan\t$index";
		$ms2scan_last_index{$index} = $premz_index_ms2scan{$premz_} if exists $premz_index_ms2scan{$premz_};
		$ms2scan_next_index{$premz_index_ms2scan{$premz_}} = $index if exists $premz_index_ms2scan{$premz_};
		
		$premz_index_ms2scan{$premz_} = $ms2scan;
	}elsif(m{accession="MS:1000505" name="base peak intensity" value="([^"]+)"} && $mslevel == 1){
		$ms1scan_basepeak{$ms1scan} = $1;
	}
}
close(FH);


my @ms1_charges;
my %charge_assaignment;
my %ms1_file_position;
my %ms2_file_position;
my %pseudomgffilehandle_position;
my %title_to_mgffile;
my %factorial_buffer;

mkdir "${dir}/$ARGV[0]";
my %filehandle;
my $line_block;
my $counter_currernt;


foreach my $key (1..$num_of_windows){
	open($filehandle{ms1}{$key},">${dir}/$ARGV[0]/index${key}_$ms1_file") or die "Can't write mgf file! ${dir}/$ARGV[0]/index${key}_$ms1_file $!\n";
	
	open($filehandle{ms2}{$key},">${dir}/$ARGV[0]/index${key}_$mgf_file") or die "Can't write mgf file! ${dir}/$ARGV[0]/index${key}_$mgf_file $!\n";
	
	foreach my $pseudo_mgffile (@pseudomgffiles){
		$pseudo_mgffile =~ s/${dir}\///;
		open($filehandle{pseudo_mgf}{$pseudo_mgffile}{$key},">${dir}/$ARGV[0]/pseudo_index${key}_$pseudo_mgffile") or die "Can't write mgf file! ${dir}/$ARGV[0]/pseudo_index${key}_$pseudo_mgffile $!\n";
	}
	
	foreach my $pepxmlfile(@pepxmlfiles){
		$pepxmlfile =~ s/${dir}\///;
		open($filehandle{pepxml}{$pepxmlfile}{$key},">${dir}/$ARGV[0]/index${key}_$pepxmlfile") or die "Can't write mgf file! ${dir}/$ARGV[0]/index${key}_$pepxmlfile $!\n";
	}
}



# read pepXML file;
#print "read pepXML file: $pepxmlfile\n";
my $end_scan;
my $retention_time_sec;
my $precursor_intensity;
my $peptide_prev_aa;
my $peptide_next_aa;
#my $alternative_protein;
my $num_matched_ions;
my $tot_num_ions;
my $num_missed_cleavages;
my $calc_neutral_pep_mass;
my $protein;
my $start_scan;
my $assumed_charge;
my $spectrum;
my $precursor_neutral_mass;
my $peptide;
my $massdiff;
my $mod_aminoacid_mass;
my $mod_aminoacid_position;
my %mod_aminoacid;
my $modified_peptide;
my $strip_peptide;
my $hyperscore;
my $nextscore;
my $expect;
my $search_data_line;
my $pepLen;
my $hit_rank;
my %search_res;
my $data_index;
my @ppms;
my %search_res;
my %multiple_charge_hits_rm;
my %ms1scan_clean_data;
my %xx;
my %uniq_pl;
my %keep_specs;
my %iso_distribution_buffer;
my $modification_tag;

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
my $cycle;

my $n;
my $tell_pos_last;


print "\n pepXML file spliting!\n";

my %pepXMLfilehandles;
foreach my $pepxmlfile(@pepxmlfiles){
	$n++;
	print "$pepxmlfile\n";
	open($pepXMLfilehandles{$pepxmlfile},"$dir/$pepxmlfile") or die "no $dir/$pepxmlfile $!\n";
	my $fh = $pepXMLfilehandles{$pepxmlfile};
	
	my $line_block;
	my $pepxml_ms2scan;
	
	while(<$fh>){
		chomp;
		my $line = $_;
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

		
		$line_block .= $line."\n";
		#print "$line_block\n-###-\n";
		if($line =~ /^<\/search_summary/){ #print "$line_block\n---\n";
			if($n == 1){
				foreach my $key (1..$num_of_windows){
					foreach my $pepxmlfile (@pepxmlfiles){
						$pepxmlfile =~ s/${dir}\///;
						my $fh_ = $filehandle{pepxml}{$pepxmlfile}{$key};
						print $fh_ $line_block,"\n";
					}
				}
			}
			$line_block = "";
		}
		if($line =~ /<\/spectrum_query/){#print "$line_block\n@@@@\n";
			if($line_block =~ /spectrum="[^."]+\.([^."]+)\.([^."]+)\.[^."]+"/){
				$pepxml_ms2scan = $1;
				#print "@@\n";
			}
			my $pepxml_premz = $ms2toms1{$pepxml_ms2scan}->[1];
			printf "\r$pepxml_ms2scan";
			my $fh_ = $filehandle{pepxml}{$pepxmlfile}{$SWATH_window{$pepxml_premz}{window}};
			print $fh_ $line_block;
			$line_block = "";
		}	
	}
}



print "\n pseudo mgf file spliting!\n";

my %pseudomgffilehandles;
foreach my $pseudomgffile(@pseudomgffiles){
	print $pseudomgffile,"\n";
	open($pseudomgffilehandles{$pseudomgffile},"${dir}/$pseudomgffile") or die "no ${dir}/$pseudomgffile $!\n";
	my $fh = $pseudomgffilehandles{$pseudomgffile};
	my $line_block;
	my $premz;
	while(<$fh>){
		chomp;
		my $Line = $_;
		printf "\r$Line" if $Line =~ /TITLE/;
		$line_block .= $_."\n";
		if (/^PEPMASS=([\d.]+)/){
			$premz = $1;
		}
		if(/END/){
			my $fh =  $filehandle{pseudo_mgf}{$pseudomgffile}{$SWATH_window{$premz}{window}};
			#print STDERR join "\t",$premz,$SWATH_window{$premz}{window},$fh,"\n";
			print $fh $line_block;
			$line_block = "";
		}
	}
}


foreach my $key (1..$num_of_windows){
	foreach my $pseudo_mgffile (@pseudomgffiles){
		$pseudo_mgffile =~ s/${dir}\///;
		close($filehandle{pseudo_mgf}{$pseudo_mgffile}{$key});
	}
	foreach my $pepxmlfile(@pepxmlfiles){
			$pepxmlfile =~ s/${dir}\///;
			close($filehandle{pepxml}{$pepxmlfile}{$key});
	}
	
}



print "\nMS1 file spliting!\n";
# MS1 files

open(FH,$ms1_file) or die "$! $ms1_file \n";
while(<FH>){
	chomp;
	my $Line = $_;
	if($Line =~ /^([SI])/){
		printf "\r$Line";
		foreach my $key (1..$num_of_windows){
			my $fh =  $filehandle{ms1}{$key};
			print $fh $Line,"\n";
		}
		$counter_currernt = 1;

	}
	if(/^[0-9]/){
		my @line = split / /,$Line;
		#printf $Line,"\n";
		my $fh =  $filehandle{ms1}{$counter_currernt};
		#print $fh $Line,"    $counter_currernt---$DIA_window_count{$counter_currernt}{lower}---$DIA_window_count{$counter_currernt}{higher}\n";
		
			
		if($line[0] > $SWATH_window_count{$counter_currernt}{higher}  ){
			#$counter_currernt = $counter_next;
			#print "$line[0] > $DIA_window_count{$counter_currernt}{higher},\n";
			while($line[0] > $SWATH_window_count{$counter_currernt}{higher} ){
				last if $counter_currernt >= $num_of_windows;
				$counter_currernt++;
			}
		}
		if($line[0] > $SWATH_window_count{$counter_currernt}{lower}  && $line[0] < $SWATH_window_count{$counter_currernt}{higher} ){
			my $fh =  $filehandle{ms1}{$counter_currernt};
			print $fh $Line,"\n";
			if( $counter_currernt < $num_of_windows){
				my $fh_1 =  $filehandle{ms1}{$counter_currernt+1};
				print $fh_1 $Line,"\n";
			}
			if( $counter_currernt > 1){
				my $fh_2 =  $filehandle{ms1}{$counter_currernt-1};
				print $fh_2 $Line,"\n";
			}
		}

		
	}
# S       000002  000002
# I       RTime   3.527
# 371.100725948884 153.93243
# 372.101898686153 57.99728
# 373.099650735678 38.47581
# 374.100456055843 10.70806
# 375.092751978369 3.28436
# 375.851551494201 4.19612
# 380.900989420792 5.18557
# 391.285290185289 19.50571
# 391.852436861849 3.43790
# 392.29003339664 4.15411

}
close(FH);
foreach my $key (1..$num_of_windows){
	close($filehandle{ms1}{$key});
}



print "\nMGF file spliting!\n";
# mgf files

my $premz;
open(FH,$mgf_file) or die "$! $mgf_file \n";
while(<FH>){
	
	chomp;
	my $Line = $_;
	printf "\r$Line" if $Line =~ /TITLE/;
	$line_block .= $_."\n";
	if (/^PEPMASS=([\d.]+)/){
		$premz = $1;

	}
	if(/END/){
		my $fh =  $filehandle{ms2}{$SWATH_window{$premz}{window}};
		
		print $fh $line_block;
		$line_block = "";
	}
	
}
close(FH);
foreach my $key (1..$num_of_windows){
	close($filehandle{ms2}{$key});
}



my $window_info = [\%SWATH_window,\%isolationlist,$num_of_windows];
my $mzML_info = [\%ms1scan_start,\%ms1scan_cnt_index,\%ms1scan_basepeak,\%ms2scan_last_index,\%ms2scan_next_index,\%premz_index_ms2scan,\%ms2toms1];

print "\nstart to parallal fork! \n";
my $pm = new Parallel::ForkManager($max_processes);  

my $i=1;

for($i=1;$i<=$num_of_windows;$i++){
	$pm->start and next;

	my $index = $i;
	print "\nParallel_index: ",$index,"\n";
	data_refinement(
	$filename,
	$dir,
	$ms1_ppm_tolerence,
	$ms2_ppm_tolerence,
	$ptm,
	$ptmfind,
	$fastafile,
	$index,
	$window_info,
	$mzML_info
	);

	$pm->finish;
}

$pm->wait_all_children;

sub data_refinement {

	SWATH_pepXML_to_pin_inparallel::pepXML_refinement(@_);
}

my @pre_ppms;
my @frag_ppms;
my %keep_pep_list;
my %keep_pep_list_expect_strip_pep;
my @tmpfiles = glob "${dir}/$filename/*_masscal_tmp.$filename";

open(TMPOUT,">${dir}/masscal_tmp.$filename");
foreach my $f (@tmpfiles){
	open(tmpf,$f);
	while(<tmpf>){

		my @lines = split /\t/,$_;
		my $expect = $lines[2];
		my $scan = $lines[0];
		my $pep = $lines[1];
		my $strip_pep = $lines[6];
		my $ppm_ = $lines[8];
		
		push @{$keep_pep_list{$scan}},$pep;
		$keep_pep_list_expect_strip_pep{$pep} = [$expect,$strip_pep];
		my $ref1 = [split /#/,$lines[9]];
		my $iso = $ref1->[2];
		my $ref3 = [split /#/,$lines[11]];
		my $frag_ppm_ = $ref3->[3];	
		my $protein = $lines[7];
		if ($expect <= 0.01 && $protein !~ /REV_/ && $iso ){
			push @pre_ppms, $ppm_ ;
			push @frag_ppms, $frag_ppm_ ;
		}
		print TMPOUT $_;
	}
	close(tmpf);
}



my $ppm_ave;
my $frag_ppm_ave;
my $ppm_sd;
my $frag_ppm_sd;

#print Dumper \@ppms;exit;

if (scalar @pre_ppms){
	$ppm_ave = ave(\@pre_ppms);
	$ppm_sd = stdev(\@pre_ppms,$ppm_ave);
	printf ("\n\033[31m%s\t%s\t\033[32m%s\033[0m\n", $filename,"pre_ppm_ave: $ppm_ave","pre_ppm_sd: $ppm_sd");
	print logfh join "\t", $filename,"pre_ppm_ave: $ppm_ave","pre_ppm_sd: $ppm_sd\n";
	if(scalar @pre_ppms >500 && ! $ptm){
		my ($prob1,$prob2);
		($ppm_ave,$ppm_sd,$prob1,$prob2) = MassAccuracy::MassAccuracy(\@pre_ppms);
		printf ("\n\033[31m%s\t%s\t\033[32m%s\033[0m\n", $filename,"pre_ppm_ave: $ppm_ave","pre_ppm_sd: $ppm_sd");
		print logfh join "\t", $filename,"pre_ppm_ave: $ppm_ave","pre_ppm_sd: $ppm_sd\n";
		print "prob1; $prob1; prob2: $prob2\n";
		print logfh "prob1; $prob1; prob2: $prob2\n";
	}else{
		print "scalar \@pre_ppms = ",(scalar @pre_ppms),"\n";
		print logfh "scalar \@pre_ppms = ",(scalar @pre_ppms),"\n";
	}
}


if (scalar @frag_ppms){
	$frag_ppm_ave = ave(\@frag_ppms);
	$frag_ppm_sd = stdev(\@frag_ppms,$frag_ppm_ave);;
	printf ("\n\033[31m%s\t%s\t\033[32m%s\033[0m\n", $filename,"frag_ppm_ave: $frag_ppm_ave","frag_ppm_sd: $frag_ppm_sd");
	printf logfh join "\t",$filename,"frag_ppm_ave: $frag_ppm_ave","frag_ppm_sd: $frag_ppm_sd\n";
	if(scalar @frag_ppms > 500 && ! $ptm){
		my ($prob1,$prob2);
		($frag_ppm_ave,$frag_ppm_sd,$prob1,$prob2) = MassAccuracy::MassAccuracy(\@frag_ppms);
		printf logfh join "\t", $filename,"frag_ppm_ave: $frag_ppm_ave","ppm_sd: $frag_ppm_sd\n";
		print "prob1; $prob1; prob2: $prob2\n";
		print logfh "prob1; $prob1; prob2: $prob2\n";
	}else{
		print "scalar \@pre_ppms = ",(scalar @frag_ppms),"\n";
		print logfh "scalar \@pre_ppms = ",(scalar @frag_ppms),"\n";
	}
}
#print STDERR Dumper \%search_res;

my $yy;
my %pepNum;
my %keep_;
my %cmped_list;

open(TMP_INPUT,"${dir}/tmp.$filename");
while(<TMP_INPUT>){
	chomp;
	my $line = $_;
	my @lines = split /\t/,$_;
	my $scan = $lines[0];
	printf "\r$scan";
	my $pep = $lines[1];
	my $protein = $lines[7];
	my $expect = $lines[2];
	my $strip_pep = $lines[6];
	my $charge = $lines[4];
	$pepNum{$filename}{raw}{$protein =~/REV_/?"Decoy":"Target"}++ if($expect <0.01 and !$uniq_pl{raw}{$pep}++);
	$yy->{raw}++;
	my $sp = $lines[3];
	my $keep = 1;	
	
	my $ppm_ = $lines[8];
	my $ref1 = [split /#/,$lines[9]];
	my $iso = $ref1->[2];
	my $ref2 = [split /#/,$lines[10]];
	my $ref3 = [split /#/,$lines[11]];
	my $mod_name = $ref3->[6];
    my $frag_ppm_ = $ref3->[3];	

	if(abs($frag_ppm_ - $frag_ppm_ave) < 3*$frag_ppm_sd && abs($ppm_ave-$ppm_) < 3 * $ppm_sd  ){
		if($mod_name){
			if(abs($ppm_ - $ppm_ave) > $ppm_sd){
				$keep = 0;
			}
		}
		if($keep){
			foreach my $pep_ (grep{$_ ne $pep && (! $cmped_list{$scan}{$_}) }@{$keep_pep_list{$scan}}){
				my $i;
				my $s;
				my $strip_pep_ = $keep_pep_list_expect_strip_pep{$pep_}->[1];
				$i = ($strip_pep =~ m{(^${strip_pep_}.+)|(.+${strip_pep_}$)});
				$keep = 0 if $i;
				if(! $ptm){
					$s = seq_similar($strip_pep,$strip_pep_);
					if($s){
						if($keep_pep_list_expect_strip_pep{$pep_}->[0] < $expect){
							$keep = 0;
						}
					}
				}
				if($strip_pep eq $strip_pep_  && $keep_pep_list_expect_strip_pep{$pep_}->[0] < $expect){
					$keep = 0;
				}
			
			}
		}
	}else{
		$keep = 0;
	}

	if($keep  ){
		$yy->{keep}++;
		my ($x,$y,$z) = @lines[3,4,5];
		$pepNum{$filename}{keep}{$protein =~ /REV_/?"Decoy":"Target"}++ if($expect < 0.01 and !$uniq_pl{keep}{$pep}++);
		$keep_{$x}{$y}{$z} = 1000 if ! defined $keep_{$x}{$y}{$z};
		if($keep_{$x}{$y}{$z} > $expect){
			$keep_{$x}{$y}{$z} = $expect;
			$keep_specs{$x}{$y}{$z} = [$ref1,$ref2,$ref3];
		}
	}else{
		$cmped_list{$scan}{$pep} = 1;
		$yy->{del}++;
		$pepNum{$filename}{del}{$protein =~ /REV_/?"Decoy":"Target"}++ if($expect < 0.01 and !$uniq_pl{del}{$pep}++);
	}
}
close(TMP_INPUT);
#unlink "${dir}/tmp.$filename";
print Dumper \%pepNum;

my %peptide_with_aa_mutation;

foreach  my $title(keys %keep_specs){
	foreach my $charge (keys %{$keep_specs{$title}}){
		foreach my $exp_mz (keys %{$keep_specs{$title}{$charge}}){
			my $iso = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[2];
			my $pep = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[3];
			my $expect = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[4];
			my $protein = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[5];
			my $strip_pep = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[6];
			my $complete_peptide =  $keep_specs{$title}{$charge}{$exp_mz}->[0]->[7];
			my $centermz =  $keep_specs{$title}{$charge}{$exp_mz}->[0]->[8];
			my $window = $isolationlist{$centermz};
			my $mod_name = $keep_specs{$title}{$charge}{$exp_mz}->[2]->[6];
			
			if($mod_name =~ /(...)[-][>](...)\[(.)\]/){
				$peptide_with_aa_mutation{$pep}{$mod_name} = 1;
			}
			$pepNum{$filename}{res}{$protein=~/REV_/?"Decoy":"Target"}++ if($expect<0.01 and !$uniq_pl{res}{$pep}++);
		}
	}
}

# retrieve the mutated peptide candidate using fasta 
my %replaced_peptide_with_aa_mutation;
if($ptmfind){
print "\n refine mutated peptide candidates by retrieve proteins in fasta\n";
foreach my $pep (keys %peptide_with_aa_mutation){
	foreach my $mod_name (keys %{$peptide_with_aa_mutation{$pep}}){ 
		$mod_name =~ /(...)[-][>](...)\[(.)\]/;
		my $aa = $3;
		my $replaced_aa = $aa_short{$2};
		#print join "\t",$1,$aa_short{$1},$2,$aa_short{$2},$3;
		#print "\n";
		my $pep_ = $pep;
		if($replaced_aa eq "I|L"){
			$pep_ =~ s/(.+)${aa}\[[^][]+\](.+)/${1}L${2}/;
			$pep_ =~ s/\[[^][]+\]//g;
			$replaced_peptide_with_aa_mutation{$pep_} = $pep;
			$pep_ = $pep;
			$pep_ =~ s/(.+)${aa}\[[^][]+\](.+)/${1}L${2}/;
			$pep_ =~ s/\[[^][]+\]//g;
			$replaced_peptide_with_aa_mutation{$pep_} = $pep;
		}else{
			$pep_ =~ s/(.+)${aa}\[[^][]+\](.+)/$1${replaced_aa}$2/;
			$pep_ =~ s/\[[^][]+\]//g;
			$replaced_peptide_with_aa_mutation{$pep_} = $pep;
		}
	}

}
}

my $id;
my %retrieved_list;
my %retrieved_done_list;
if($ptmfind){
#print Dumper \%replaced_peptide_with_aa_mutation;
open(FASTA,"${dir}/$fastafile") or die "$! no $fastafile \n";
while(<FASTA>){
	my $line = $_;
	if($line =~ />/ && $line !~ /REV_/){
		my $id = $line;
		$id =~ s/>//;
		chomp($id);
		my $seq;
		($seq = <FASTA>);
	
		foreach my $pep_ (keys %replaced_peptide_with_aa_mutation){
			if($seq =~ /$pep_/ && ! exists $retrieved_list{$pep_}){
				$retrieved_list{$replaced_peptide_with_aa_mutation{$pep_}} = [$pep_,$id];
				$retrieved_done_list{$pep_} = 1;
			}
		}
	}
}
close(FASTA);
#print Dumper \%retrieved_list;
}

my @files = glob "$dir/$filename/*";
unlink $_ for @files;

rmdir "$dir/$filename";

my %outfile1;
my %outfile2;
my %m;
my %iontype;
my %uniq_pl;

open(pinout,">${dir}/${filename}.pin") or die  "${dir}/${filename}.pin $!";
open(mgfout,">${dir}/$outfilename") or die  "${dir}/$outfilename $!";
open(PTMOut,">${dir}/${filename}.PTM.tsv") or die  "${dir}/${filename}.PTM.tsv $!";
print pinout join "\t", qw/SpecId	Label	ScanNr	ExpMass	CalcMass	hyperscore deltahyperscore lnExpect	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	enzInt	dM	absdM	Iso1 Iso2 Iso3 Iso_ge_4 LnPreInt RTINSECONDS PreSN/;
print PTMOut join "\t", qw/SpecId	Label	ScanNr	ExpMass	CalcMass	hyperscore deltahyperscore lnExpect	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	enzInt	dM	absdM	Iso1 Iso2 Iso3 Iso_ge_4 LnPreInt RTINSECONDS PreSN/;
print pinout "\t";
print PTMOut "\t";
if(scalar @PTM_variable){
	print pinout join "\t",@PTM_variable; 
	print PTMOut join "\t",@PTM_variable; 
	print pinout "\t";
	print PTMOut "\t";
}

print pinout join "\t",qw /Undefined_mod Y_ions n_Y_ions	B_ions n_B_ions matched_coverage n_matched_coverage matched_int matched_int_sd FragAveppm abs_FragAveppm FragSdppm Peptide	Proteins/;
print PTMOut join "\t",qw /Undefined_mod Y_ions n_Y_ions	B_ions n_B_ions matched_coverage n_matched_coverage matched_int matched_int_sd FragAveppm abs_FragAveppm FragSdppm PTMScore Peptide	Proteins	PhosScore	PhosProbs PhosphoSiteloc Phosphosite_loc_or_not PTMSites /;
print PTMOut "\n";
print pinout "\n";

my %m_;
my @pseudomgffiles = glob "$dir/w_[0-9]*_".$filename.".mgf";
foreach my $mgffile(@pseudomgffiles){
	open(MGF,$mgffile) or die "$dir/$mgffile $!\n";
	my $p = 0;
	my $title;

	my $charge;
	my $line1;
	my $line2;
	my $fh1;
	my $fh2;
	my @ion;

	while(<MGF>){
		my $txt = $_;
		if(/BEGIN/../RTINSECONDS/){
			$line1 .= $txt;
			if($txt =~ /TITLE=([^.]+\.[^.]+\.[^.]+)\..*/){
				$title = $1;
				#print "$title\n";
			}
		}
		
		if(/PEPMASS/../END/){
			next if /CHARGE|PEPMASS/;
			$line2 .= $txt;
			if($txt =~ /^[0-9]/){
				chomp($txt);
				push @ion,[split / /,$txt];
			}
			if(/END/){
				if(exists $keep_specs{$title}){
					printf "\r$title";
					foreach my $charge (keys %{$keep_specs{$title}}){
						if(exists $keep_specs{$title}{$charge}){
						foreach my $exp_mz (keys %{$keep_specs{$title}{$charge}}){	
							if(exists $keep_specs{$title}{$charge}{$exp_mz} && scalar @{$keep_specs{$title}{$charge}{$exp_mz}->[0]} 
								&& scalar @{$keep_specs{$title}{$charge}{$exp_mz}->[1]} && 
								scalar @{$keep_specs{$title}{$charge}{$exp_mz}->[2]}){
							my $iso = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[2];
							my $cal_premz =  $keep_specs{$title}{$charge}{$exp_mz}->[0]->[0];

							my $pep = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[3];
							my $expect = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[4];
							my $protein = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[5];
							my $strip_pep = $keep_specs{$title}{$charge}{$exp_mz}->[0]->[6];
							my $complete_peptide =  $keep_specs{$title}{$charge}{$exp_mz}->[0]->[7];
							my $centermz =  $keep_specs{$title}{$charge}{$exp_mz}->[0]->[8];
							my $window = $isolationlist{$centermz};

							my ($matched_y,$matched_b);
							my ($matched_mzs_type,$matched_int,$matched_int_sd,$ppm_m,$match_type,$matched_coverage,$ppm_s,$ppm);
							my ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$mod_name,$PTMSiteScore,$PTMSiteProb,$unimod_mass);
							my $mod_info = $keep_specs{$title}{$charge}{$exp_mz}->[2];
								($matched_y,$matched_b,$matched_int_sd,$ppm_m,$matched_int,$matched_coverage,$mod_name,$PTMSiteScore,$PTMSiteProb,$unimod_mass,$ppm_s,$ppm) = @{$mod_info};
							my $deisotope_peaks = deisotope([@ion]);

							my $peptide_mod = $pep;#print STDERR "\r$title\t$peptide_mod\n";
							if(exists $retrieved_list{$peptide_mod}){
								$protein = $retrieved_list{$peptide_mod}->[1];
								$peptide_mod = $retrieved_list{$peptide_mod}->[0];
								$complete_peptide =~ /(^..)(.+)(..$)/;
								$complete_peptide = $1.$peptide_mod.$3;
								($mod_name,$PTMSiteScore,$PTMSiteProb,$unimod_mass) = (undef, undef, undef, undef);
							}
							my $PhosSite_num = ($peptide_mod =~ s/\[166.9984\]|\[181.0140\]|\[243.0297\]//g);
							#print Dumper $peptide,$peptide_mod, $site_num;
							#print STDERR "\r$title\t$peptide_mod,  $PhosSite_num,  $expect\n";
							my ($Phosscore,$PhosSiteProb,$PhosphoSiteloc,$phosphosite_loc_or_not);
							if($PhosSite_num > 0){

								my $Site = "STY";
								printf "\r$title\t$peptide_mod";
								($peptide_mod,$Phosscore,$PhosSiteProb,$PhosphoSiteloc,$phosphosite_loc_or_not) = PhosphoSiteScore($peptide_mod,$charge,$Site,$deisotope_peaks,$PhosSite_num,$strip_pep);

								my $pep_for_fragmentation = $peptide_mod;
								$pep_for_fragmentation =~ s/S\[Phos\]/S[166.9984]/g;
								$pep_for_fragmentation =~ s/T\[Phos\]/T[181.0140]/g;
								$pep_for_fragmentation =~ s/Y\[Phos\]/Y[243.0297]/g;
								my $theoretical_mz_table = peptide::fragmentation($pep_for_fragmentation,$charge);
								($matched_mzs_type,$matched_int,$matched_int_sd,$ppm_m,$ppm_s) = ms2spectra_match($theoretical_mz_table,[@ion]);
								($match_type,$matched_coverage) = tag($matched_mzs_type);
								($matched_y,$matched_b) = ($match_type->{"Y"},$match_type->{"B"});
								
								$complete_peptide =~ /(^..)(.+)(..$)/;
								$peptide_mod = $1.$peptide_mod.$3;
								#print STDERR "pep:","\t",$complete_peptide,"\t",$peptide_mod,"\t",$Phosscore||"-","\t",$PhosSiteProb||"-","\n";
							}else{
								$peptide_mod = $complete_peptide;
							}
								$m_{$filename}{2}{$protein=~/REV_/?"Decoy":"Target"}++ if($expect<0.01 and !$uniq_pl{2}{$pep}++);
								   # print "$ppm_m\t$matched_int_sd\n";
									if($ppm_m ne "NA" && $matched_int_sd ne "NA"){
										$m_{$filename}{1}{$protein=~/REV_/?"Decoy":"Target"}++ if($expect<0.01 and !$uniq_pl{1}{$pep}++);
										if($PhosSiteProb || $mod_name){
											print PTMOut join "\t", @{$keep_specs{$title}{$charge}{$exp_mz}->[1]};
											print PTMOut "\t";
										}else{
											print pinout join "\t", @{$keep_specs{$title}{$charge}{$exp_mz}->[1]};
											print pinout "\t";
										}
										
										
										if(scalar @PTM_variable){
											if($PhosSiteProb || $mod_name){
												print PTMOut join "\t",map{$PTM_Types{$_}{$complete_peptide} + 0}@PTM_variable;
												print PTMOut "\t";
											}else{
												print pinout join "\t",map{$PTM_Types{$_}{$complete_peptide} + 0}@PTM_variable;
												print pinout "\t";
											}
										}
										if($PhosSiteProb || $mod_name){
											print PTMOut join "\t",( $mod_name ? 1:0), $matched_y+0,($matched_y+0)/length($strip_pep),$matched_b+0,($matched_b+0)/length($strip_pep),$matched_coverage,($matched_coverage+0)/length($strip_pep),$matched_int ? log($matched_int):0,$matched_int_sd,$ppm_m,abs($ppm_m),$ppm_s,$PTMSiteScore||"-", $peptide_mod.( $mod_name ? ":(".$mod_name."[".$unimod_mass."]".")":""), $protein,$Phosscore||"-",$PhosSiteProb||"-",$PhosphoSiteloc,$phosphosite_loc_or_not,($PTMSiteProb||"-");
											print PTMOut "\n" ;
										}else{
											print pinout join "\t",( $mod_name ? 1:0), $matched_y+0,($matched_y+0)/length($strip_pep),$matched_b+0,($matched_b+0)/length($strip_pep),$matched_coverage,($matched_coverage+0)/length($strip_pep),$matched_int ? log($matched_int):0,$matched_int_sd,$ppm_m,abs($ppm_m),$ppm_s, $peptide_mod.( $mod_name ? ":(".$mod_name."[".$unimod_mass."]".")":""), $protein;
											print pinout "\n" ;
										}
										
									}
									
									if($expect < 1 ){
										if ($PhosSiteProb){$mod_name = "Phos"};
										print mgfout $line1;
										print mgfout "CalPEPMASS=$keep_specs{$title}{$charge}{$exp_mz}->[0]->[0] $keep_specs{$title}{$charge}{$exp_mz}->[0]->[1]\n";
										print mgfout "PEPMASS=$exp_mz $keep_specs{$title}{$charge}{$exp_mz}->[0]->[1]\n";
										print mgfout "CHARGE=${charge}+\n";
										print mgfout "MOD=$mod_name\n" if($PhosSiteProb || $mod_name);;
										print mgfout $line2;
									}
								#}
							
							#}
						}
						}
						}
					}
				}
				$line1 = "";
				$line2 = "";
				undef @ion; 
			}
			
		}
	}
	close(MGF);
}

# foreach my $iso (2..$min_n_iso_required){
	# my $fh1 = $outfile1{$iso};
	# my $fh2 = $outfile2{$iso};
	# close($fh1);
	# close($fh2);
# }

# foreach my $pseudomgffile(@pseudomgffiles){
	# my $fh = $pseudomgffilehandles{$pseudomgffile};
	# close($fh);
# }
close(mgfout);
close(pinout);
close(PTMOut);
close(MS2);
close(MS1);


print "\nDone\n";



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


sub stdev{ # Standard deviation calculation
	my ($data,$ave_) = (@_);
	if($#$data){
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
			#print $mod_title," ",$neutral_loss,"\n";
			push @tmp,[$site,$position];
		}
		if(m{<umod:NeutralLoss mono_mass="([^"]+)"}){
			$neutral_loss{$site}{$position} = $1 if $1 != 0;
			
		}
		if(m{<umod:delta mono_mass="([^"]+)"}){
			$mod_delta = $1;
			#print $mod_title," ",$neutral_loss,"\n";
			push @mod_mass,$mod_delta if !$uniq{$mod_delta}++;
			# push @{$mods{$mod_delta}},$mod_title;
			foreach my $i (@tmp){
				($site,$position) = @{$i};
				push @{$mods{$mod_delta}{$site}{$position}},[$neutral_loss{$site}{$position},$mod_title];
			}
			undef @tmp;
			undef %neutral_loss;
			#$neutral_loss = "";
			#print join "\t",$mod_delta,$site,$position,$mod_title,"\n";
		}
	}
	close(MOD);
	#exit;
	#print Dumper \%mods;exit;
	return [\%mods,\@mod_mass];
}



sub iso_peak_check{
	my ($prems1scan,$cal_premz,$assumed_charge,$iso_require,$ms2scan,$peptide) = @_;
	
	my ($ms1scan_mz_list,$lowest_sig,$highest_sig);
	if(exists $ms1scan_clean_data{$prems1scan}{$ms2scan}){
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = @{$ms1scan_clean_data{$prems1scan}{$ms2scan}};
	}else{
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = MS1file($prems1scan,$ms2scan);
		$ms1scan_clean_data{$prems1scan}{$ms2scan} =[$ms1scan_mz_list,$lowest_sig,$highest_sig];
	}
	my $len = $#$ms1scan_mz_list; #print Dumper $len;
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
	
	return ($tag_iso, $obs_mono_mz, $obs_mono_int,$measured_ppm_mono,$measured_pre_SN) ;
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
	
	# $ms2scan_last_index{$scan} = $ms2scan
	# $ms2scan_next_index{$scan} = $ms2scan
	# $ms2toms1{$ms2scan} = [sprintf("%06.f",$ms1scan),$premz_];
	
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
			#$prems1scan = $ms2toms1{$ms2scan}->[0];
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
		#my $cmp_spec = $ms2spectra_data{$premz}{$ind};
		$com_mz = Spectrum_cmp($curr_ms1scan_mz_list,$cmp_spec,$com_mz);
		#print STDERR Dumper '##ind\t'.$ind,$cmp_spec->[1];
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
			#$ms2peaks = readMS2(sprintf("%06.f",$rms2scan));
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
			#$ms2peaks = readMS2(sprintf("%06.f",$rms2scan));
		}else{
			$ms2scan = $ms2scan_next_index{$ms2scan};
		}
		$ms2peaks = readMS2(sprintf("%06.f",$ms2scan));
		push @ms2_data_list,$ms2peaks;
	}
	foreach my $cmp_spec (@ms2_data_list){
		#my $cmp_spec = $ms2spectra_data{$premz}{$ind};
		$com_mz = Spectrum_cmp($currpeaks,$cmp_spec,$com_mz);
		#print STDERR Dumper '##ind\t'.$ind,$cmp_spec->[1];
	}
	
	return $com_mz;
}



sub MS1file{
	#** 需要input：指针位置，scan起始和终止，scan mz window **
	my ($ms1scan_,$ms2scan) = (@_);
	#print Dumper ($ms1scan_,$ms2scan);
	#$ms2toms1{$ms2scan} = [$ms1scan,$premz_]
	my $premz = $ms2toms1{$ms2scan}->[1];
	#my ($ms1scan_,$premz_) = (@{$ms2toms1{$ms2scan}});
	my ($mz_start,$mz_end) = ($SWATH_window{$premz}{lower} - 3, $SWATH_window{$premz}{higher} + 3);
	
	my $curr_ms1scan_pos = $ms1_file_position{sprintf("%06.f",$ms1scan_)};
	# my $curr_ms1scan_cnt = $ms1scan_cnt_index{$ms1scan_};
	# my $pre_ms1scan_pos  = $ms1_file_position{sprintf("%06.f",$ms1scan_ - 1)};
	# my $next_ms1scan_pos = $ms1_file_position{sprintf("%06.f",$ms1scan_ + 1)};

	my $cur_scan = $ms1scan_ + 0;
	my $pre_scan = $ms1scan_ - 1;
	my $next_scan = $ms1scan_ + 1;

	my ($curr_ms1scan_data,$curr_ms1scan_lowsig_data,$highest_sig) = readMS1($curr_ms1scan_pos,$mz_start,$mz_end);
	# my ($pre_scan_data,$pre_ms1scan_lowsig_data) = readMS1($pre_ms1scan_pos,$mz_start,$mz_end,$ms1scan_basepeak{$pre_scan});
	# my ($next_scan_data,$next_ms1scan_lowsig_data) = readMS1($next_ms1scan_pos,$mz_start,$mz_end,$ms1scan_basepeak{$next_scan});
	
	# my $com_mz;
	# $com_mz = Spectrum_cmp($curr_ms1scan_lowsig_data,$pre_scan_data,$com_mz);
	# $com_mz = Spectrum_cmp($curr_ms1scan_lowsig_data,$curr_ms1scan_lowsig_data,$com_mz);
	# $com_mz = Spectrum_cmp($curr_ms1scan_lowsig_data,$next_scan_data,$com_mz);
	
	
	# my $filtered_ms1scan_data;
	# my $lowest_sig = 10;
	# foreach my $ref (@{$curr_ms1scan_data}){
		# my ($mz,$i) = @{$ref};
		
		# if(defined $com_mz->{$mz}){
			# my @ints = (@{$com_mz->{$mz}});
			# if($ints[0] > 0 && $ints[2] > 0){
				# push @{$filtered_ms1scan_data},[$mz,$ints[1]];
				# $lowest_sig = $ints[1] if $lowest_sig > $ints[1];
			# }
		# }else{
			# push @{$filtered_ms1scan_data},$ref;
		# }
	# }
	
	return ($curr_ms1scan_data,$curr_ms1scan_lowsig_data, $highest_sig);	

}

sub readMS1{
	my ($pos,$mz_start,$mz_end) = (@_);	
	#print Dumper $bp;
	#my $pos = $ms1_file_position{$ms1scan_};
	#print join "#",$ms1scan_,$mz_start,$mz_end,$pos,"\n";exit;
	my @ms1_peaks_;
	my @ms1_low_sig;
	my $lowest_sig = 1000000 ;
	my $highest_sig = 0;
	seek(MS1,$pos,0);
	my $exit;
	while(<MS1>){
		#print $_;
		chomp;
		if(/^\d+/){
			my @line = split / /,$_;
			if ($line[0] >= $mz_start && $line[0] <= $mz_end){
				push @ms1_peaks_,[@line[0,1]] if $line[1] > $ms1_min_intensity;
				#push @ms1_low_sig,[@line[0,1]] if $line[1] < 10;
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
	#print Dumper \@ms1_peaks_,\@ms1_low_sig;exit;
	#return ([@ms1_peaks_],[@ms1_low_sig], $highest_sig);
	return ([@ms1_peaks_],$lowest_sig, $highest_sig);
}

sub readMS2{
	my ($ms2scan) = (@_);
	my @ms2_peaks_;
	my $pos = $ms2_file_position{$ms2scan};
	
	#print Dumper $ms2scan,$pos;
	seek(MS2,$pos,0);
	my $exit;
	#my $lowest_sig = 100000000;
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
	
	#print Dumper \%pseudomgffilehandle_position;
	#print Dumper $title;
	my @ms2_peaks_;
	my $pos = $pseudomgffilehandle_position{$title};
	my $fh = $pseudomgffilehandles{$title_to_mgffile{$title}};
	#print Dumper $fh,$pos,$title,$title_to_mgffile{$title};
	seek($fh,$pos,0);
	my $exit;
	#my $lowest_sig = 100000000;
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
				#print STDERR join "\t",$n,"@@",$mz1,$int1,"\n";#  if $n == 9;
				push @{$ref->{$mz1}},0;
				$sp1_index++;
				last if $sp1_index > $sp1_len;
			}
			last;
		}else{
			my ($mz1,$int1) = (@{$sp1->[$sp1_index]});
			my ($mz2,$int2) = (@{$sp2->[$sp2_index]}); 
			#print STDERR join "\t",$n ,"##",$mz1,$int1,$mz2,$int2,abs($mz1-$mz2),"\n";# if $n == 9;
			
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
	#print Dumper $res;
	return $ref;
}

sub boundary_index{
    my $data = shift;
	my $bp = shift;
	my $m = $#$data/2;
	#my $m = max_ind($data);
	
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
	#print join "\t",$start,$end,"\n";
    return ($start,$end);
}

sub find_mod(){
	#my ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs) = find_mod($prems1scan,$cal_premz,$assumed_charge,$start_scan,$peptide);
	my ($prems1scan,$cal_premz,$assumed_charge,$ms2scan,$pep,$pepLen,$strip_pep,$spectrum,$ms1profile,$ms2profile,$peptide_prev_aa,$peptide_next_aa) = (@_);
	#print Dumper \@_;exit;
	my $peaks ;
	my $deisotope_peaks; 
	
#if ($ms2scan == 45635){
	
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
	#print STDERR join "\t","\ntest:","prems1scan:$prems1scan,cal_premz:$cal_premz,assumed_charge:$assumed_charge,ms2scan:$ms2scan,pep:$pep,pepLen:$pepLen,strip_pep:$strip_pep","\n"; #exit;
	
	my ($ms1scan_mz_list,$lowest_sig,$highest_sig);
	if(exists $ms1scan_clean_data{$prems1scan}{$ms2scan}){
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = @{$ms1scan_clean_data{$prems1scan}{$ms2scan}};
	}else{
		($ms1scan_mz_list,$lowest_sig,$highest_sig) = MS1file($prems1scan,$ms2scan);
		$ms1scan_clean_data{$prems1scan}{$ms2scan} =[$ms1scan_mz_list,$lowest_sig,$highest_sig];
	}
	
	#my $ms1profile = MS1_profile($prems1scan,$ms2scan);
	#my $ms2profile = MS2_profile($ms2scan);
	
	my $len = $#$ms1scan_mz_list; #print Dumper $len;
	my $obs_mono_mz;
	my $obs_mono_int;
	my $tag_iso;
	my $measured_ppm_mono;
	my $corr_distribution;
	
	my $ms1_start_mz = $ms1scan_mz_list->[0]->[0];
	my $ms1_end_mz = $ms1scan_mz_list->[$len]->[0];
	#print  STDERR Dumper $ms1scan_mz_list;
	if(! defined $ms1_start_mz){
		return (0);
		next;
	}
	
	# print STDERR "MS1 start:$ms1_start_mz  MS1 end:$ms1_end_mz \n";
	
	#my @mod_mass_list = sort{$a <=> $b}grep {($cal_premz + $_/$assumed_charge) > $ms1_start_mz && ($cal_premz + $_/$assumed_charge) < $ms1_end_mz} (keys %unimods);
	#$peaks = readMS2(sprintf("%06.f",$ms2scan)) if ! defined $peaks;
	#print Dumper $ms2scan,sprintf("%06.f",$ms2scan),$peaks;
	#my $peaks = readMS2(sprintf("%06.f",$ms2scan+0));
	#print Dumper sprintf("%06.f",$ms2scan+0),$peaks;exit;
	#print Dumper $cal_premz,$assumed_charge,$ms1_start_mz,$ms1_end_mz,\@mod_mass_list;exit;
	
	#my $len_of_matched_by_without_mod;
	my $theoretical_mz_table = peptide::fragmentation($pep,$assumed_charge);
	my ($matched_mzs_type,$matched_int_without_mod,$matched_int_sd,$m,$match_mzs_topN5) = ms2spectra_match($theoretical_mz_table,$peaks);
	my ($match_tag,$match_coverage_without_mod) = tag($matched_mzs_type);
	#my  = tag($matched_mzs_type);
	#print STDERR Dumper  $matched_mzs_type,$matched_int_without_mod,$matched_int_sd,$m,$match_mzs_topN5,$match_tag,$match_coverage_without_mod;
	
	
	# my $pep_ = $pep;
	# $pep_ =~ s/\[[^][]+\]|c\([^()]+\)|n\([^()]+\)//g;
 
	my @theoretical_iso_distribution; 
	if(exists $iso_distribution_buffer{$strip_pep}){
		@theoretical_iso_distribution = @{$iso_distribution_buffer{$strip_pep}};
	}else{
		@theoretical_iso_distribution = @{IsotopeDistribution::CAPTT($strip_pep)};
		$iso_distribution_buffer{$strip_pep} =[@theoretical_iso_distribution];
	}
	#print Dumper \@theoretical_iso_distribution,$strip_pep;exit;
	
	# # print STDERR $match_mzs_topN5;
		# # map{print STDERR Dumper  $_, normalization_int($ms2profile->{$_}) }@{$match_mzs_topN5};
		# # print STDERR Dumper normalization_int($ms1profile->{$curr_mz});
	
	#print STDERR Dumper $ms2profile,$ms1profile;
	# Calculate average
	#*************************
	
	
	#my $s = tag($matched_mzs_type,"type");
	#print Dumper $match_coverage,$s;
	
	# print Dumper $theoretical_mz_table,$match_type,$matched_mzs_type; exit;
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
	#my ($matched_int_sd_mod,$ppm_ave_mod,$matched_int_reported,$match_type_report);
	my ($matched_int_sd_mod,$ppm_ave_mod,$ppm_sd_mod,$matched_int_reported,$match_type_report);
	my %last_candidate_trace_mz;
	my $pre_SN;
	
	while($index <= $len){
		
		#print Dumper $ms1scan_mz_list->[$index]->[0];
		my $curr_mz = $ms1scan_mz_list->[$index]->[0];
		my $curr_int = $ms1scan_mz_list->[$index]->[1];
		$pre_SN = $curr_int / (1 + $lowest_sig );
		# print STDERR "$curr_int / (1 + $lowest_sig );\n";
		# print STDERR "--check1 mz: $curr_mz pre_SN:$pre_SN\n";
		# print STDERR join "##", exists $trace_mz{$curr_mz}, exists $last_candidate_trace_mz{$curr_mz},$pre_SN, "\n";
		if (exists $trace_mz{$curr_mz} || exists $last_candidate_trace_mz{$curr_mz} || $pre_SN < 5 ){
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
			# foreach my $i (keys %tmp_trace_mz_high){
				# $trace_mz{$i} = 1;
			# }
			next;
		}

		my $non_zero_in_MS1;
		
		my @ms1arr = @{$ms1profile->{$curr_mz}};
		
		my @obs_iso_distribution;
		push @obs_iso_distribution,$curr_int;
		# add testing if the false positive hits by other charged precursor
		
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
		#my $isotopic_ints_highest_peak = $curr_int; 
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
		
		
		#print STDERR "check mz: $curr_mz\n";
		my @obs_mz;
		while($index_++ < $len){
			my $curr_mz_ = $ms1scan_mz_list->[$index_]->[0];
			my $curr_int_ = $ms1scan_mz_list->[$index_]->[1];
			
			last if $curr_mz_ > $predict_ios_peaks[$#predict_ios_peaks] + 0.1;
			if($curr_int_ < 3 * $lowest_sig ){
				next;
			}

			my $pred_iso = $predict_ios_peaks[$index_predict_ion_peaks];
			#print "2\n";
			my $measured_ppm_iso = ppm($curr_mz_,$pred_iso);
			#print "$pred_iso,  $curr_mz_,  $measured_ppm_iso @@@\n";
			#print STDERR "check mz: $curr_mz Curr_mz_:$curr_mz_ pred_iso:$pred_iso measured_ppm_iso:$measured_ppm_iso\n";
			if(abs( $measured_ppm_iso ) <= $ms1_ppm_tolerence ){
				# if($control == 1 && $curr_int_ > $last_int){
					# last;
				# }
				
				if($second_last_int > $last_int && $last_int < $curr_int_){
					$index_--;
					last;
				}
				#push @obs_ppm,$measured_ppm_iso;
				push @obs_mz,[$curr_int_,$curr_mz_];
				if($iso_top_intensity < $curr_int_){
					$iso_top_intensity = $curr_int_;
					$iso_top_mz = $curr_mz_;
				}
				# if($2nd_last_int < $last_int && $last_int > $curr_int_){
					# $control = 1;
				# }
				
				push @obs_iso_distribution,$curr_int_;
				$n_of_ios_matched++;
				$index_predict_ion_peaks++;
				#push @obs_iso_distribution,$curr_int_;
				$isotopic_ints += $curr_int_;
				$tmp_trace_mz{$curr_mz_} = 1;
				$second_last_int = $last_int;
				$last_int = $curr_int_;
			}
			
			my $pred_iso_high = $predict_ios_peaks_high[$index_predict_ion_peaks_high];
			#print "3\n";
			my $measured_ppm_iso_high = ppm($curr_mz_,$pred_iso_high);
			#print "$pred_iso,  $curr_mz_,  $measured_ppm_iso @@@\n";
			if(abs( $measured_ppm_iso_high ) <= $ms1_ppm_tolerence ){
				#print "Curr_mz:$curr_mz_, candi: $pred_iso\n";
				#push @obs_ppm,$measured_ppm_iso;
				#push @obs_mz,$curr_mz_;
				$n_of_ios_matched_high++;
				$index_predict_ion_peaks_high++;
				#push @obs_iso_distribution,$curr_int_;
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
		  
		
		#$position,
		# => 
		#$unimod_mono_mass,
		#$site_,
		#$mod_name_,
		#$NL
		#my %sites_candidate = map{$_->[1],[$_->[0],$_->[2],$_->[3],$_->[4]]}@{$site}; # 
		# my %sites_candidate;
		# foreach my $i (@{$site}){
			# push @{$sites_candidate{$i->[1]}},[$i->[0],$i->[2],$i->[3],$i->[4]];
		# }
		#print Dumper \%sites_candidate;
		
		my ($matched_mzs_type,$matched_int,$matched_ppms,$match_mzs_topN5);
		my ($match_type,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$ppm_sd_mod_,$matched_int_reported_,$matched_coverage_,$match_mzs_topN5,$matched_mod_type_,$mod_matched_coverage_);
		#my $theoretical_mz_table_mod;
		#my ($mod_matched_res);
		
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
		#my ($PTMSiteScore_,$PTMSiteProb_,$position_,$unimod_mass_,$NL_,$mod_name_,$mod_AA);
		#foreach my $unimod_name_in (keys %{$PTMsiteScore_modification_Res}){
		my ($PTMSiteScore_,$PTMSiteProb_,$position_,$unimod_mass_,$NL_,$mod_name_,$mod_AA) = @{$PTMsiteScore_modification_Res};
		
		# my $x;
		# ($match_type,$x,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$matched_int_reported_,$matched_coverage_,$match_mzs_topN5,$matched_mod_type_,$mod_matched_coverage_,$position,$mod_name,$unimod_mass) = @{$mod_matched_res};
		
		
		#foreach my $site_position (keys %sites_candidate){
			#print Dumper $pep,$assumed_charge,$mod_mass,$site_mod_test,$sites_candidate{$site_mod_test};
			#foreach my $m (@{$sites_candidate{$site_position}}){
		#my ($unimod_mono_mass,$site_,$mod_name_,$NL) = @{$m};
		my $theoretical_mz_table_mod = peptide::fragmentation_mod($pep,$assumed_charge,$unimod_mass_,$position_,$mod_name_,$NL_,$pepLen);

		#print STDERR "--check6 mz: $curr_mz $site_:\n";
		#print Dumper ($theoretical_mz_table_mod,$peaks,$ms1profile->{$iso_top_mz},$ms2profile);
		($matched_mzs_type,$matched_int,$matched_ppms,$match_mzs_topN5) = ms2spectra_match_mod($theoretical_mz_table_mod,$peaks,$ms1profile->{$iso_top_mz},$ms2profile);
		#print Dumper ($matched_mzs_type,$matched_int,$matched_ppms,$match_mzs_topN5);
		# print STDERR "--check4.2 mz: $curr_mz\n";
		if(! defined $matched_mzs_type ){
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}
			$index++;
			next;
		}
		($match_type,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$matched_int_reported_,$matched_coverage_,$match_mzs_topN5,$matched_mod_type_,$mod_matched_coverage_,$ppm_sd_mod_)= tag_mod($matched_mzs_type,$matched_int,$matched_ppms);
		#print Dumper ($match_type,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$matched_int_reported_,$matched_coverage_,$match_mzs_topN5,$matched_mod_type_,$mod_matched_coverage_);
		#print STDERR Dumper $matched_mzs_type,$matched_int,$matched_ppms,$match_mzs_topN5,$site_;
		# print STDERR "--check5 mz: $curr_mz\n";
		#print Dumper  $match_type,$matched_int_mod; 
		my ($y,$b) =  ($match_type->{"Y"}, $match_type->{"B"}); 
		# print STDERR Dumper "match_mzs_topN5",$match_mzs_topN5;
		#print STDERR Dumper $theoretical_mz_table_mod,$matched_mzs_type;
		# map{print STDERR Dumper  $_, normalization_int($ms2profile->{$_}) }@{$match_mzs_topN5};
		# print STDERR Dumper normalization_int($ms1profile->{$curr_mz});
		
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
		#$ave_corr = ave([@corrlist]);
		$median_corr = median([@corrlist]);
		# print STDERR "median_corr: $median_corr; median_corr_keep: $median_corr_keep  @corrlist\n";
		# print STDERR "--check7 mz: $curr_mz\n";
		if($median_corr < 0.95  || scalar @corrlist < 4 ){
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}
			$index++;
			next;
		}
		
		
		# print STDERR "comp: $curr_mz $cal_premz $mod_mass $position $mod_name \n";
		# print STDERR "comp: 
							# matched_int: $matched_int_with_mod_ > $matched_int_with_mod  ... 
							# matched_coverage: $matched_coverage_ >= $match_coverage
							# isotopic_ints: $isotopic_ints > $matched_isotops_int
							# median_corr: $median_corr > $median_corr_keep
							# mod_matched_coverage: $mod_matched_coverage = $mod_matched_coverage_
							# matched_mod_type_num: $#$matched_mod_type + 1 >=  $matched_mod_type_num
							# $y  > $len_of_matched_y_with_mod [$len_of_matched_y_without_mod]  
							# $b  > $len_of_matched_b_with_mod [$len_of_matched_b_without_mod]   \n";
	
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
			#$#$matched_mod_type + 1 >=  $matched_mod_type_num #&& 
			($median_corr > $median_corr_keep && $isotopic_ints > $matched_isotops_int)
		){
			($PTMSiteScore,$PTMSiteProb,$position,$unimod_mass,$NL,$mod_name) = ($PTMSiteScore_,$PTMSiteProb_,$position_,$unimod_mass_,$NL_,$mod_name_."[".$mod_AA."]");
			#$len_of_matched_by_with_mod = $keep_max_len_of_matched_by;
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
			
			
			
			#$n_of_ios_matched = $n_of_ios_matched > $iso_require ? $iso_require : $n_of_ios_matched; # in case two more peak within 0.01Da, use '>=', rather '=='
			#$tag_iso = $tag_iso < $n_of_ios_matched +1 ? $n_of_ios_matched + 1 : $tag_iso;
			
			$mod_Iso_tag = $n_of_ios_matched +1;
			#$matched_isotops_int = $isotopic_ints;
			$mod_exp_premz = $curr_mz;
			$mod_exp_int = $curr_int;
			$mod_ppm_obs = 0;
			#$mod_mass = $mod_mass;
			$massdiff = $mod_mass;
			#$mod_pos = $keep_mod_pos;
			#print STDERR join "\t","---",$mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$mod_pos,$len_of_matched_by_without_mod,$len_of_matched_by_with_mod,"\n";
			$index = $index_;
		}else{
			foreach my $i (keys %tmp_trace_mz){
				$trace_mz{$i} = 1;
				# print STDERR $curr_mz,"-", $i,"\n";
			}
		}
		$index++;
	} 
	#exit;
	my $mod_peptide;
	if($mod_Iso_tag >= 3 ){
		#printf ("\033[31m%s\t%s\t\033[32m%s\033[0m\n", $filename,"ppm_ave: $ppm_ave","ppm_sd: $ppm_sd");
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
		#exit;print Dumper ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$match_type_report->{'Y'},$match_type_report->{'B'},$matched_int_sd_mod,$ppm_ave_mod,$matched_int_reported,$match_coverage_with_mod,$matched_mod_type,$matched_int_with_mod);
		return ($mod_Iso_tag,$mod_exp_premz,$mod_exp_int,$mod_ppm_obs,$massdiff,$match_type_report->{'Y'},$match_type_report->{'B'},$matched_int_sd_mod,$ppm_ave_mod,$matched_int_reported,$match_coverage_with_mod,$matched_mod_type,$matched_int_with_mod,$pre_SN,$position,$mod_name,$unimod_mass,$mod_peptide,$PTMSiteScore,$PTMSiteProb,$new_calmz,$ppm_sd_mod);	
	}else{
		# print STDERR "\n---FAILED\n";
		#if($cal_premz > $ms1_start_mz && $cal_premz < $ms1_end_mz){
		return (0,$cal_premz);
		#}
		
	}

#}
	#exit;$matched_int_without_mod,$matched_int_sd,$m
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
	# if(abs($mass) < 0.001 && $model_tag eq "re_loc"){
		# $mod_name = [[0,undef,undef,"not-modificated",undef]];
	# }else{
	foreach my $unimod_mono_mass (@unimod_mass_list){
		if($unimod_mono_mass < $mass - 0.1){
			next;
		}elsif($unimod_mono_mass > $mass + 0.1){
			last;
		}else{#print "6\n";
			my $ppm_ = ppm($cal_premz+$unimod_mono_mass/$assumed_charge,$cal_premz+$mass/$assumed_charge);
			# print STDERR join "\t","unimod_site_test:",$unimod_mono_mass,$mass,$ppm_,"\n";
			
			#print STDERR Dumper \%fixed_variable;
			my $tag;
			foreach my $fixed_mod_aa (keys %fixed_variable){
				# print STDERR Dumper $fixed_mod_aa;
				# print STDERR Dumper "abs($fixed_variable{$fixed_mod_aa} + $unimod_mono_mass )";
				
				if(abs($fixed_variable{$fixed_mod_aa} + $unimod_mono_mass ) < 0.01){
					while($peptide =~ /$fixed_mod_aa/g){
						my $pos_fixed_aa = pos($peptide);
						push @{$mod_name},[0,$pos_fixed_aa,$fixed_mod_aa,"no-fixed_mod",undef];
						$tag = 1;
					}
				}
			}
			#print STDERR Dumper "###",$mod_name;
			if($mod_name && $tag){next;}
			# [
            # '79.966331',
            # 2,
            # 'H',
            # 'Phospho',
            # undef
			# ]
			
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
										#push @{$mod_name},$unimods->{$unimod_mono_mass}->{$site_};
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
	#print Dumper $mod_name,$mass,$peptide,$peptide_prev_aa,$peptide_next_aa;
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

sub ms2spectra_match_mod{
	my ($pred,$obs,$ms1profile_i,$ms2profile_i) = (@_); # predicted, observed mz
	#my $ms2_tolence = $ms2_ppm_tolerence;
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	#my @pred_list = sort{$a <=> $b}grep{!$uniq{$_}++}(values %{$pred});
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	#print Dumper ($len1,$len2);
	#my $matched_mzs;
	my $matched_mzs_type;
	my $matched_mzs;
	my $match_mzs_topN5;
	my $ppms_frag;
	my $matched_int;
	my @matched_ints;
	my $matched_type;
	# my $not_matched_tic;
	my %checked_obs_mz;
	
	my ($s,$e) = boundary_index(normalization_int($ms1profile_i)); 
		# foreach my $ion_mz (@{$match_mzs_topN5}){
			# next if ! defined $ms2profile->{$ion_mz};
			# my $n_non_zero;
			# map{if($_){$n_non_zero++}}@{$ms2profile->{$ion_mz}};
			# if ($n_non_zero <= 2 || !( $ms2profile->{$ion_mz}->[1] && $ms2profile->{$ion_mz}->[2] && $ms2profile->{$ion_mz}->[3])){
				# next;
			# }
			# print STDERR Dumper (normalization_int($ms2profile->{$ion_mz}),normalization_int($ms1profile->{$iso_top_mz}),$s,$e);
			# push @corrlist, pearson_corr([@{$ms2profile->{$ion_mz}}[$s..$e]],[@{$ms1profile->{$iso_top_mz}}[$s..$e]]);

		# }	

		
	while($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
		
		#print join "\t",($pred_list[$pred_index1],@{$obs->[$obs_index2]}),"\n";
		if(defined $checked_obs_mz{$obs_index2}){
			$obs_index2++;
		}
		if($pred_index1 > $pred_len1){
			while($obs_index2 <= $obs_len2){
				#push @{$not_matched_mzs},$obs->[$obs_index2] if ! $checked_obs_mz{$obs_index2};
				#$not_matched_tic += $obs->[$obs_index2]->[1];
				$obs_index2++;
			}
			last;
		}
		if($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
			my ($pred_mz1,$obs_mz2,$obs_int2) = ($pred_list[$pred_index1],@{$obs->[$obs_index2]});
			#print "7\n";
			my $measured_ppm = ppm($obs_mz2,$pred_mz1);
			# print STDERR join "\t",$obs_mz2,$pred_mz1,"$measured_ppm\n";
			if(abs($measured_ppm) <= $ms2_ppm_tolerence ){
				#push @{$matched_mzs},$obs->[$obs_index2];
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				my $corr_i = pearson_corr([@{$ms2profile_i->{$obs_mz2}}[$s..$e]],[@{$ms1profile_i}[$s..$e]]);
				#print "====\t", join "\t",$corr_i,$obs_mz2,$obs_int2,$predicted_mz_ion_type,"\n";
				#if($predicted_mz_ion_type =~ /mod/){
					
				if($corr_i > 0.8){
					push @{$matched_mzs},[$obs_mz2,$obs_int2];
					#push @{$matched_type},$predicted_mz_ion_type;
				
					push @{$matched_mzs_type},$predicted_mz_ion_type;
					$ppms_frag->{$predicted_mz_ion_type} = $measured_ppm;
					$matched_int->{$predicted_mz_ion_type} = [$obs_mz2,$obs_int2];
				}
				#push @matched_ints, log($obs_int2);
				$checked_obs_mz{$obs_index2} = 1;
				$pred_index1++;
				$obs_index2++;
			}else{
				if($pred_mz1 < $obs_mz2){
					$pred_index1++;
				}else{
					#print $obs->[$obs_index2],"\n";
					#print Dumper $obs->[$obs_index2];
					#push @{$not_matched_mzs},$obs->[$obs_index2] if ! $checked_obs_mz{$obs_index2};
					#$not_matched_tic += $obs_int2;
					$obs_index2++;
				}
			}
		}
	}
	
	#print STDERR Dumper "mod",$matched_mzs;
	my @sorted = map{$_->[0]}sort{$b->[1] <=> $a->[1]}@{$matched_mzs};
	$match_mzs_topN5 = [@sorted[0..4]];
	#return ($matched_mzs,$matched_mzs_type,$not_matched_tic);
	# my $ppm_ave = ave(\@ppms_frag);
	# my $m_int_ave;
	# my $m_int_sd;
	# if(scalar @matched_ints > 1){
		# $m_int_ave = ave(\@matched_ints);
		# $m_int_sd = stdev(\@matched_ints,$m_int_ave);
	# }
	#my $ppm_sd; 
	#$ppm_sd = stdev(\@ppms_frag,$ppm_ave) if $#ppms_frag >= 1;
	#print STDERR "matched_type: ",join "\t",@{$matched_type},"\n" ;
	
	return ($matched_mzs_type,$matched_int,$ppms_frag,$match_mzs_topN5);
}



sub ms2spectra_match{
	my ($pred,$obs) = (@_); # predicted, observed mz
	#print Dumper ($pred,$obs);
	#my $ms2_tolence = $ppm_tolerence;
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	#my @pred_list = sort{$a <=> $b}grep{!$uniq{$_}++}(values %{$pred});
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	#print Dumper ($len1,$len2);
	my $matched_mzs;
	my $matched_mzs_type;
	my @ppms_frag;
	my $matched_int;
	my @matched_ints;
	# my $not_matched_mzs;
	# my $not_matched_tic;
	my %checked_obs_mz;
	while($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
		
		#print join "\t",($pred_list[$pred_index1],@{$obs->[$obs_index2]}),"\n";
		if(defined $checked_obs_mz{$obs_index2}){
			$obs_index2++;
		}
		if($pred_index1 > $pred_len1){
			while($obs_index2 <= $obs_len2){
				#push @{$not_matched_mzs},$obs->[$obs_index2] if ! $checked_obs_mz{$obs_index2};
				#$not_matched_tic += $obs->[$obs_index2]->[1];
				$obs_index2++;
			}
			last;
		}
		if($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
			my ($pred_mz1,$obs_mz2,$obs_int2) = ($pred_list[$pred_index1],@{$obs->[$obs_index2]});
			
			my $measured_ppm = ppm($obs_mz2,$pred_mz1);
			# print STDERR join "\t",$obs_mz2,$pred_mz1,"$measured_ppm\n";
			if(abs($measured_ppm) <= $ms2_ppm_tolerence){
				#push @{$matched_mzs},$obs->[$obs_index2];
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				push @{$matched_mzs},$obs_mz2;
				push @{$matched_mzs_type},$predicted_mz_ion_type;
				push @ppms_frag,$measured_ppm if $predicted_mz_ion_type =~ /y\+|b\+/;
				$matched_int += $obs_int2 if $predicted_mz_ion_type =~ /y\+|b\+/;
				push @matched_ints, log($obs_int2) if $predicted_mz_ion_type =~ /y\+|b\+/;
				$checked_obs_mz{$obs_index2} = 1;
				$pred_index1++;
				$obs_index2++;
			}else{
				if($pred_mz1 < $obs_mz2){
					$pred_index1++;
				}else{
					#print $obs->[$obs_index2],"\n";
					#print Dumper $obs->[$obs_index2];
					#push @{$not_matched_mzs},$obs->[$obs_index2] if ! $checked_obs_mz{$obs_index2};
					#$not_matched_tic += $obs_int2;
					$obs_index2++;
				}
			}
		}
	}
	#return ($matched_mzs,$matched_mzs_type,$not_matched_tic);
	my $ppm_md; 
	my $ppm_sd;
	
	
	my $m_int_ave;
	my $m_int_sd;
	if(scalar @matched_ints > 1){
		$m_int_ave = ave(\@matched_ints);
		$m_int_sd = stdev(\@matched_ints,$m_int_ave);
		$ppm_md = median(\@ppms_frag);
		my $ppm_ave = ave(\@ppms_frag);
		$ppm_sd = stdev(\@ppms_frag,$ppm_ave);
	}
	#my $ppm_sd; 
	#$ppm_sd = stdev(\@ppms_frag,$ppm_ave) if $#ppms_frag >= 1;
	#print Dumper ($matched_mzs_type,$matched_int,$m_int_sd || "NA",$ppm_md,$matched_mzs,$ppm_sd);
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
		#print $ion;
		#$iontype{$pro=~/REV_/?"Decoy":"Target"}{$fragment_type}++;
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
	#print   Dumper "match",\%uniq;
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
		# if ($fragment_type =~ /mod/){
			# $matched_ints_mod += $uniq{$i}->[0];
		# }
		
	}
	#print Dumper \%mod_types;#exit;
	my %tag_len;
	foreach my $type (keys %ion_list){
		my @n = sort{$a <=> $b}@{$ion_list{$type}};
		#print Dumper $type,\@n;
		my $tag=1;
		my $max_len;
		foreach my $i(0..$#n-1){
			if($n[$i]+1 == $n[$i+1]){
				$tag++;
				#print $n[$i],"\t",$tag,"\n"
			}elsif($n[$i]+2 == $n[$i+1]){
				$tag++;
				$i++;
			}else{
				#$max_len = $tag if $max_len < $tag;
				$max_len += $tag if $tag >= 3;
				$tag=1;
			}
		}
		$max_len = $tag if $max_len < $tag;
		#print Dumper $type, \@tags;
		$tag_len{$type} = $max_len; 
	}
	#print Dumper "ionlist",\%ion_list;
	#print Dumper \%mod_types,\%ion_list;exit;
	
	
	my @n = sort{$a <=> $b}(keys %mod_types);
	#print Dumper $type,\@n;
	my $tag=1;
	my $max_len;
	foreach my $i(0..$#n-1){
		if($n[$i]+1 == $n[$i+1]){
			$tag++;
			#print $n[$i],"\t",$tag,"\n"
		}elsif($n[$i]+2 == $n[$i+1]){
			$tag++;
			$i++;
		}else{
			#$max_len = $tag if $max_len < $tag;
			$max_len += $tag if $tag >= 3;
			$tag=1;
		}
	}
	$max_len = $tag if $max_len < $tag;
	#print Dumper $type, \@tags;
	my $mod_coverage_n = $max_len;
	
	
	my @n = sort{$a <=> $b}(keys %ion_for_coverage);
	#print Dumper $type,\@n;
	my $tag=1;
	my $max_len;
	foreach my $i(0..$#n-1){
		if($n[$i]+1 == $n[$i+1]){
			$tag++;
			#print $n[$i],"\t",$tag,"\n"
		}elsif($n[$i]+2 == $n[$i+1]){
			$tag++;
			$i++;
		}else{
			#$max_len = $tag if $max_len < $tag;
			$max_len += $tag if $tag >= 3;
			$tag=1;
		}
	}
	$max_len = $tag if $max_len < $tag;
	#print Dumper $type, \@tags;
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
	
	
	#print Dumper \@matched_mods;exit;
	
	my @sorted = map{$_->[0]}sort{$b->[1] <=> $a->[1]}@matched_mods;
	$match_mods_topN5 = [@sorted[0..4]];
	# print STDERR "match_mods_type: ",join "\t",@{$matched_mod_type},"\n" if defined $matched_mod_type;
	# print STDERR "match_topN5_type: ",join "\t",(map{if($_){$_->[2]}}(sort{$b->[1] <=> $a->[1]}@matched_mods)[0..4]),"\n" if scalar @matched_mods;
	# print STDERR "mod_coverage_n: $mod_coverage_n\n";
	#my $coverage_n = $tag_len{'Y'} + $tag_len{'B'};
	#$match_type,$matched_int_with_mod_, $matched_int_sd_mod_,$ppm_ave_mod_,$matched_int_reported_
	return (\%tag_len,$matched_ints_mod,$m_int_sd || "NA",$ppm_md,$matched_ints_mod_reported,$coverage_n,$match_mods_topN5,$matched_mod_type,$mod_coverage_n,$ppm_sd);
}



sub pearson_corr{
    my ($ref_a, $ref_b) = @_;
	#print Dumper ($ref_a, $ref_b);
    #my @x = map{$_ ? log($_)/log(10) : 0}@{$ref_a};
    my @x = @{$ref_a};
    #my @y = map{$_ ? log($_)/log(10) : 0}@{$ref_b};
    my @y = @{$ref_b};
	my $mean_x = sum(@x)/($#x+1);
	my $mean_y = sum(@y)/($#y+1);
	#my $len = $#x + 1;
	my $correlation;
	#if($#x==$#y){
		 my $N = $#x;
		 my $sum_sq_x = 0;
		 my $sum_sq_y = 0;
		 my $sum_coproduct = 0;
		 #my $mean_x = $x[1];
		 #my $mean_y = $y[1];
		 for(my $i=0;$i<=$N;$i++){
				 # my $sweep = ($i - 1.0) / $i;
				  my $delta_x = $x[$i] - $mean_x;
				  my $delta_y = $y[$i] - $mean_y;
				  $sum_sq_x += $delta_x * $delta_x ;
				  $sum_sq_y += $delta_y * $delta_y ;
				  $sum_coproduct += $delta_x * $delta_y;
				  #$mean_x += $delta_x / $i;
				  #$mean_y += $delta_y / $i;
		 }
		 my $pop_sd_x = sqrt( $sum_sq_x );
		 my $pop_sd_y = sqrt( $sum_sq_y );
		 my $cov_x_y = $sum_coproduct;
		 $correlation = $cov_x_y / ($pop_sd_x * $pop_sd_y) if $pop_sd_x * $pop_sd_y;
		 #print OUT "$name1 $name2: $correlation\n";
	#}
	#return ($correlation, [$mean_x, $sum_sq_x, $mean_x, $pop_sd_x], [$mean_y, $sum_sq_y, $mean_y, $pop_sd_y]);
	#print Dumper ($correlation);#exit;
	return $correlation;
}


sub deisotope{
	my $peaks = shift;
	
	my @deisotope_peaks;
	# charge_range {1..2};
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
				#print join "\t",$mz_next,$ind_next,$mz_z2,$ppm_next,"-2-2\n";
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
				#print join "\t",$mz_next,$ind_next,$mz_z1,$ppm_next,"-1-1\n";
				if(abs($ppm_next) <= $ms2_ppm_tolerence && $ind_next < $int_last * 1.5){
					($mz_last,$int_last) = ($mz_next,$ind_next);
					$rm_inds{$ind_} = 1;
				}
			}

		push @deisotope_peaks, $peaks->[$ind];
		# ($mz_last,$int_last) = ($mz_next,$ind_next);
	}
	return [@deisotope_peaks];
}

sub peakpick{
	my $peaks = shift;
	#print Dumper $peaks;
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
		#print Dumper $window{$w}, [@sorted_int];
		foreach my $i (0..$max_peakdepth -1 ){
			#print join "---",$i,$sorted_int[$i]->[0],"\n";
			if($sorted_int[$i]->[0]){
				$sorted_peaks{$sorted_int[$i]->[0]} = [$i+1,$w];
				map{$extracted_ions_num{$_}++}$i+1..$max_peakdepth;
			}
		}
		#$sorted_peaks{$w} = [@sorted_int[0..9]];
	}
	#print Dumper \%extracted_ions_num,\%sorted_peaks;exit;
	return (\%sorted_peaks,\%extracted_ions_num);
}


sub ms2spectra_Site{
	my ($pred,$obs,$top8_mz) = (@_); # predicted, observed mz
	#my $ms2_tolence = $ms2_ppm_tolerence;
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	#my @pred_list = sort{$a <=> $b}grep{!$uniq{$_}++}(values %{$pred});
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	#print Dumper ($len1,$len2);
	#my $matched_mzs;
	my $matched_mzs_type;
	my $matched_mzs;
	my @ppms_frag;
	my $matched_int;
	my @matched_ints;
	my $therotical_peak_n;
	# my $not_matched_mzs;
	# my $not_matched_tic;
	my %checked_obs_mz;
	while($pred_index1 <= $pred_len1 && $obs_index2 <= $obs_len2){
		
		#print join "\t",($pred_list[$pred_index1],@{$obs->[$obs_index2]}),"\n";
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
			# print STDERR join "\t",$obs_mz2,$pred_mz1,"$measured_ppm\n";
			if(abs($measured_ppm) <= $ms2_ppm_tolerence ){
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				my $top_n = $top8_mz->{$obs_mz2}->[0];
				#my $window = $top10_mz->{$obs_mz2}>[1];
				#print STDERR join "##",$top_n,$predicted_mz_ion_type,"\n";
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
	# print STDERR "PeakdepthOptimizationScore\n";
	# print STDERR Dumper @_[1..6];
	my $matched_info = shift;
	my $w = shift;
	my $n = shift;
	my $site_pos = shift;
	my $Score = shift;
	my $Npeak = shift;
	my $_score = shift;
	# $max_score_of_peakdepth;
	my $P;
	#my %Score;
	#print STDERR Dumper $matched_info;
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
				#$Score->{$i}->{$site_pos} = 0;
				$s = 0;
				$Score->{$peakdepth}->{$site_pos} = $s;
				next;
			}
			my $k = scalar @{$matched_info->{$peakdepth}};
			#my $k = matched_k($matched_info->{$i});
			my $p = $Npeak->{$peakdepth} * 0.02 / $w;
			# print STDERR Dumper "PeakdepthOptimizationScore", ($n,$k,$p);
			if($p == 0){
				#$P = 0;
				$s = 0;
				#$Score->{$i}->{$site_pos} = 0;
			}else{	
				$P = Cumulative_binominal($n,$k,$p);
				#$P = binominal($n,$k,$p);
				if($P == 0){
					$s = 0;
					#$Score->{$i}->{$site_pos} = 0;
				}else{
					$s = -10 * log($P)/log(10);
					#$Score->{$i}->{$site_pos} = -10 * log($P)/log(10);
				}
				#print -10 * log($P)/log(10),"\n";
			}
			#print Dumper "PeakdepthOptimizationScore", $P;
		}else{
			$s = 0;
			#$Score->{$i}->{$site_pos} = 0;
		}
		# print STDERR "score = $s\n";
		$Score->{$peakdepth}->{$site_pos} = $s;
		$_score->{$site_pos}+= $s;
		
		# if($max_score_of_peakdepth->[0] < $Score->{$i}->{$site_pos}){
			# $max_score_of_peakdepth = [$Score->{$i}->{$site_pos}, $site_pos,$i];
		# }
	}
	#print STDERR Dumper ($Score,$_score);
	return ($Score,$_score);
}

# sub matched_k {
	# my $matched_ion_types = shift;
	# # 9.b++[mod]
	# my %slim_ion;
	# foreach my $ion (@{$matched_ion_types}){
		# $ion =~ /(^\d+\.[by]).+/;
		# $slim_ion{$1} = 1;
	# }
	# my $k = scalar (keys %slim_ion);
	# return $k;
# }

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
	#my $biprob = $n == $k ?  ( $p ** $k ) * ((1 - $p) ** ($n - $k)) : (factorial($n) / (factorial($k) * factorial($n - $k))) * ( $p ** $k ) * ((1 - $p) ** ($n - $k));
	my $biprob = $n == $k ?  ( $p ** $k ) * ((1 - $p) ** ($n - $k)) : comb_p($n,$k) * ( $p ** $k ) * ((1 - $p) ** ($n - $k));
	if($biprob eq "Inf"){
		$biprob = 0;
	}
	#print $biprob,"====44\n";
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
	#my $max_score_of_peakdepth;
	my @SitePos;
	my %modfied_sites;
	while($strip_pep =~ m{[$Site]}g){
		my $pos = pos($strip_pep);
		push @SitePos,$pos;
		$modfied_sites{$pos} = $&;
	}
	#print Dumper @SitePos;
	
	my @CombSitePos = @{sitecombination([@SitePos],$site_num)};
	#print STDERR Dumper \@CombSitePos;
	#print Dumper @CombSitePos;exit;
	my %start_points_adj;
	foreach my $c(@CombSitePos){
		#print join "#",($peptide_mod,$charge,(join "#",@$c),length($peptide_mod)),"\n";
		my $putative_frags = peptide::fragmentation_phosSite($peptide_mod,$charge,$c,length($strip_pep));
		#print Dumper $putative_frags;
		my ($matched_frags,$therotical_peak_n) = ms2spectra_Site($putative_frags,$deisotope_peaks,$TOP6_peaks);
		#print STDERR Dumper "positon: ",$c,$matched_frags,$therotical_peak_n;
		#$score = PeakdepthOptimizationScore($matched_frags,$deisotope_peaks->[$#$deisotope_peaks]->[0],$therotical_peak_n,(join "#",@$c),$score,$extracted_ions_n);
		($score,$score_for_peakdepth) = PeakdepthOptimizationScore($matched_frags,$deisotope_peaks->[$#$deisotope_peaks]->[0],$therotical_peak_n,(join "#",@$c),$score,$extracted_ions_n,$score_for_peakdepth);
		#print STDERR Dumper $matched_frags,$deisotope_peaks->[$#$deisotope_peaks]->[0],$therotical_peak_n,(join "#",@$c),$score,$extracted_ions_n,$score_for_peakdepth;exit;
		#print Dumper $score,$score_for_peakdepth,"\n";
		#push @{$start_points_adj{$max_score_of_peakdepth->[1]}},$c;
	}
	#exit;
	if(! defined $score){
		return ($peptide_mod,"-","-","-",0);
	}else{
		my ($peakdepth_of_max_score,$pos_of_peakdepth) = max_diff_score($score);
		
		#my $max_score = max(keys %{$score_for_peakdepth});
		#my $peakdepth_of_max_score = min(keys %{$score_for_peakdepth->{$max_score}});
		#my ($pos_of_peakdepth) = (keys %{$score_for_peakdepth->{$max_score}->{$peakdepth_of_max_score}});
			
		#print STDERR "peakdepth_of_max_score",$peakdepth_of_max_score,"\n";
		#print STDERR Dumper $score,"\n";
		my $modification_sitepos = $pos_of_peakdepth;
		my %PTMSiteScore_res;
		my $sumScore;
		#print STDERR Dumper \@CombSitePos;
		
		foreach my $c(@CombSitePos){
			my $pos = join "#",@{$c};
			#print STDERR Dumper "pos",$pos;
			foreach my $site (@{$c}){
				$PTMSiteScore_res{$site} += 1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;
				print STDERR "$score->{$peakdepth_of_max_score}->{$pos} > 0;\n"
				#print join "\t",$peptide_mod,$site,$score->{$peakdepth}->{$pos},$pos,$peakdepth,"\n";
			}
			$sumScore +=  1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;
		}
		#print STDERR "sumScore: $sumScore\n";
		#print STDERR Dumper "PTMSiteScore_res:", \%PTMSiteScore_res;
		#exit;
		my $PTM_peptide = $peptide_mod;
		my $PTM_score_loc;
		my $ptm;
		my $max_prob;
		my $localized_site;
		my $localized_or_not;
		if($sumScore){
			print STDERR "PTMsiteScore\t",$peptide_mod,":\t", $peakdepth_of_max_score, "\t",$score->{$peakdepth_of_max_score}->{$modification_sitepos},"\t";
			print STDERR join ";",map{sprintf "%d(%.2f%%)",$_,100*$PTMSiteScore_res{$_}/$sumScore}@SitePos;
			
			my @siteprob;
			map{
				#if($max_prob < 100*$PTMSiteScore_res{$_}/$sumScore){
					push @siteprob,[$_,100*$PTMSiteScore_res{$_}/$sumScore];
					#$localized_site = $_;
				#}
			}@SitePos;
			my @sorted_siteprob = sort{$b->[1] <=> $a->[1]}@siteprob;
			$localized_site = join ";",map{$_->[0] . "(" . $_->[1] . ")"}@sorted_siteprob[0..$site_num-1];
			$ptm = join ";",map{sprintf "%d(%.2f%%)",$_,100*$PTMSiteScore_res{$_}/$sumScore}@SitePos;
			map{$PTM_peptide = peptide_phos_add($PTM_peptide, $_->[0])}@sorted_siteprob[0..$site_num-1];	
			$localized_or_not = $sorted_siteprob[$site_num-1]->[1] >= 75 ? 1 : 0;
			$PTM_score_loc = $score->{$peakdepth_of_max_score}->{$modification_sitepos};
			print STDERR "\n";
			# my $ind;
			
			# my @res_modificaiton_sites = split /#/,$modification_sitepos;
			# foreach my $site (@res_modificaiton_sites){
				# #$ind++;
				# my $ptmscore = sprintf ("%.2f",$PTMSiteScore_res{$site}/$sumScore);
				# #if( $ptmscore >= 0.75){
				# $PTM_peptide = peptide_phos_add($PTM_peptide, $site, $ptmscore);	
				# #}
				# push @{$ptm},$modfied_sites{$site}.$site."($ptmscore)";
			# }
			
			#return $PTM_peptide;
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
	#print STDERR Dumper $score_;
	foreach my $depth(keys %{$score_}){
		my @sorted_pos = (sort{$score_->{$depth}->{$b} <=> $score_->{$depth}->{$a}}keys %{$score_->{$depth}});
		push @{$max_score_depth{$score_->{$depth}->{$sorted_pos[0]}}},$depth;
		if($max_score < $score_->{$depth}->{$sorted_pos[0]}){
			$max_score = $score_->{$depth}->{$sorted_pos[0]};
		}
		my $diff; 
		# if($score_->{$depth}->{$sorted_pos[1]} == 0 && (scalar @sorted_pos) > 1){
			# next;
		# }
		if((scalar @sorted_pos) == 1){
			$diff = $score_->{$depth}->{$sorted_pos[0]};
			#next;
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
	#print STDERR Dumper \%score_diff,\%diff_depth,\%depth_pos;
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
		#print join "-",$a,$n,$t,$mod_peptide,"\n";
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
	#print STDERR \%mod_aa;
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
			#print STDERR "re_loc", "$mod_aa{$mod_AA} - $mono_aa{$mod_AA} + $unimod_mass", "\t",$new_mod_mass,"\n";
			my $new_site = unimod_site_test($new_mod_mass,$charge,$cal_premz,$strip_pep,$ms1_ppm_tolerence,$peptide_prev_aa,$peptide_next_aa);
			#print STDERR Dumper $new_site;
			# if($new_site->[0]->[3] eq "not-modificated"){
				# my $new_uniq_key = join ":",$pos,0,undef,$mod_AA,"re_loc";
				# #push @{$CombSitePos->{$uniq_key}}, [$pos,$unimod_mass,$NL,$unimod_name,$mod_AA];
				# $Site_uniq->{$new_uniq_key}{"not-modificated"} = 1;
				# next;
			# }
			foreach my $new_i (@{$new_site}){
				my $new_pos = $new_i->[1];
				my $new_unimod_mass = $new_i->[0];
				my $new_NL = $new_i->[-1];
				my $new_unimod_name = $new_i->[-2];
				my $new_mod_AA = $new_i->[2];
				if($new_mod_AA ne $mod_AA){next;}
				my $new_uniq_key = join ":",$new_pos,$new_unimod_mass,$new_NL,$new_mod_AA,"re_loc";
				#push @{$CombSitePos->{$uniq_key}}, [$pos,$unimod_mass,$NL,$unimod_name,$mod_AA];
				$Site_uniq->{$new_uniq_key}{$new_unimod_name} = 1;
			}
			next;
			#print Dumper $site;
		}
		
		my $uniq_key = join ":",$pos,$unimod_mass,$NL,$mod_AA,"loc";
		#push @{$CombSitePos->{$uniq_key}}, [$pos,$unimod_mass,$NL,$unimod_name,$mod_AA];
		$Site_uniq->{$uniq_key}{$unimod_name} = 1;
		my $aa = $i->[2];
		$modfied_sites{$pos} = $aa;
	}
	#print STDERR Dumper $Site_uniq;
	#my @CombSitePos = [@SitePos]; # current version only support single mod
	foreach my $key (keys %{$Site_uniq}){
		my ($pos,$unimod_mass,$NL,$mod_AA,$loc) = split /:/,$key;
		my $unimod_name = join ";",sort{$a cmp $b}(keys %{$Site_uniq->{$key}});
		push @{$CombSitePos->{$unimod_mass}}, [$pos,$unimod_mass,$NL,$unimod_name,$mod_AA,$loc];
	}
	
	#print Dumper $CombSitePos;#exit;
	my $ptmscore_res;
	
	#print Dumper $CombSitePos;
	foreach my $unimod_mass(keys %{$CombSitePos}){
		#print join "#",($peptide_mod,$charge,(join "#",@$c),length($peptide_mod)),"\n";
		my $score;
		my $score_for_peakdepth;
		my %start_points_adj;
		#my $max_score_of_peakdepth;
		my $pos_cnt  = scalar @{$CombSitePos->{$unimod_mass}};
		#print Dumper $CombSitePos->{$unimod_mass};
		foreach my $i (@{$CombSitePos->{$unimod_mass}}){
			
			my ($site_position,$unimod_mono_mass,$NL,$mod_name_,$aa,$loc) = @{$i};
			# if($pos_cnt > 1 && ($site_position eq "N-" or $site_position eq "C-")){
				# next;
			# }
			#print STDERR join "\t",($site_position,$unimod_mono_mass,$NL,$mod_name_,$aa) ,"\n";
			my $new_peptide_mod;
			if($loc eq "re_loc"){
				
				my $old_peptide_mod = $peptide_mod;
				my ($n,$t,$a);
				while($old_peptide_mod =~ m{([A-Z])(\[([^][]+)\])?}g){
					$a = $1;
					$n += length($2);
					$t = $3;
					my $pos_ = pos($old_peptide_mod)-$n;
					#print join "-",$a,$n,$t,$mod_peptide,"\n";
					if($pos_ == $site_position){
						$new_peptide_mod .= $a;
						
					}else{
						$new_peptide_mod .= $&;
					}
				}
				#print STDERR Dumper $peptide_mod,$new_peptide_mod;
				#$peptide_mod = $new_peptide_mod;
			}else{
				$new_peptide_mod = $peptide_mod;
			}
			my $putative_frags = peptide::fragmentation_mod($new_peptide_mod,$charge,$unimod_mono_mass,$site_position,$mod_name_,$NL,length($strip_pep));
			#my $putative_frags = peptide::fragmentation_phosSite($peptide_mod,$charge,$c,length($peptide_mod),1);
			#print Dumper $putative_frags;
			my ($matched_frags,$therotical_peak_n) = ms2spectra_Site($putative_frags,$deisotope_peaks,$TOP6_peaks);
			# print STDERR Dumper $matched_frags;
			($score,$score_for_peakdepth) = PeakdepthOptimizationScore($matched_frags,$deisotope_peaks->[$#$deisotope_peaks]->[0],$therotical_peak_n,$site_position,$score,$extracted_ions_n,$score_for_peakdepth);
		}
		#print STDERR Dumper "score",$score;
		if(! defined $score){
			$ptmscore_res = undef;
		}else{
			my ($peakdepth_of_max_score,$pos_of_peakdepth) = max_diff_score($score);
			#my $max_score = max(keys %{$score_for_peakdepth});
			#my $peakdepth_of_max_score = min(keys %{$score_for_peakdepth->{$max_score}});
			#my ($pos_of_peakdepth) = (keys %{$score_for_peakdepth->{$max_score}->{$peakdepth_of_max_score}});
			
			#print join "\t",$max_score,$peakdepth_of_max_score,$pos_of_peakdepth,"\n";
			
			#my $peakdepth;
			my $modification_sitepos;
			my $modification_unimod_mono_mass;
			my $modification_NL;
			my $modification_mod_name;
			my $modification_mod_AA;
			my $max_diff_score;
			
			# foreach my $i (@{$CombSitePos->{$unimod_mass}}){
				# if($i->[0] eq $pos_of_peakdepth){
					# ($modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA) = (@{$i});
					# print join "\t",($modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA),"\n";
					# last;
				# }	
			# }
			my %PTMSiteScore_res;
			my $sumScore;
			foreach my $c(@{$CombSitePos->{$unimod_mass}}){
				my $pos = $c->[0];
				$PTMSiteScore_res{$pos} = 1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;

				$sumScore +=  1 / (10 ** (-$score->{$peakdepth_of_max_score}->{$pos}/10)) if $score->{$peakdepth_of_max_score}->{$pos} > 0;
				if($pos eq $pos_of_peakdepth){
					($modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA) = (@{$c});
					#print join "\t",($modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA),"\n";
					#last;
				}	
			}
			#exit;
			my $PTM_peptide;
			my $ptmscore;
			my $PTMSiteScore;
			my $new_calmz;
			if($sumScore){
				# print STDERR $peptide_mod,":\t unimod_mass: $unimod_mass\t","peakdepth_of_max_score: ", $peakdepth_of_max_score, "\t", $score->{$peakdepth_of_max_score}->{$modification_sitepos},"\t";
				# print STDERR join ";",map{sprintf "%d(%.2f%%)",$_->[0],100*$PTMSiteScore_res{$_->[0]}/$sumScore}@{$CombSitePos->{$unimod_mass}};
				# print STDERR "\n";
				my $ind;
				
				#my @res_modificaiton_sites = split /#/,$modification_sitepos;
				#foreach my $site (@res_modificaiton_sites){
					#$ind++;
				$ptmscore = sprintf ("%.2f",$PTMSiteScore_res{$modification_sitepos}/$sumScore);
				if( $ptmscore >= 0.75){
					$PTM_peptide = $peptide_mod;
					$PTMSiteScore = $score->{$peakdepth_of_max_score}->{$modification_sitepos};		
				}
				#push @{$ptm},$modfied_sites{$site}.$site."($ptmscore)";
				#}
				#return $PTM_peptide;
			}
			#print Dumper join "#",$PTMSiteScore,$ptmscore,$modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA,"\n";
			if($ptmscore_res->[0] < $PTMSiteScore){
				$ptmscore_res = [$PTMSiteScore,$ptmscore,$modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name,$modification_mod_AA];
			}elsif($ptmscore_res->[0] == $PTMSiteScore && $ptmscore_res->[2] == $modification_sitepos){
				$ptmscore_res = [$PTMSiteScore,$ptmscore,$modification_sitepos,$modification_unimod_mono_mass,$modification_NL,$modification_mod_name.";".$ptmscore_res->[5],$modification_mod_AA];
			}
			#
		}
		#exit;
		#print Dumper $ptmscore_res;
		
	}
	return ($ptmscore_res);
}
