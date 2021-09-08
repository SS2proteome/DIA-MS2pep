package pseudo_ms2_InParallel;

# charge_script [$file] [$charge_dir] [ms2 ppm] [data points] [min_pearson] [min_res_pearson] [min_peaks_splited]

# 20191031 1. one mgf files, 2. no precharge adjustment, 3. keep the mgf with the basepeak of the highest intensity within 30s.

# 201912010 re-write the processing workflow of data spliting, decide which signal to be centric, then do split 

use strict;
use Data::Dumper;
use POSIX qw(floor);
use List::Util qw(max min sum);
use POSIX qw(floor ceil round);
use File::Basename;
use lib dirname(__FILE__);
use peptide;
my $PI = 4 * atan2(1, 1);
my $massdiff_C12_C13 = 1.0033548378;

use Statistics::KernelEstimation;


$| = 1;



my @ion_arr;
my @low_signals;
my %lowest_intenstiy;
my %basepeak_mz;
my %bp_int;
my %spectra_data;
my %spectra_noise;

my $buffer1_ms2scan_data;
my $buffer2_ms2scan_data;
my $buffer3_ms2scan_data;

my $buffer1_ms2scan_low_signal_data;
my $buffer2_ms2scan_low_signal_data;
my $buffer3_ms2scan_low_signal_data;

my $last1_ms2scan;
my $last2_ms2scan;
my $last1_title_line;
my $last2_title_line;
my $last1_rt_line;
my $last2_rt_line;
my $last1_pepmass_line;
my $last2_pepmass_line;
my $rt_;

my $mz_buffer_for_alignment;

my $ms2scan;
my %DIA_window;
my %isolationlist;
my $number_of_windows;

my %ms2toms1;
my $scan;
my $ms1scan;
my $ms2scan;
my %basepeak_intensity;
my %basepeak_intensity_whole;
my $number_of_scan;
my $scantime;
my $lastscantime;
my $cycletime;
my $cyclecount;
my @basepeak_ms1_hash;
my $basepeak_ms1;
my $currbasepeak_ms1;
my $lastbasepeak_ms1;
my $basepeak_ms1_cnt;
my $premz_target;
my $mslevel;

my @total_ms2scan;

my $title_line;
my $rt_line;
my $pepmass_line;
my $ion_lines;
my $premz;
my $scan_premz_index;
my %ms2scan_of_index;
my %ms1scan_of_index;

my %ms2spectra_data;
my %ms2spectra_TIC;
my %ms1spectra_data;
my %ms2scan_header;
my %stat;
my %stat_basepeak;
my $min_peaks_splited = 20;
my $ms2_low_mz = 140;

my $min_data_points = 3;

my $min_charge = 1;
my $max_charge = 5;
my %filehandle;
my $num_of_data_points;
my $mgf_file;
my $charge_dir;
my $ms2_precison;
my $window_file;

sub demultiplexing {

	$mgf_file = "@_[0]";
	$charge_dir = "@_[1]";
	$ms2_precison = "@_[2]";# ppm;
	$num_of_data_points = "@_[3]";

	$window_file = "@_[4].mzML.DIA_acquisition_window.txt";


	open(FH,$window_file) or die "$!\n $window_file";
	while(<FH>){
		chomp;
		next if $. == 1;
		my @line = split /\t/,$_;
		if($line[-1] >100){$line[-1] = 100;}
		$DIA_window{$line[-2]}{size} = $line[-1];
		$DIA_window{$line[-2]}{lower} = $line[1];
		$DIA_window{$line[-2]}{higher} = $line[2];
		$isolationlist{$line[-1]} = 1;
		$number_of_windows = $line[0];
	}
	close(FH);
	#print Dumper \%DIA_window; exit;




	#  READ MGF FILE
	print "\n\n*****Read mgf_file file: $mgf_file!\n";

	open(FH,"${charge_dir}/$mgf_file") or die "$! ${charge_dir}/$mgf_file \n";




	#my %filehandle_num;
	foreach my $key (keys %isolationlist){
		open($filehandle{$key},">${charge_dir}/w_${key}_$mgf_file") or die "Can't write mgf file! ${charge_dir}/w_${key}_$mgf_file $!\n";
	}



	while(<FH>){
		chomp;
		my $Line = $_;
		if(/^(TITLE=[^.]+\.([^.]+)\.[^.]+)/){
			$title_line = $Line;
			$ms2scan = $2;
			printf  "\r$ms2scan";
			undef @ion_arr;	
			undef @low_signals;	
		}elsif(/^RTINSECONDS=([\d.]+)/){
			$rt_line = $Line;
			$rt_ = $1;
		}elsif (/^PEPMASS=([\d.]+)/){
			$pepmass_line = $Line;
			$premz = $1;

		}elsif(/^\d+/){
			my @line = split /\s+/,$_;
			if ($line[0] > $ms2_low_mz  &&
				($line[0] < $DIA_window{$premz}{lower} - 3  || $line[0] > $DIA_window{$premz}{higher} + 3 ) #&&
				){
				push @ion_arr,[@line];
				if($basepeak_intensity{$ms2scan}->[1] < $line[1]){
					$basepeak_intensity{$ms2scan} = [@line];	
				}
			}
		}elsif(/^END/){
				
			my $filtered_ion_arr;
			if($buffer1_ms2scan_data->{$premz}->[0]){
				$filtered_ion_arr->{$premz} = $buffer1_ms2scan_data->{$premz}
			}
			
			if (defined $filtered_ion_arr->{$premz} ){
					$scan_premz_index->{$premz}++;
					$ms2scan_of_index{$premz}{$scan_premz_index->{$premz}} = $last2_ms2scan->{$premz};		
					$ms2scan_header{$last2_ms2scan->{$premz}} = "BEGIN IONS\n$last2_title_line->{$premz}\n$last2_rt_line->{$premz}\n$last2_pepmass_line->{$premz}\n";
					
					for(1..$num_of_data_points-1){
						$ms2spectra_data{$premz}{$_} = $ms2spectra_data{$premz}{$_+1};
					}
					$ms2spectra_data{$premz}{$num_of_data_points} = $filtered_ion_arr->{$premz};
					
					
					
					if($scan_premz_index->{$premz} >= $num_of_data_points){
						Spectrum_ms2_seed_finder(
							$scan_premz_index->{$premz} - 1,
							$premz
						);
					
					}
			}
			
			#print Dumper "ddd",@ion_arr[0];
			
			$buffer1_ms2scan_data->{$premz} = $buffer2_ms2scan_data->{$premz};
			$buffer2_ms2scan_data->{$premz} = $buffer3_ms2scan_data->{$premz};
			$buffer3_ms2scan_data->{$premz} = [@ion_arr];
			#print Dumper $buffer1_ms2scan_data->[0],$buffer2_ms2scan_data->[0],$buffer3_ms2scan_data->[0];
			
			
			
			
			
			
			# $buffer1_ms2scan_low_signal_data->{$premz} = $buffer2_ms2scan_low_signal_data->{$premz};
			# $buffer2_ms2scan_low_signal_data->{$premz} = $buffer3_ms2scan_low_signal_data->{$premz};
			# $buffer3_ms2scan_low_signal_data->{$premz} = [@low_signals];
			
			$last2_ms2scan->{$premz} = $last1_ms2scan->{$premz};
			$last1_ms2scan->{$premz} = $ms2scan;
			$last2_title_line->{$premz} = $last1_title_line->{$premz};
			$last1_title_line->{$premz} = $title_line;
			$last2_rt_line->{$premz} = $last1_rt_line->{$premz};
			$last1_rt_line->{$premz} = $rt_line;
			$last2_pepmass_line->{$premz} = $last1_pepmass_line->{$premz};
			$last1_pepmass_line->{$premz} = $pepmass_line;
			#$last_premz->{$premz} = $premz;
			#$last_filtered_ion_arr->{$premz} = $filtered_ion_arr->{$premz};	
		}
	}
	close(FH);

	foreach my $key (keys %isolationlist){
			close($filehandle{$key});
	}

}
#close(MS1);
#close($filehandle{'0'});

sub pesudoms2{
	my ($center_ms2scan,$psms2,$seed,$premz,$title,$rest) = (@_);
	#print Dumper \@_;
	my $header = $ms2scan_header{$center_ms2scan};
	$header =~ s/TITLE=/$&${title}_/;
	#printf "\r\t",$header."\n";
	#if($#$psms2 + 1 >= $min_peaks_splited && $center_ms2scan == 31327){
	#print  "$#$psms2 + 1 >= $min_peaks_splited","\n";
	if($#$psms2 + 1 >= $min_peaks_splited){
		#my ($ions,$ratio_of_intensity_for_charge) = sort_ion_mz($psms2,$rest,$DIA_window{$premz}{lower},$DIA_window{$premz}{higher});
		my ($ions) = sort_ion_mz($psms2,$rest,$DIA_window{$premz}{lower},$DIA_window{$premz}{higher});
		#length ($premz_i) ?  do {$header =~ s/PEPMASS=[\d.]+/PEPMASS=${premz_i}/ } : do {$header =~ s/PEPMASS=[\d.]+/PEPMASS=${premz}/ };
		#print join "\t", $ms2scan, $premz_i, "#\n";
		#print $center_ms2scan, " ", $ratio_of_intensity_for_charge,"\n" if $center_ms2scan == 23083;
		
		# count the number of splited mgf
		my $key = $DIA_window{$premz}{size};
		
		# my @filesize = stat($filehandle{$key});
		
		print_mgf(
			$ms2scan_of_index{$premz}{$scan_premz_index->{$premz}-1},
			$ions,
			$header,
			$seed,
			#$ratio_of_intensity_for_charge > 0.95 ? [1] : [2..$max_charge],
			#length($premz_i) ? $filehandle{'0'} : $filehandle{$DIA_window{$premz}{size}}
			$filehandle{$key}
		);
		
	}
	#exit if $center_ms2scan > 31327;
}


sub print_mgf{
	my $scan = shift;
	my $ion_mz_ = shift;
	my $header = shift;
	my $seed = shift;
	my $fh = shift;
	
	# if(! defined $pre_charges){
		# $pre_charges = [$min_charge..$max_charge];
	# }
	#foreach my $z (@{$pre_charges}){
		#my $z = 2;
		
		$header =~ s/(TITLE=[^.]+\.\d+\.\d+)\.(\d+)?/${1}.0/;
		#print Dumper $h,$header;
		#my $fh = $h->{$z};
		print $fh $header;
		print $fh "SEED=$seed\n";
		map{
			print $fh  join " ",@{$_};
			print $fh  "\n";
		}@{$ion_mz_};
		print $fh "END IONS\n";
	#}
}



sub Spectrum_ms2_seed_finder{
	#my ($spec_ref,$sc_ind_curr,$ms1spec_ref,$premz) = @_;
	my ($sc_ind_curr,$premz) = @_;
	#print Dumper $spec_ref,$sc_ind_curr;
	my $center_sc = $num_of_data_points % 2 ? $ms2scan_of_index{$premz}{$sc_ind_curr - ($num_of_data_points - 1)/2} : $ms2scan_of_index{$premz}{$sc_ind_curr -$num_of_data_points/2};
	#my $ms1_center_sc = $ms1scan_of_index{$premz}{$sc_ind_curr - ($num_of_data_points - 1)/2};
	#my $center_ms2_noise = $spectra_noise{$center_sc};
	#next if $center_sc < 32096;
	my $center_ms2_tic = $ms2spectra_TIC{$center_sc};
	
	#if ($center_sc == 54279){
	
	#print $center_sc,"##";
	#if ($center_sc == 21091){
	#next unless $center_sc == 37002;
	my $center_spec = $num_of_data_points % 2 ? $ms2spectra_data{$premz}{($num_of_data_points + 1)/2} : $ms2spectra_data{$premz}{$num_of_data_points /2};
	
	#my $ms1_center_spec = $ms1spectra_data{$premz}{($num_of_data_points + 1)/2};
	#if($center_sc == 35062){print Dumper $center_spec;exit};
	#next;
	#print STDERR Dumper $center_spec if($center_sc == 25879);
	my $center_bp_mz = $basepeak_mz{$center_sc};
	my $center_bp_int = $basepeak_intensity{$center_sc}->[1];
	my $com_mz;
	
	foreach my $ind (1..$num_of_data_points){
		my $cmp_spec = $ms2spectra_data{$premz}{$ind};
		$com_mz = Spectrum_cmp($center_spec,$cmp_spec,$com_mz);
		#print STDERR Dumper '##ind\t'.$ind,$cmp_spec->[1];
	}
	#print "---$inx\n";
	#print Dumper $com_mz;print "\n##\n";
	#if($com_mz){
	#}
	#print "\t length of arr: $#$center_spec\t";
	my (%seed_list,%pearson_data,%candidates,%thres_p_data);	
	my $inx;
	my %low_sig;
	my %low_sig_for_each;
	my $f;
	foreach my $ref (@{$center_spec}){
		my ($mz,$int) = @{$ref};
		# if( $int < $center_bp_int * 0.01){
			# $f++;
			# #$low_sig_for_each{$mz} = $int if $int > $center_bp_int * 0.01;
			# next;
		# }
		
		$low_sig{$mz} = $int ;
		if($mz > 300 && defined $com_mz->{$mz}){
			#my $norm_ = $norm_data->{$mz};
			my ($boundary_start, $boundary_end,$n_nonzero) = boundary_index($com_mz->{$mz},$center_bp_int);
			#my $n_nonzero = abs($boundary_start-$boundary_end) + 1;
			#print STDERR join "\t",$mz,$n_nonzero,@{normalization_int($com_mz->{$mz})},"\n" if($center_sc == 25897);
			#if($n_nonzero >= 3 && $int > $center_bp_int * 0.01){
			# if($candadate_filter){#print "XX ";
				# if($n_nonzero >= 3 && $int > $basepeak_intensity_whole{$center_sc} * 0.01 ){ #print "$int > $center_bp_int * 0.01 ";
					# #my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$mz});	
					# $seed_list{$center_sc}{$mz} = [$int,$n_nonzero];
					# $inx++;
					
				# }
			# }else{
				if($n_nonzero >=  ($num_of_data_points - 2)){
					#my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$mz});	
					$seed_list{$center_sc}{$mz} = [$int,$n_nonzero];
					$inx++;
					#print STDERR join "\t","keep:",$mz,$int,@{$com_mz->{$mz}},"\n";
				}
				# else{
					# #print STDERR join "\t","rm:", $mz,$int,@{$com_mz->{$mz}},"\n";
				# }			
			# }
		}
	}
	#print "rm: $f peaks";
	# print Dumper $center_bp_int, \%seed_list if $center_sc == 19530;
	# exit if $center_sc == 19530;
	
	# foreach my $mz(sort{$a<=>$b}keys %{$seed_list{$center_sc}}){
		# print STDERR join "\t",$mz,@{$com_mz->{$mz}},"\n";
	# }
	#exit;
	#if($center_sc == 32096){
		#print STDERR Dumper \%seed_list,$com_mz,$center_spec;exit;
		#}
	
	my %finished_list;
	#my %mz_list_high_score;
	my %com_sig;
	my @seedlist;
	my %next_candidate_seeds;
	my $sec;
	my $dss;
	#print "$inx\n";
	#print join ":",(localtime())[1,0],"\n";
	if($inx > 20){
	my $cycle = 0;
	#my %thr = ( 1 => 20, 2=> 15, 3 => 5, 4 => 5);
		#print $dss++,"\n";
		while(1){
			#print $cycle,"\n";
			$cycle++;
			my @sorted_seed_list_by_intensity;
			if($sec){
				@sorted_seed_list_by_intensity = sort{
				#$seed_list{$center_sc}{$b}->[1] <=> $seed_list{$center_sc}{$a}->[1] ||
				$seed_list{$center_sc}{$b}->[0] <=> $seed_list{$center_sc}{$a}->[0]
				}grep{$low_sig{$_} > 0 and exists $next_candidate_seeds{$_}}(keys %{$seed_list{$center_sc}});
			}else{
				@sorted_seed_list_by_intensity = sort{
				#$seed_list{$center_sc}{$b}->[1] <=> $seed_list{$center_sc}{$a}->[1] ||
				$seed_list{$center_sc}{$b}->[0] <=> $seed_list{$center_sc}{$a}->[0]
				}grep{$low_sig{$_} > 0}(keys %{$seed_list{$center_sc}});
				my $seed = $sorted_seed_list_by_intensity[0];
				my $seed_int  = $seed_list{$center_sc}{$seed}->[0];
				#last if ($seed_int < $center_bp_int * 0.05 && $candadate_filter);
			}
			
			last if (scalar @sorted_seed_list_by_intensity) < 20 ;
			my $seed = $sorted_seed_list_by_intensity[0];
			my $seed_int  = $seed_list{$center_sc}{$seed}->[0];
			my $seed_n_nonzero = $seed_list{$center_sc}{$seed}->[1];
			#print STDERR Dumper $seed if($center_sc == 30415);
			# if($seed_int < $center_bp_int * 0.05 ){
				# #last;
				# undef %next_candidate_seeds;
				# delete $seed_list{$center_sc}{$seed};
				# next;
			# }
			
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$seed},$center_bp_int);		
			my $seed_int_list = [@{$com_mz->{$seed}}[$boundary_start..$boundary_end]];

			
			
			# if (! ifdecreaingOrincreasing($seed_int_list)){
				# #$fragment_ion_bp_not_matched_list{$rest_bp_mz} = $rest_bp_int_;
				# #print $cycle,"\t","#3\t",join "\t",@{normalization_int($seed_int_list)},"\n";
				# #print STDERR $seed."#4\n" if($center_sc == 30415);
				# delete $seed_list{$center_sc}{$seed};
				# undef %next_candidate_seeds;
				# next;
			# };			
			
			my @candidate_list;
			my @high_corr;
			my %p;
			my $p_cnt;
			undef %next_candidate_seeds;
			my $s = Statistics::KernelEstimation->new_gauss();
			
			#print scalar @sorted_seed_list_by_intensity,"\n";
			foreach my $mz (@sorted_seed_list_by_intensity){
				
				my $item_int_list = [@{$com_mz->{$mz}}[$boundary_start..$boundary_end]];
				my $mz_int = $seed_list{$center_sc}{$mz}->[0];
				my $mz_n_nonzero = $seed_list{$center_sc}{$mz}->[1];
				if($mz_int > $seed_int){
					next;
				}
				# if($seed_n_nonzero < $mz_n_nonzero){
					# next;
				# }
				#my $n_nonzero = scalar (grep{$_}@{$item_int_list});
				# if($n_nonzero < 3){
					# next;
				# }
				#if(! defined $pearson_data{$center_sc}{$seed}{$mz}){
					#print Dumper "\n\n",$com_mz if (! defined $com_mz->{$seed}) || (! defined $com_mz->{$mz}) ;
					#pearson_data{$center_sc}{$seed}{$mz} = MSE($seed_int_list,$item_int_list);
					$pearson_data{$center_sc}{$seed}{$mz} = pearson_corr($seed_int_list,$item_int_list);
					#$pearson_data{$center_sc}{$mz}{$seed} = $pearson_data{$center_sc}{$seed}{$mz};
				#}
				#if($cycle == 1){print STDERR $center_sc,"\t",$pearson_data{$center_sc}{$seed}{$mz},"\n";}
				
				
				# if($pearson_data{$center_sc}{$seed}{$mz} > 0.8 #0.9 25Da 0.85 10Da		
				# ){
					# push @candidate_list, $mz;
					#if($pearson_data{$center_sc}{$seed}{$mz} > 0){
						$p{$mz} = $pearson_data{$center_sc}{$seed}{$mz}; 
						$s->add_data($pearson_data{$center_sc}{$seed}{$mz}) if $pearson_data{$center_sc}{$seed}{$mz} >= 0.5;
						$p_cnt ++ if $pearson_data{$center_sc}{$seed}{$mz} > 0.9;
					#}
				# }elsif($pearson_data{$center_sc}{$seed}{$mz} < 0.5){
					# $next_candidate_seeds{$mz} = 1;
				# }
			}
			#print $seed,"\t",$p_cnt,"\n";
			
			if ($p_cnt < 10 ){
				#$fragment_ion_bp_not_matched_list{$rest_bp_mz} = $rest_bp_int_;
				#print $cycle,"\t","#3\t",join "\t",@{normalization_int($seed_int_list)},"\n";
				#print STDERR $seed."#4\n" if($center_sc == 30415);
				delete $seed_list{$center_sc}{$seed};
				undef %next_candidate_seeds;
				next;
			};	
			
			#print STDERR Dumper $s;
			my $w = 0.05;
			
			
			my $n;
			
			#my %pdf_data;
			my $min_pdf = 10000;
			my $min_pdf_p;
			
			my @sorted_mz_by_p = sort{$p{$a} <=> $p{$b}}grep{$p{$_}>=0.8}keys %p;
			
			my $max_p;
			my $max_pdf;
			my %pdf_data;
			# foreach my $mz(@sorted_mz_by_p){
				# $n = $p{$mz};
				# #my $pdf = $s->pdf($n,$w);
				# #print STDERR join "\t",$mz,$n,$pdf,"\n";
				# my $pdf = _pdf($s,$n,$w);
				# print STDERR join "\t","-",$mz,$n,$pdf,"\n";
				# $pdf_data{$n} = $pdf;
				# last if $n < 0;
				# #print STDERR join "\t",$center_sc, $seed, $mz, $n, $pdf,"\n";
				# if($max_pdf < $pdf && $n > 0.8){
					# $max_pdf = $pdf;
					# $max_p = $n;
				# }
			# }
			
		 
		#p 1 0.05*0.5 
		# ... now use bandwidth from above to find result
		 
			
			
			#print "max_p, $max_p","\n";
			my $_0point8 = $p{$sorted_mz_by_p[0]};
			my $pdf_0point8 = $s->pdf($_0point8,$w);
			#print STDERR "pdf_0point8: $_0point8   $pdf_0point8\n";
			foreach my $mz(@sorted_mz_by_p){
				$n = $p{$mz};
				
				my $pdf = $s->pdf($n,$w);
				#print  STDERR join "\t",$n,$pdf,$pdf_0point8,"\n";
				#print  STDERR join "\t",$center_sc, $seed, $mz, $n, $s->pdf($n,$w),"\n";
				if ($pdf_0point8 >= $pdf){
					$min_pdf = $pdf;
					$min_pdf_p = $n;
				}
				if($n > 0.9){
					last;
				}
					
			}
			#print join ":",(localtime())[0,1],"\n";exit;
			
			my $thres_p = $min_pdf_p < 0.8 ? 0.8 : $min_pdf_p;
			$thres_p_data{$seed} = $thres_p;
			 
			
			foreach my $mz (@sorted_seed_list_by_intensity){
				if($p{$mz} >= $thres_p){
					push @candidate_list,$mz;
				}elsif($p{$mz} < 0.5){
					$next_candidate_seeds{$mz} = 1;
				}
			}
			
			#print $thres_p,"\t",$#sorted_mz_by_p,"\t",$#candidate_list,"\n";
			if (scalar @candidate_list >= 10 ){
				# foreach my $m(keys %p){
					# print STDERR join "\t", $center_sc,$cycle,$seed,$m,$p{$m},"\n";
				# }
				#print STDERR  join "\t",$center_sc,$w,$max_p,$thres_p,"\n";
				$candidates{$center_sc}{$seed} = [@candidate_list] ;
				push @seedlist,$seed;
				$sec = 1;
				$cycle++;
				#print $cycle,"\n";
				for (@candidate_list){
					$low_sig{$_} = 0 ;
				}
				#last;
			}else{
				#undef %next_candidate_seeds;
				$sec = 0;
				delete $seed_list{$center_sc}{$seed};
			}
			if( $seed_int < $center_bp_int * 0.001  && $cycle >= 5){
				#$f++;
				#$low_sig_for_each{$mz} = $int if $int > $center_bp_int * 0.01;
				last ;
			}
		}
		#print "cycle: $cycle";
		
		foreach my $i (keys %low_sig){
			unless ($low_sig{$i} == 0){
				push @{$com_sig{$center_sc}},[$i,$low_sig{$i}];
			}
		}
	}else{
		$seed_list{$center_sc}{$center_bp_mz} = $center_bp_int;
		$candidates{$center_sc}{$center_bp_mz} = [];
		push @seedlist,$center_bp_mz;
		$com_sig{$center_sc} = [[0,0]];
	}
	
	#print "\tcom_sig_length : ",scalar @{$com_sig{$center_sc}};
	#print "\n";
	#print join ":",(localtime())[1,0],"\n";
	#print Dumper "seed",keys (%{$candidates{$center_sc}}),"\n";
	#print Dumper \@seedlist;
	#exit;
	#@@@@@@@@@
	#
	# mgf generation
	#
	#@@@@@@@@@
	
	
	if(scalar @seedlist == 0){
		$seed_list{$center_sc}{$center_bp_mz} = $center_bp_int;
		$candidates{$center_sc}{$center_bp_mz} = [];
		$com_sig{$center_sc} = [[0,0]];
		push @seedlist,$center_bp_mz;
	}
	#print Dumper \@seedlist;
	my $label;
	foreach my $seed (@seedlist){
	#foreach my $seed (keys %{$candidates{$center_sc}}){
		#print join "\t",$seed,scalar @{$candidates{$center_sc}{$seed}},defined $com_mz->{$seed},"\n";
		if(scalar @{$candidates{$center_sc}{$seed}} and defined $com_mz->{$seed}){
			my @peaklists_candidate = @{$candidates{$center_sc}{$seed}};
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$seed},$center_bp_int);	
			my $seed_mz_ints = $com_mz->{$seed};
			$seed_mz_ints = [@{$seed_mz_ints}[$boundary_start..$boundary_end]];			
			$label++;
			my %list;
			my $data_;
			my $seed_int;
			map{$list{$_} = 1}(@peaklists_candidate);
			foreach my $ref (@{$center_spec}){
				my ($mz,$int) = @{$ref};
				if($seed == $mz){
					$seed_int = $int;
				}
				if($list{$mz}){
					push @{$data_},$ref;
				}
			}
			

			my %u;
			foreach my $candidate (@peaklists_candidate){
				# my $n;
				my $candidate_ints = $com_mz->{$candidate};
				$candidate_ints = [@{$candidate_ints}[$boundary_start..$boundary_end]];
				foreach my $ref (@{$com_sig{$center_sc}}){
					my ($mz,$int) = @{$ref};
					next if $int > $seed_int && $mz > 300;
					next unless (defined $com_mz->{$mz} );
					next if defined $u{$mz};
					#if(! defined $pearson_data{$center_sc}{$mz}{$candidate}){
						my $item_int_list = [@{$com_mz->{$mz}}[$boundary_start..$boundary_end]];
						$pearson_data{$center_sc}{$candidate}{$mz} = pearson_corr($candidate_ints,$item_int_list);
						#$pearson_data{$center_sc}{$mz}{$candidate} = $pearson_data{$center_sc}{$candidate}{$mz};
						my $ppm_ = ppm($candidate + $massdiff_C12_C13,$mz);
					#}
					#if($pearson_data{$center_sc}{$candidate}{$mz} > 0.8){
					if($pearson_data{$center_sc}{$candidate}{$mz} > $thres_p_data{$seed} || ( abs($ppm_) < $ms2_precison && $pearson_data{$center_sc}{$candidate}{$mz} > 0.3 )){
						#print STDERR join "\t",$center_sc,":",$label,$candidate,$mz,"\n";
						push @{$data_}, [$mz,$int]; 
						$low_sig{$mz} = 0;
						$u{$mz} = 1;
					}elsif($pearson_data{$center_sc}{$candidate}{$mz} < 0.3){
						$u{$mz} = 1;
					}					
				}
			}
			#pesudoms2($center_sc,$data_,1+abs($boundary_start-$boundary_end),$premz,"res" . $label,[map{[50+$_*10,$seed_mz_ints->[$_]]}0..$#$seed_mz_ints]);
			# my $xx;$low_sig_for_each{$mz} = $int 
			foreach my $i (keys %low_sig_for_each){
				#unless ($low_sig{$i} == 0){
					push @{$data_},[$i,$low_sig_for_each{$i}];
				#}
			}
			#pesudoms2($center_sc,$xx,1,$premz,"xx" . $label,[[0,0]]);		
			#print join "\t",$center_sc,$seed,$data_,1,$premz,"res" . $label,[[0,0]],"\n";
			#print Dumper $center_sc,$data_,1,$premz,"res" . $label,[[0,0]];
			pesudoms2($center_sc,$data_,$seed,$premz,"res" . $label,[[0,0]]);				
		}else{
			pesudoms2($center_sc,$center_spec,1,$premz,"raw",[[0,0]]) ;
		}
	}
		#print join ":",(localtime())[1,0],"\n";
	#exit;
	#}
	
	#exit if($center_sc == 25879);
}		

sub _pdf(){
	my $d = shift;
	my $x = shift;
	my $w = shift;
	my $count  = $d->{sum_cnt};
	my $y = 0;
	for my $q ( @{ $d->{data} } ) {
	$y += $q->{cnt} * _gauss_pdf( $x, $q->{pos}, $w );
	}

	return $y/$count;
}
sub _gauss_pdf {
  my ( $x, $m, $s ) = @_;
  my $z = ($x - $m)/$s;
  return exp(-0.5*$z*$z)/( $s*sqrt( 2.0*$PI ) );
}

sub _box_pdf() {
	my ( $x, $m, $s ) = @_;
	if( $x < $m-0.5*$s || $x > $m+0.5*$s ) { return 0.0; }
	return 1.0/$s;
}

sub normalization_int{
	my $data = shift;
	my $max = max(@{$data});
	my @n = map{ $max ? $_ / $max : 0}@{$data};
	return ([@n]);
}

sub normal_function{ # p(x|w1) probability calculation
	my ($ave_,$std_,$x) = (@_);
	my $prob = exp(-($x-$ave_)**2/(2*($std_**2)))/(sqrt(2*$PI)*$std_);
	#print "exp(-($x-$ave_)/(2*($std_**2)))/(sqrt(2*$PI)*$std_)\n";
	#print $x,"\t",$prob,"\n";
	return $prob;
}

sub ppm{
	my ($observed_mz,$theoretical_mz) = (@_);
	return 1000000*($observed_mz-$theoretical_mz)/$theoretical_mz;
}


sub ifdecreaingOrincreasing{
	my $data = shift;
	my $normdata = normalization_int($data);
	my $min_ = min(@{$normdata});
	my $n;
	my $m;
	my $max_ind;
	# for(0..$#$normdata){
		# if($normdata->[$_] == 1){
			# $max_ind = $_;
			# last;
		# }
	# }
	# foreach my $i (1..$max_ind){
		# if($normdata->[$i] > $normdata->[$i-1]){
			# $n++;
		# }
	# }
	# foreach my $i ($max_ind+1..$#$normdata){
		# if($normdata->[$i] < $normdata->[$i-1]){
			# $m++;
		# }
	# }
	
	#my %h = ( 3 => 0.3, 4=> 0.4, 5 => 0.5, 6 => 0.6, 7=> 0.6); 
	#if($m+$n >= $#$normdata-1 && $min_ <  $h{$#$data + 1}){
	if($min_ < 0.8){
	#if($min_ < 0.5){
		# if($#$data % 2 && ($max_ind == ($#$data + 1)/2 || $max_ind == ($#$data - 1)/2 )){
			# return 1;
		# }elsif($#$data / 2 == $max_ind){
			return 1;
		# }
	}else{
		return 0;
	}
}


sub boundary_index{
    my $data = shift;
	my $bp = shift;
	my $m = $#$data/2;
	my $non_zero;
	#my $m = max_ind($data);
	
	my ($start,$end);
	
	for($m+1..$#$data){
		$non_zero++ if $data->[$_] > 0;
		if($data->[$_] == 0){
			$end = $_;
			last;
		}
	}
	
	for(0..$m-1){
		$non_zero++ if $data->[$_] > 0;
		if($data->[$_] == 0){
			$start = $_;
		}
	}
	$start = defined $start ? $start : 0;
    $end = defined $end ? $end : $#$data;
	#print join "\t",$start,$end,"\n";
    return ($start,$end,$non_zero+1);
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
			
			if(abs(ppm($mz1,$mz2)) <= $ms2_precison){
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

sub sort_ion_mz{
	my ($d1,$d2,$lower,$higher) = (@_);
	my %hash;
	my @ref_data = sort{$a->[0] <=> $b->[0]}grep{!$hash{$_->[0]}++ && length($_->[0]) }(@{$d1},@{$d2});
	return ([@ref_data]);
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

1;
