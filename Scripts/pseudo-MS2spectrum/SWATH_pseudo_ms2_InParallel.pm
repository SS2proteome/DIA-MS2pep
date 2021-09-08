package SWATH_pseudo_ms2_InParallel;

# charge_script [$file] [$charge_dir] [ms2 ppm] [data points] [min_pearson] [min_res_pearson] [min_peaks_splited]

use strict;
use Data::Dumper;
use List::Util qw(max min sum);
use POSIX qw(floor ceil round);
use File::Basename;
use lib dirname(__FILE__);
use peptide;
use Statistics::KernelEstimation;
my $PI = 4 * atan2(1, 1);
my $massdiff_C12_C13 = 1.0033548378;


$| = 1;

my $min_peaks_splited = 20;
my $mgf_file;
my $ms1_file;
my $mzML_file;
my $charge_dir;
my $ms2_precison;	
my $num_of_data_points;
my $window_file;
my %ms2toms1;
my $start_premz ;
my $end_premzmz ;
my $windowsize ;
my $ms2_low_mz = 140;
my $ms2_min_intensity = 2;
my $ms2scan;
my %SWATH_window;

my %isolationlist;
my $number_of_windows;

my %ms2toms1;
my $start_premz ;
my $end_premzmz ;
my $windowsize ;

my $scan;
my $ms1scan;
my $ms2scan;
my %premz_list;
#my %ms1toms2;
my $cycle;
my $index;
my $experiment;
my $isolationwindow;
#my %index_ms2scan;
##my %index_ms1scan;
#my $index_ms1;
my %basepeak_intensity;
my $mslevel;
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


my $title_line;
my $rt_line;
my $pepmass_line;
my $ion_lines;
my $premz;
my $experiment;
my $cycle;
my $cnt;
#my %basepeak_intensity;
my %basepeak_mz;
my %filehandle;

my $basepeak;
my $median_peak;
my @peaks;
my @low_signals;
my $scan_premz_index;
my %ms2scan_of_index;
my %ms1scan_of_index;
my %ms2spectra_bin_data;
my %ms2spectra_data;
my %ms2spectra_TIC;
my %ms1spectra_data;
my %ms2scan_header;
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


sub demultiplexing {
	$mgf_file = "@_[0]";;
	$charge_dir = "@_[1]";
	$ms2_precison = "@_[2]";
	$num_of_data_points = "@_[3]";
	$window_file = "@_[4].mzML.SWATH_acquisition_window.txt";


	open(FH,$window_file) or die "$!\n $window_file";
	while(<FH>){
		chomp;
		next if $. == 1;
		my @line = split /\t/,$_;
		if($line[-1] >100){$line[-1] = 100;}
		$SWATH_window{$line[-2]}{size} = $line[-1];
		$SWATH_window{$line[-2]}{lower} = $line[1];
		$SWATH_window{$line[-2]}{higher} = $line[2];
		$isolationlist{$line[-1]} = 1;
		$number_of_windows = $line[0];
	}
	close(FH);
 
	#  READ MGF FILE


	open(FH,"${charge_dir}/$mgf_file") or die "$! ${charge_dir}/$mgf_file \n";


	#my %filehandle_num;
	foreach my $key (keys %isolationlist){
			open($filehandle{$key},">${charge_dir}/w_${key}_$mgf_file") or die "Can't write mgf file! ${charge_dir}/w_${key}_$mgf_file $!\n";
	}



	while(<FH>){
		chomp;

		if(/^TITLE=([^.]+).([^.]+).([^.]+).([^.]+)/){#TITLE=Site10_DynRge_S1_Day_1.1.1.2
			undef @peaks;
			undef @low_signals;
			$title_line = $_;
			$ms2scan = $2;#.".".$3.".".$4;
			($experiment,$cycle) = ($4,$3);
			$cnt = 0;
			printf "\r$ms2scan";
		}elsif(/^RTIN/){
			$rt_line = $_;
		}elsif(/^PEPMASS=([\d.]+)/){
			$pepmass_line = $_;
			$premz = $1;
		}elsif(/^\d+/){
			my @line = split /\s+/,$_;
			if($line[1] < 1){
				push @low_signals, [@line];
			}		
			# if($basepeak_intensity{$ms2scan} < $line[1]){
					# $basepeak_intensity{$ms2scan} = $line[1];
					# $basepeak_mz{$ms2scan} = $line[0];
			# }
			if ($line[0] > $ms2_low_mz && 
			($line[0] < $SWATH_window{$premz}{lower} -1  || $line[0] > $SWATH_window{$premz}{higher} + 1)# &&
			#$line[1] > $ms2_min_intensity
			){
				#$basepeak = $basepeak < $line[1] ? $line[1] : $basepeak;
				push @peaks,[@line];
				$cnt++;
				
			}
		}elsif(/END\sIONS/){
			#$median_peak = peak_feature(\@peaks);
			#if($cnt >= 50 #&& #($median_peak / $basepeak) < 1/3  &&
			#){
				#next if $cnt < 30 ;
			
				#my ($filtered_ion_arr,$noise_) = noise_filtering(\@peaks,$basepeak_intensity{$ms2scan});
				my $filtered_ion_arr;
				my $noise_ = 100000000;
				
				if($buffer1_ms2scan_data->{$premz}->[0]){
					foreach my $ref (@{$buffer2_ms2scan_data->{$premz}}){
						my ($mz,$i) = @{$ref};
							push @{$filtered_ion_arr->{$premz}},$ref;
					}
				}
				
				if (defined $filtered_ion_arr->{$premz}){
				
					#print Dumper $ms2toms1{$ms2scan};
					#my $ms1peaks = MS1file($ms2scan);
					#$spectra_noise{$last_ms2scan} = $last_noise_;
					$scan_premz_index->{$premz}++;
					$ms2scan_of_index{$premz}{$scan_premz_index->{$premz}} = $last2_ms2scan->{$premz};		
					$ms2scan_header{$last2_ms2scan->{$premz}} = "BEGIN IONS\n$last2_title_line->{$premz}\n$last2_rt_line->{$premz}\n$last2_pepmass_line->{$premz}\n";
					
					for(1..$num_of_data_points-1){
						$ms2spectra_data{$premz}{$_} = $ms2spectra_data{$premz}{$_+1};
						#$ms1spectra_data{$premz}{$_} = $ms1spectra_data{$premz}{$_+1};
					}
					$ms2spectra_data{$premz}{$num_of_data_points} = $filtered_ion_arr->{$premz};
					if($scan_premz_index->{$premz} >= $num_of_data_points){# && $ms2scan_of_index{$premz}{$scan_premz_index->{$premz} - ($num_of_data_points - 1)/2} == 35272){
						Spectrum_ms2_seed_finder(
						#print Dumper (
							#$ms2spectra_data{$last_premz},
							$scan_premz_index->{$premz},
							#$ms1spectra_data{$premz},
							$premz
						);
						#exit;
					}
				}
					$buffer1_ms2scan_data->{$premz}  = $buffer2_ms2scan_data->{$premz} ;
					$buffer2_ms2scan_data->{$premz}  = $buffer3_ms2scan_data->{$premz} ;
					$buffer3_ms2scan_data->{$premz}  = [@peaks];
					
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
					
		}
	}

	close(FH);

	foreach my $key (keys %isolationlist){
		close($filehandle{$key});
	}
}



sub pesudoms2{
	my ($center_ms2scan,$psms2,$seed,$premz,$title,$rest) = (@_);
	my $header = $ms2scan_header{$center_ms2scan};
	$header =~ s/TITLE=/$&${title}_/;
	#print $header."\n";
	
	if($#$psms2 + 1 >= $min_peaks_splited){
		my ($ions) = sort_ion_mz($psms2,$rest,$SWATH_window{$premz}{lower},$SWATH_window{$premz}{higher});
		#length ($premz_i) ?  do {$header =~ s/PEPMASS=[\d.]+/PEPMASS=${premz_i}/ } : do {$header =~ s/PEPMASS=[\d.]+/PEPMASS=${premz}/ };
		#print join "\t", $ms2scan, $premz_i, "#\n";
		
		my $key = $SWATH_window{$premz}{size};
		print_mgf(
			$ms2scan_of_index{$premz}{$scan_premz_index->{$premz}-1},
			$ions,
			$header,
			$seed,
			#$ratio_of_intensity_for_charge > 0.95 ? [1] : [2..$max_charge],
			$filehandle{$key}
		);
	}
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
		$header =~ s/(TITLE=[^.]+\.\d+\.\d+)\.(\d+)?/$1.0/;
		#print Dumper $h,$header;
		
		#my $fh = $h->{$z};
		print $fh $header;
		print $fh "SEED=$seed\n";
		#print $fh "CHARGE=${z}+\n";
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
#print Dumper ($sc_ind_curr,$premz);exit;
	my $center_sc = $num_of_data_points % 2 ? $ms2scan_of_index{$premz}{$sc_ind_curr - ($num_of_data_points - 1)/2} : $ms2scan_of_index{$premz}{$sc_ind_curr -$num_of_data_points/2};

	
	my $center_spec = $num_of_data_points % 2 ? $ms2spectra_data{$premz}{($num_of_data_points + 1)/2} : $ms2spectra_data{$premz}{$num_of_data_points /2};

	my $center_bp_mz = $basepeak_mz{$center_sc};
	my $center_bp_int = $basepeak_intensity{$center_sc};
	my $com_mz;
	foreach my $ind (1..$num_of_data_points){
		my $cmp_spec = $ms2spectra_data{$premz}{$ind};
		$com_mz = Spectrum_cmp($center_spec,$cmp_spec,$com_mz);
		#print STDERR Dumper '##ind\t'.$ind,$center_spec,$cmp_spec,$com_mz;
	}
	
	my (%seed_list,%pearson_data,%candidates,%thres_p_data);	
	my $inx;
	my %low_sig;
	#my $center_noise = $center_spec->[$#$center_spec]->[1];
	#print "center_noise:$center_noise\n";
	foreach my $ref (@{$center_spec}[0..$#$center_spec-1]){
		my ($mz,$int) = @{$ref};
		$low_sig{$mz} = $int;
		if($mz > 300 && defined $com_mz->{$mz}){
			#my $norm_ = $norm_data->{$mz};
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$mz});
			my $n_nonzero = abs($boundary_start-$boundary_end) + 1;
			#if($n_nonzero >= 3 && $int > $center_bp_int * 0.01 ){
			if($n_nonzero >= 3 && $int > 3){
				#my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$mz});	
				$seed_list{$center_sc}{$mz} = [$int,$n_nonzero];
				$inx++;
			}
			
		}
	}
	
	my %finished_list;
	#my %mz_list_high_score;
	my %com_sig;
	my %next_candidate_seeds;
	my $sec;
	my $cy;
	if($inx > 20){
	
		
		while(1){
			# my @sorted_seed_list_by_intensity = sort{
				# #$seed_list{$center_sc}{$b}->[1] <=> $seed_list{$center_sc}{$a}->[1] ||
				# $seed_list{$center_sc}{$b}->[0] <=> $seed_list{$center_sc}{$a}->[0]
			# }grep{$low_sig{$_} > 0}(keys %{$seed_list{$center_sc}});
			$cy++;
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
				#last if ($seed_int < 5 );
			}
			
			#print Dumper \@sorted_seed_list_by_intensity;
			last if (scalar @sorted_seed_list_by_intensity) < 10 ;
			
			my $seed = $sorted_seed_list_by_intensity[0];
			my $seed_int  = $seed_list{$center_sc}{$seed}->[0];
			
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$seed});		
			my $seed_int_list = [@{$com_mz->{$seed}}[$boundary_start..$boundary_end]];

			if ($boundary_end - $boundary_start + 1 < 3 ){
				#$fragment_ion_bp_not_matched_list{$rest_bp_mz} = $rest_bp_int_;
				#print $cycle,"\t","#3\n";
				undef %next_candidate_seeds;
				delete $seed_list{$center_sc}{$seed};
				next;
			};		
			
			if (! ifdecreaingOrincreasing($seed_int_list)){
				#$fragment_ion_bp_not_matched_list{$rest_bp_mz} = $rest_bp_int_;
				#print $cycle,"\t","#3\t",join "\t",@{normalization_int($seed_int_list)},"\n";
				undef %next_candidate_seeds;
				delete $seed_list{$center_sc}{$seed};
				next;
			};			
			
			my @candidate_list;
			my %p;
			undef %next_candidate_seeds;
			#++$cycle;
			my $p_cnt;
			my $s = Statistics::KernelEstimation->new_gauss();
			
			foreach my $mz (@sorted_seed_list_by_intensity){
				
				my $item_int_list = [@{$com_mz->{$mz}}[$boundary_start..$boundary_end]];
				my $mz_int = $seed_list{$center_sc}{$mz}->[0];
				if($mz_int > $seed_int){
					next;
				}
					$pearson_data{$center_sc}{$seed}{$mz} = pearson_corr($seed_int_list,$item_int_list);
					$p{$mz} = $pearson_data{$center_sc}{$seed}{$mz}; 
					$s->add_data($pearson_data{$center_sc}{$seed}{$mz}) if $pearson_data{$center_sc}{$seed}{$mz} > 0.5;;
					$p_cnt ++ if $pearson_data{$center_sc}{$seed}{$mz} > 0.8;
			}
			
			if ($p_cnt < 10 ){
				delete $seed_list{$center_sc}{$seed};
				undef %next_candidate_seeds;
				next;
			};	
			
			my $w = 0.05;
			my $n;
			
			my %pdf_data;
			my $min_pdf = 10000;
			my $min_pdf_p;
			
			my @sorted_mz_by_p = sort{$p{$b} <=> $p{$a}}keys %p;
			my $max_p;
			my $max_pdf_all_ions;
			my $max_pdf;
			foreach my $mz(@sorted_mz_by_p){
				$n = $p{$mz};
				my $pdf = $s->pdf($n,$w);
				$pdf_data{$n} = $pdf;
				
				if($max_pdf < $pdf && $n > 0.8){
					$max_pdf = $pdf;
					$max_p = $n;
				}
				
				if($max_pdf_all_ions < $pdf){
					$max_pdf_all_ions = $pdf;
				}
			}
			#print "max_p, $max_p","\n";
			
			unless ($max_pdf_all_ions == $max_pdf ){
				#$fragment_ion_bp_not_matched_list{$rest_bp_mz} = $rest_bp_int_;
				#print $cycle,"\t","#3\t",join "\t",@{normalization_int($seed_int_list)},"\n";
				#print STDERR $seed."#4\n" if($center_sc == 30415);
				delete $seed_list{$center_sc}{$seed};
				undef %next_candidate_seeds;
				next;
			};	
			
			foreach my $mz(@sorted_mz_by_p){
				#my $x = $n->{"pos"};
				$n = $p{$mz};
				
				my $pdf = $pdf_data{$n};
				#print  STDERR join "\t",$center_sc, $seed, $mz, $n, $s->pdf($n,$w),"\n";
				if($n >= $max_p){
					next;
				}
				
				#$pdf_data{$pdf} = $mz;
				if ($min_pdf > $pdf && $n > 0.8){
					$min_pdf = $pdf;
					$min_pdf_p = $n;
				}
				if($n < 0.8){
					last;
				}
			}
			
			#print STDERR $center_sc,"\t",$w, "\t",$thres_p,"\n";
			my $thres_p = $min_pdf_p < 0.8 ? 0.8 : $min_pdf_p;
			
			foreach my $mz (@sorted_mz_by_p){
				$n = $p{$mz};
				#print STDERR join "\t",$center_sc, $seed, $mz, $n, $pdf_data{$n},"\n";
				if($p{$mz} >= $thres_p){
					push @candidate_list,$mz;
				}elsif($p{$mz} < 0.5){
					$next_candidate_seeds{$mz} = 1;
				}
			}
			
			if (scalar @candidate_list >= 10 ){
				# foreach my $m(keys %p){
					# print STDERR join "\t", $center_sc,$seed,$m,$p{$m},"\n";
				# }
				#print STDERR $center_sc,"\t",$w, "\t",$thres_p,"\n";
				
				$candidates{$center_sc}{$seed} = [@candidate_list] ;
				$sec = 1;
				#$cycle++;
				for (@candidate_list){
					$low_sig{$_} = 0 ;
				}
				#last;
			}else{
				$sec = 0;
				delete $seed_list{$center_sc}{$seed};
			}
		}
		foreach my $i (keys %low_sig){
			unless ($low_sig{$i} == 0){
				push @{$com_sig{$center_sc}},[$i,$low_sig{$i}];
			}
		}
	}else{
		$seed_list{$center_sc}{$center_bp_mz} = $center_bp_int;
		$candidates{$center_sc}{$center_bp_mz} = [];
		$com_sig{$center_sc} = [[0,0]];
	}
	
	#print "\t\tC::$cy     ";
	#@@@@@@@@@
	#
	# mgf generation
	#
	#@@@@@@@@@
	my $label;
	foreach my $seed (keys %{$candidates{$center_sc}}){
		if(scalar @{$candidates{$center_sc}{$seed}} and defined $com_mz->{$seed}){
			my @peaklists_candidate = @{$candidates{$center_sc}{$seed}};
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$seed});	
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
						my $ppm_ = ppm($candidate + $massdiff_C12_C13,$mz);
						#$pearson_data{$center_sc}{$mz}{$candidate} = $pearson_data{$center_sc}{$candidate}{$mz};
					#}
					if($pearson_data{$center_sc}{$candidate}{$mz} > $thres_p_data{$seed} || ( abs($ppm_) < $ms2_precison && $pearson_data{$center_sc}{$candidate}{$mz} > 0.3 )){
						push @{$data_}, [$mz,$int]; 
						$u{$mz} = 1;
					}elsif($pearson_data{$center_sc}{$candidate}{$mz} < 0.3){
						$u{$mz} = 1;
					}					
				}
			}
			#pesudoms2($center_sc,$data_,1+abs($boundary_start-$boundary_end),$premz,"res" . $label,[map{[50+$_*10,$seed_mz_ints->[$_]]}0..$#$seed_mz_ints]);				
			pesudoms2($center_sc,$data_,$seed,$premz,"res" . $label,[[0,0]]);				
		}else{
			pesudoms2($center_sc,$center_spec,1,$premz,"raw",[[0,0]]) ;
		}
	}
}

sub ppm{
	my ($observed_mz,$theoretical_mz) = (@_);
	return 1000000*($observed_mz-$theoretical_mz)/$theoretical_mz;
}


sub normalization_int{
	my $data = shift;
	my $max = max(@{$data});
	my @n = map{ $max ? $_ / $max : 0}@{$data};
	return [@n];
}

sub ifdecreaingOrincreasing{
	my $data = shift;
	my $normdata = normalization_int($data);
	my $min_ = min(@{$normdata});
	# my $n;
	# my $m;
	# my $max_ind;
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
	#if($m+$n >= $#$normdata-2 && $min_ < 0.5){
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
			
			if(abs(2000000*($mz1-$mz2)/($mz2+$mz1)) <= $ms2_precison){
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
	#my $bp = shift;
	my $m = $#$data/2;
	#my $m = max_ind($data);
	
	my ($start,$end);
	
	for($m+1..$#$data){
		if($data->[$_] == 0){
			$end = $_;
			last;
		}
	}
	
	for(0..$m-1){
		if($data->[$_] == 0){
			$start = $_;
		}
	}
	$start = defined $start ? $start : 0;
    $end = defined $end ? $end : $#$data;
	#print join "\t",$start,$end,"\n";
    return ($start,$end);
}

sub sort_ion_mz{
	my ($d1,$d2,$lower,$higher) = (@_);
	my %hash;
	my @ref_data = sort{$a->[0] <=> $b->[0]}grep{!$hash{$_->[0]}++}(@{$d1},@{$d2});
	# my $int1;
	# my $cnt1;
	# my $int2;
	# my $cnt2;
	# map{
		# if($_->[0] < $lower){
			# $int1 += $_->[1];
			# $cnt1++;
		# }
		# if($_->[0] > $higher){
			# $int2 += $_->[1];
			# $cnt2++;
		# }
	# }@$d1;
	
	#return ([@ref_data], ($int1/($int1+$int2)));
	return ([@ref_data]);
	#return (\@ref_data, $cnt1." ".$cnt2." ".$lower." ".$higher);
}

sub pearson_corr{
    my ($ref_a, $ref_b) = @_;
    my @x = @{$ref_a};
    my @y = @{$ref_b};
	my $mean_x = sum(@x)/($#x+1);
	my $mean_y = sum(@y)/($#y+1);
	#my $len = $#x + 1;
	my $correlation;
	if($#x==$#y){
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
	}
	return $correlation;
}

1;
