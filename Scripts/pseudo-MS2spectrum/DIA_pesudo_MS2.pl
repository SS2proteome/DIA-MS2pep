#!/usr/bin/env perl

use strict;
use Data::Dumper;
use POSIX qw(floor);
use List::Util qw(max min sum);
use POSIX qw(floor ceil round);
use File::Basename;
use lib dirname(__FILE__);
use peptide;
use Statistics::KernelEstimation;

my $PI = 4 * atan2(1, 1);
my $massdiff_C12_C13 = 1.0033548378;

$| = 1;

# ARGV
# Usage:
# DIA_pesudo_MS2.pl [filename] [dir] [ms2ppm]

my $mgf_file = "$ARGV[0].mgf";
my $mzML_file = "$ARGV[0].mzML";
my $dir = $ARGV[1];
my $ms2_tolerance = $ARGV[2];# ppm
my $window_file = "$ARGV[0].mzML.DIA_acquisition_window.txt";
my $ms2_low_mz = 140;
my $num_of_data_points;

if ( ! -d $dir ){
	mkdir( $dir ) or die "Can't creat $dir  $!\n";
}
my %DIA_window;
my %isolationlist;
my $number_of_windows;
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
	$number_of_windows++;
}
close(FH);


my $scan;
my $ms1scan;
my $ms2scan;
my %basepeak_intensity;
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

# READ mzML file
print "*****Read mzML file: $mzML_file!\n";
open(FH,$mzML_file) or die "$! $mzML_file\n";
while(<FH>){
	if(/<spectrum index=".*scan=(\d+).*/){
		$scan = $1;
		$number_of_scan = $scan;
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="([^"]+)"}){
		$scantime = $1;
		
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="([^"]+)" }){
		if($mslevel == 2){
			$basepeak_ms1->{$premz_target} = int($1);
		}
	}elsif(/name="ms level" value="(\d+)"/){
		$mslevel = $1;
		if($mslevel == 2){
			$ms2scan = $scan; 
			if($lastbasepeak_ms1->{$premz_target} == $basepeak_ms1->{$premz_target}){
				$basepeak_ms1_cnt->{$premz_target}++;
			}else{
				push @basepeak_ms1_hash,$basepeak_ms1_cnt->{$premz_target};
				$basepeak_ms1_cnt->{$premz_target} = 1;
			}
			$lastbasepeak_ms1->{$premz_target} = $basepeak_ms1->{$premz_target};
		}else{
			$ms1scan = $scan;
			#$currbasepeak_ms1 = $basepeak_ms1;
			$cycletime += $scantime - $lastscantime;
			$lastscantime = $scantime;
			
			$cyclecount++;
			
		}
	}elsif(m{accession="MS:1000827" name="isolation window target m/z" value="([\d.]+)" }){
		$premz_target = $1;
		printf "\r$ms1scan";
	}
}
close(FH);


print "\n";
print join "\t",$ARGV[0],$number_of_scan,$number_of_windows,$cycletime/$cyclecount,$cyclecount, $scantime,ceil(sum(@basepeak_ms1_hash)/($#basepeak_ms1_hash + 1)),"\n";

$num_of_data_points = ceil(sum(@basepeak_ms1_hash)/($#basepeak_ms1_hash + 1)) + 1;
$num_of_data_points = $num_of_data_points < 3 ? 3 : $num_of_data_points;

#  READ MGF FILE
print "\n\n*****Read mgf_file file: $mgf_file!\n";
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

open(FH,$mgf_file) or die "$! $mgf_file \n";

my %filehandle;
foreach my $key (keys %isolationlist){
	open($filehandle{$key},">${dir}/w_${key}_$mgf_file") or die "Can't write mgf file! $!\n";
}

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
		
		$buffer1_ms2scan_data->{$premz} = $buffer2_ms2scan_data->{$premz};
		$buffer2_ms2scan_data->{$premz} = $buffer3_ms2scan_data->{$premz};
		$buffer3_ms2scan_data->{$premz} = [@ion_arr];
		
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

sub pesudoms2{
	my ($center_ms2scan,$psms2,$seed,$premz,$title,$rest) = (@_);

	my $header = $ms2scan_header{$center_ms2scan};
	$header =~ s/TITLE=/$&${title}_/;

		my ($ions) = sort_ion_mz($psms2,$rest,$DIA_window{$premz}{lower},$DIA_window{$premz}{higher});
		my $key = $DIA_window{$premz}{size};
		print_mgf(
			$ms2scan_of_index{$premz}{$scan_premz_index->{$premz}-1},
			$ions,
			$header,
			$seed,
			$filehandle{$key}
		);
}


sub print_mgf{
	my $scan = shift;
	my $ion_mz_ = shift;
	my $header = shift;
	my $seed = shift;
	my $fh = shift;
	
	$header =~ s/(TITLE=[^.]+\.\d+\.\d+)\.(\d+)?/${1}.0/;
	print $fh $header;
	print $fh "SEED=$seed\n";
	map{
		print $fh  join " ",@{$_};
		print $fh  "\n";
	}@{$ion_mz_};
	print $fh "END IONS\n";

}

print "\n$mgf_file Finished!\n";


sub Spectrum_ms2_seed_finder{
	my ($sc_ind_curr,$premz) = @_;
	my $center_sc = $num_of_data_points % 2 ? $ms2scan_of_index{$premz}{$sc_ind_curr - ($num_of_data_points - 1)/2} : $ms2scan_of_index{$premz}{$sc_ind_curr -$num_of_data_points/2};
	my $center_ms2_tic = $ms2spectra_TIC{$center_sc};
	my $center_spec = $num_of_data_points % 2 ? $ms2spectra_data{$premz}{($num_of_data_points + 1)/2} : $ms2spectra_data{$premz}{$num_of_data_points /2};
	
	my $center_bp_mz = $basepeak_mz{$center_sc};
	my $center_bp_int = $basepeak_intensity{$center_sc}->[1];
	my $com_mz;
	
	foreach my $ind (1..$num_of_data_points){
		my $cmp_spec = $ms2spectra_data{$premz}{$ind};
		$com_mz = Spectrum_cmp($center_spec,$cmp_spec,$com_mz);
	}
	my (%seed_list,%pearson_data,%candidates,%thres_p_data);	
	my $inx;
	my %low_sig;
	my %low_sig_for_each;
	my $f;
	foreach my $ref (@{$center_spec}){
		my ($mz,$int) = @{$ref};
		if($int < $center_bp_int * 0.01){
			$f++;
			next;
		}
		
		$low_sig{$mz} = $int ;
		if($mz > 300 && defined $com_mz->{$mz}){
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$mz},$center_bp_int);
			my $n_nonzero = abs($boundary_start-$boundary_end) + 1;
				if($n_nonzero >= 3){
					$seed_list{$center_sc}{$mz} = [$int,$n_nonzero];
					$inx++;
				}

		}
	}

	my %finished_list;
	my %com_sig;
	my @seedlist;
	my %next_candidate_seeds;
	my $sec;
	my $dss;

	if($inx > 20){
	my $cycle = 0;
		while(1){
			#print $cycle,"\n";
			$cycle++;
			my @sorted_seed_list_by_intensity;
			if($sec){
				@sorted_seed_list_by_intensity = sort{
				$seed_list{$center_sc}{$b}->[0] <=> $seed_list{$center_sc}{$a}->[0]
				}grep{$low_sig{$_} > 0 and exists $next_candidate_seeds{$_}}(keys %{$seed_list{$center_sc}});
			}else{
				@sorted_seed_list_by_intensity = sort{
				$seed_list{$center_sc}{$b}->[0] <=> $seed_list{$center_sc}{$a}->[0]
				}grep{$low_sig{$_} > 0}(keys %{$seed_list{$center_sc}});
				my $seed = $sorted_seed_list_by_intensity[0];
				my $seed_int  = $seed_list{$center_sc}{$seed}->[0];
			}
			
			last if (scalar @sorted_seed_list_by_intensity) < 20 ;
			my $seed = $sorted_seed_list_by_intensity[0];
			my $seed_int  = $seed_list{$center_sc}{$seed}->[0];
			my $seed_n_nonzero = $seed_list{$center_sc}{$seed}->[1];
			
			my ($boundary_start, $boundary_end) = boundary_index($com_mz->{$seed},$center_bp_int);		
			my $seed_int_list = [@{$com_mz->{$seed}}[$boundary_start..$boundary_end]];

			my @candidate_list;
			my @high_corr;
			my %p;
			my $p_cnt;
			undef %next_candidate_seeds;
			my $s = Statistics::KernelEstimation->new_gauss();
			
			foreach my $mz (@sorted_seed_list_by_intensity){
				
				my $item_int_list = [@{$com_mz->{$mz}}[$boundary_start..$boundary_end]];
				my $mz_int = $seed_list{$center_sc}{$mz}->[0];
				my $mz_n_nonzero = $seed_list{$center_sc}{$mz}->[1];
				if($mz_int > $seed_int){
					next;
				}
				$pearson_data{$center_sc}{$seed}{$mz} = pearson_corr($seed_int_list,$item_int_list);
				$p{$mz} = $pearson_data{$center_sc}{$seed}{$mz}; 
				$s->add_data($pearson_data{$center_sc}{$seed}{$mz}) if $pearson_data{$center_sc}{$seed}{$mz} >= 0.5;
				$p_cnt ++ if $pearson_data{$center_sc}{$seed}{$mz} > 0.9;

			}

			if ($p_cnt < 10 ){

				delete $seed_list{$center_sc}{$seed};
				undef %next_candidate_seeds;
				next;
			};	

			my $w = 0.05;
			my $n;
			
			my $min_pdf = 10000;
			my $min_pdf_p;
			
			my @sorted_mz_by_p = sort{$p{$a} <=> $p{$b}}grep{$p{$_}>=0.8}keys %p;
			
			my $max_p;
			my $max_pdf;
			my %pdf_data;
			my $_0point8 = $p{$sorted_mz_by_p[0]};
			my $pdf_0point8 = $s->pdf($_0point8,$w);
			foreach my $mz(@sorted_mz_by_p){
				$n = $p{$mz};
				
				my $pdf = $s->pdf($n,$w);
				if ($pdf_0point8 >= $pdf){
					$min_pdf = $pdf;
					$min_pdf_p = $n;
				}
				if($n > 0.9){
					last;
				}
					
			}
			
			my $thres_p = $min_pdf_p < 0.8 ? 0.8 : $min_pdf_p;
			$thres_p_data{$seed} = $thres_p;
			 
			
			foreach my $mz (@sorted_seed_list_by_intensity){
				if($p{$mz} >= $thres_p){
					push @candidate_list,$mz;
				}elsif($p{$mz} < 0.5){
					$next_candidate_seeds{$mz} = 1;
				}
			}
			
			if (scalar @candidate_list >= 10 ){

				$candidates{$center_sc}{$seed} = [@candidate_list] ;
				push @seedlist,$seed;
				$sec = 1;

				for (@candidate_list){
					$low_sig{$_} = 0 ;
				}
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
		push @seedlist,$center_bp_mz;
		$com_sig{$center_sc} = [[0,0]];
	}

	if(scalar @seedlist == 0){
		$seed_list{$center_sc}{$center_bp_mz} = $center_bp_int;
		$candidates{$center_sc}{$center_bp_mz} = [];
		$com_sig{$center_sc} = [[0,0]];
		push @seedlist,$center_bp_mz;
	}
	my $label;
	foreach my $seed (@seedlist){
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
						my $item_int_list = [@{$com_mz->{$mz}}[$boundary_start..$boundary_end]];
						$pearson_data{$center_sc}{$candidate}{$mz} = pearson_corr($candidate_ints,$item_int_list);
						my $ppm_ = ppm($candidate + $massdiff_C12_C13,$mz);
						if($pearson_data{$center_sc}{$candidate}{$mz} > $thres_p_data{$seed} || ( abs($ppm_) < $ms2_tolerance && $pearson_data{$center_sc}{$candidate}{$mz} > 0.3 )){

						push @{$data_}, [$mz,$int]; 
						$low_sig{$mz} = 0;
						$u{$mz} = 1;
					}elsif($pearson_data{$center_sc}{$candidate}{$mz} < 0.3){
						$u{$mz} = 1;
					}					
				}
			}

			foreach my $i (keys %low_sig_for_each){
					push @{$data_},[$i,$low_sig_for_each{$i}];
			}
			pesudoms2($center_sc,$data_,1,$premz,"res" . $label,[[0,0]]);				
		}else{
			pesudoms2($center_sc,$center_spec,1,$premz,"raw",[[0,0]]) ;
		}
	}
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


sub ppm{
	my ($observed_mz,$theoretical_mz) = (@_);
	return 1000000*($observed_mz-$theoretical_mz)/$theoretical_mz;
}

sub boundary_index{
    my $data = shift;
	my $bp = shift;
	my $m = $#$data/2;
	
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
    return ($start,$end);
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

			if(abs(ppm($mz1,$mz2)) <= $ms2_tolerance){
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

sub sort_ion_mz{
	my ($d1,$d2,$lower,$higher) = (@_);
	my %hash;
	my @ref_data = sort{$a->[0] <=> $b->[0]}grep{!$hash{$_->[0]}++ && length($_->[0])}(@{$d1},@{$d2});
	return ([@ref_data]);
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

