#!/usr/bin/env perl

use strict;
use Data::Dumper;
use File::Basename;
use POSIX qw (floor);
use lib dirname(__FILE__);
use peptide;

$|=1;

my $filename = shift;
my $ms2_ppm_tolerence = shift;

my @mgfs = glob "*${filename}_pseudo.mgf";#"Hela1ug_QC_Middle_DDA_150226_01.mgf";
if(scalar @mgfs == 0){
@mgfs = glob "*Calmz_${filename}.mgf";#"Hela1ug_QC_Middle_DDA_150226_01.mgf";
}
#@mgfs = glob "Calmz_${filename}*mgf" if ! (scalar @mgfs);#"Hela1ug_QC_Middle_DDA_150226_01.mgf";
#my @pinfiles = glob "${filename}*pin";# *pin
my @percolatorfiles = glob "${filename}.pin.target.pep.tsv";# *target.pep.tsv
my $qvalue = 0.01;

print Dumper \@mgfs,\@percolatorfiles,$ms2_ppm_tolerence; 
#exit;


my %res;
foreach my $percolatorfile (@percolatorfiles){
	open(RES,$percolatorfile)or die "$!\n $percolatorfile";
	print "Read $percolatorfile\n";
	while(<RES>){
		chomp;
		s/\r//;
		my ($t,$score,$q,$p,$pep,$pro) = split /\t/,$_;
		print "\r$t";
		if($q < $qvalue){
			my ($title,$charge) = ($t =~ /(.+[._](\d+))[_.]\d+$/);
			$title =~ s/.$//g;
			#$title =~ s/res\d+_|raw_|.$//g;
			#$title =~ s/^(.+_)(\d+)(_\d+)$/${1}.(${2}+1).$3/e;
			$res{$pep}{$charge}{"q"} ||= 1000;
			if($res{$pep}{$charge}{"q"} > $q){
				$res{$pep}{$charge}{"q"} = $q;
				$res{$pep}{$charge}{"title"} = $title;
				$res{$pep}{$charge}{"protein"} = $pro;
			}
		}
	}
	close(RES);
}

print "\n";
my %list_title;
foreach my $p(keys %res){
	foreach my $z(keys %{$res{$p}}){
		print "\r$p\t$z";
		$list_title{$res{$p}{$z}{"title"}} = [$p,$z,$res{$p}{$z}{"q"},$res{$p}{$z}{"protein"}];
	}
}
#print Dumper \%list_title;exit;
print "\n";

my $title;
my @ion_arr;

open(OUT,">$filename.DIANN_Lib.tsv");
#print OUT join "\t",("Peptide","Charge","q.value","Title","n_Y","n_B","PepLen","#match","#not_match");
print OUT join "\t",("ModifiedPeptide","PrecursorCharge","PrecursorMz","Tr_recalibrated","ProductMz","LibraryIntensity","ProteinID","ProteinName","FragmentCharge","FragmentType","FragmentSeriesNumber","FragmentLossType");
print OUT "\n";
foreach my $mgf (@mgfs){
print "Read: $mgf\n";
open(DIAMGF,$mgf);
my $tag;
my $rt;
my $mass;
while(<DIAMGF>){
	chomp;
	if(/^TITLE=(.*)/){
		$title = $1;
		$title =~ s/0$//;
		$tag = 0;
		undef @ion_arr;
		if(exists $list_title{$title}){
			print STDERR "\r$title";
			$tag = 1;
		}
	}
	if(/RTINSECONDS=([\d.]+)/){
		$rt = $1;
	}
	if(/PEPMASS=([\d.]+)/){
		$mass = $1;
	}
	if(/^\d+/ && $tag){
		my @line = split /\s+/,$_;
		push @ion_arr,[@line];
	}
	if(/END/ && $tag){
		my ($pep,$charge,$q,$pro) = @{$list_title{$title}};
		my ($proteinid,$proteinname) = ($pro =~ /^[^|]+\|([^|]+)\|([^ ]+).+/);
		$pep =~ s/:\(.+$//;
		$pep =~ s/^..|..$//g;
		
		my $mz = peptide::calmz($pep,$charge);
		$pep =~ s/C\[[^][]+\]/C(UniMod:4)/g;
		$pep =~ s/M\[[^][]+\]/C(UniMod:35)/g;

		#$pep =~ s/\.[\D].+$//g;
		#print $pep," ",$charge,"\n";
		my $theoretical_mz_table = peptide::fragmentation($pep,$charge);
		#print Dumper $list_title{$title};
		#print Dumper $theoretical_mz_table;
		my $deisotope_ion_arr = deisotope([@ion_arr]);
		#print Dumper $deisotope_ion_arr;
		my ($matched_mzs_type,$matched_int,$matched_mz,$topIntensity) = ms2spectra_match($theoretical_mz_table,$deisotope_ion_arr);
		foreach my $i (0..$#$matched_mzs_type){
			my ($SeriesNumber,$FragmentType,$FragmentCharge) = ($matched_mzs_type->[$i] =~ /(\d+)\.([by][^+]*)([+]+)/);
			my $loss = $FragmentType =~ /H2O/ ? "H2O" : $FragmentType =~ /NH3/ ? "NH3" : "noloss"; 
			$FragmentType =~ s/_(H2O|NH3)//;
			#print $FragmentCharge,"\t";
			$FragmentCharge = ($FragmentCharge =~ s/[+]//g);
			#print $FragmentCharge,"\n";
			#"ModifiedPeptide","PrecursorCharge","PrecursorMz","Tr_recalibrated","ProductMz","LibraryIntensity","ProteinID","ProteinName","FragmentCharge","FragmentType","FragmentSeriesNumber"
			print OUT join "\t",$pep,$charge,$mz,$rt/60,$matched_mz->[$i],$matched_int->[$i]/$topIntensity,$proteinid,$proteinname,$FragmentCharge,$FragmentType,$SeriesNumber,$loss;
			#print  join "\t",$pep,$charge,$mz,$rt/60,$matched_mz->[$i],$matched_int->[$i]/$topIntensity,$proteinid,$proteinname,$FragmentCharge,$FragmentType,$SeriesNumber,$loss;
			print OUT "\n";
			#print   "\n";
		}
		#print Dumper $pep,$matched_mzs_type;exit;
=head tmpName
$VAR1 = [
          'K.HAVSEGTK.A',
          '2',
          '5.58878e-05',
          'sp|Q96A08|H2B1A_HUMAN Histone H2B type 1-A OS=Homo sapiens OX=9606 GN=H2BC1 PE=1 SV=3'
        ];

$VAR1 = 'HAVSEGTK';
$VAR2 = [
          '1.y+',
          '2.b+',
          '2.y+',
          '3.y+',
          '4.b+',
          '4.y+',
          '5.y+',
          '5.b+',
          '6.b+',
          '6.y+',
          '7.y+'
        ];
$VAR3 = [
          '39212.80078125',
          '99045.84375',
          '17370.509765625',
          '33427.34765625',
          '6521.1240234375',
          '13455.931640625',
          '33404.61328125',
          '5949.8662109375',
          '13313.935546875',
          '42731.29296875',
          '60640.984375'
        ];


=cut 

	}
}
print STDERR "\n";
close(DIAMGF);
}
close(OUT);
print "\nOUTfile: $filename.DIANN_Lib.tsv\n";

sub deisotope{
	my $peaks = shift;
	
	my $proton = 1.00727646677;
	my $massdiff_C12_C13 = 1.0033548378;
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


sub ms2spectra_match{
	my ($pred,$obs) = (@_); # predicted, observed mz
	my @pred_list = sort{$a <=> $b}(keys %{$pred});
	my %uniq;
	my ($pred_index1,$obs_index2) = (0,0);
	my ($pred_len1,$obs_len2) = ($#pred_list,$#$obs);
	my $matched_mzs_type;
	my @ppms_frag;
	my $matched_ints;
	my $matched_mzs;
	my $topInt;
	#my @matched_ints;
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
			#print join "\t", $obs_mz2,$pred_mz1,$measured_ppm,"\n";
			if(abs($measured_ppm) <= $ms2_ppm_tolerence){
				
				my $predicted_mz_ion_type = $pred->{$pred_mz1};
				#if($predicted_mz_ion_type !~ /H2O|NH3/){	
					push @{$matched_mzs_type},$predicted_mz_ion_type;
					push @{$matched_ints},$obs_int2;
					push @{$matched_mzs},$obs_mz2;
					if($topInt < $obs_int2){
						$topInt = $obs_int2;
					}
				#}
				# $checked_obs_mz{$obs_index2} = 1;
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
	return ($matched_mzs_type,$matched_ints,$matched_mzs,$topInt);
}


sub ppm{
	my ($observed_mz,$theoretical_mz) = (@_);
	return 1000000*($observed_mz-$theoretical_mz)/$theoretical_mz;
}

sub tag{
	my $ref = shift;
	#my $pepLen = shift;
	my $return_type = shift; # type or coverage
	#my $pro = shift;
	my %ion_list;
	my %ion_for_coverage;
	my %uniq;
	
	#if($return_type eq "type"){
		foreach my $ion (@{$ref}){    
			my ($index,$fragment_type) = split /\./,$ion;
			#print $ion;
			#$iontype{$pro=~/REV_/?"Decoy":"Target"}{$fragment_type}++;
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
					$max_len = $tag if $max_len < $tag;
					$tag=1;
				}
			}
			$max_len = $tag if $max_len < $tag;
			#print Dumper $type, \@tags;
			$tag_len{$type} = $max_len; 
		}
		
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
			$max_len = $tag if $max_len < $tag;
			$tag=1;
		}
	}
	$max_len = $tag if $max_len < $tag;
	#print Dumper $type, \@tags;
	my $coverage_n = $max_len;
		
	return (\%tag_len,$coverage_n);	
}

