package peptide;

# 20180621

use strict;
use List::Util qw(max min sum);
use Data::Dumper;

# mono_mass_list
my $O = 15.99491463;
my $H = 1.007825035;
my $N = 14.003074;
my $proton = 1.00727646677;
my $H2O = $H*2 + $O;
my $NH3 = $H*3 + $N;
my $Phos_loss = 79.9663 + $H2O;
my $massdiff_C12_C13 = 1.0033548378;


my %mono_aa;
$mono_aa{'A'} = 71.037114;
$mono_aa{'R'} = 156.101111;
$mono_aa{'N'} = 114.042927;
$mono_aa{'D'} = 115.026943;
$mono_aa{'C'} = 103.009185;
$mono_aa{'C_'} = 160.030649; # Carboxyamid
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
$mono_aa{'M_'} = 131.040485 + $O; # Oxidation
$mono_aa{'F'} = 147.068414;
$mono_aa{'P'} = 97.052764;
$mono_aa{'S'} = 87.032028;
$mono_aa{'S_'} = 166.9984; # Phos
$mono_aa{'T'} = 101.047679; 
$mono_aa{'T_'} = 181.0140; # Phos
$mono_aa{'U'} = 150.95363;
$mono_aa{'W'} = 186.079313;
$mono_aa{'Y'} = 163.06332;
$mono_aa{'Y_'} = 243.0297; # Phos
$mono_aa{'V'} = 99.068414;

#<aminoacid_modification aminoacid="S" massdiff="79.9663" mass="166.9984" variable="Y"/>
#<aminoacid_modification aminoacid="T" massdiff="79.9663" mass="181.0140" variable="Y"/>
#<aminoacid_modification aminoacid="Y" massdiff="79.9663" mass="243.0297" variable="Y"/>

sub amino{
	return \%mono_aa;
}


sub fragmentation {

	my $sequence = shift;
	my $precursor_charge = shift;
	my $pep = $sequence;
	
	#my $sequence = "GHRPLDK";
	
	#$sequence =~ s/M\[147.0354\]/M_/g;
	#$sequence =~ s/C\[160.0307\]/C_/g;
	$sequence =~ s/S\[166.9984\]/S_/g;
	$sequence =~ s/T\[181.0140\]/T_/g;
	$sequence =~ s/Y\[243.0297\]/Y_/g;
	
	while($sequence =~ /(.)\[([^][]+)\]/g){
		my $aa = $1;
		my $massdiff = $2;
		$mono_aa{${aa}."_"} = $massdiff;
		
		$sequence =~ s/$aa\[[^][]+\]/${aa}_/;
	}
	#print $sequence,"\n";
	my %mz;
	my $mz_b;
	my $mz_y;
	my @seq_arr = split /(?=[[:upper:]]_?)/,$sequence;
	#print Dumper \@seq_arr;
	my $n;
	my $y_pep;
	my $b_pep;
	foreach my $index (0..$#seq_arr){
		$b_pep .= $seq_arr[$index];
		
		$mz_b += $mono_aa{$seq_arr[$index]};
		$mz{"b+"}{++$n} = $mz_b;
		
		$y_pep .= $seq_arr[$#seq_arr - $index];
		$mz_y += $mono_aa{$seq_arr[$#seq_arr - $index]};
		$mz{"y+"}{$n} = $mz_y + $O + $H;
		
		#(fragM + (charge*proton_mass))/float(charge)
		$mz{"b+"}{$n} += $proton;
		$mz{"y+"}{$n} += $H + $proton;	
		
		if ($precursor_charge >=3){
			$mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
			$mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
		}
		if ($precursor_charge >=5){
			#$mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
			$mz{"b+++"}{$n} = ($mz{"b+"}{$n} + $proton )/3 ;
			#$mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
			$mz{"y+++"}{$n} = ($mz{"y+"}{$n} + $proton )/3;
		}
		# Isotopic ions 
		# mz <1010 M+1 M+2 M+3
		# mz >=1010 M+1 M+2 M+3 M+4
		# delta mz =~ 1.0	
		my $delta_mz = $massdiff_C12_C13;
		
		# if($mz{b}{$n} >= 1010){
		# $mz{b_m1}{$n} = $mz{b}{$n} + $delta_mz;
		# $mz{b_m2}{$n} = $mz{b}{$n} + 2*$delta_mz;
		# $mz{b_m3}{$n} = $mz{b}{$n} + 3*$delta_mz;
		# $mz{b_m4}{$n} = $mz{b}{$n} + 3*$delta_mz;
		# }else{
		# $mz{b_m1}{$n} = $mz{b}{$n} + 1;
		# }
		# if($mz{y}{$n} >= 1010){
		# $mz{y_m1}{$n} = $mz{y}{$n} + 1;
		# $mz{y_m2}{$n} = $mz{y}{$n} + 2;
		# }else{
		# $mz{y_m1}{$n} = $mz{y}{$n} + 1;
		# }
		
		$mz{"b_H2O+"}{$n} = $mz{"b+"}{$n} - $H2O;
		$mz{"b_NH3+"}{$n} = $mz{"b+"}{$n} - $NH3;
		$mz{"y_H2O+"}{$n} = $mz{"y+"}{$n} - $H2O;
		$mz{"y_NH3+"}{$n} = $mz{"y+"}{$n} - $NH3;
		
		if($pep =~ /S\[166.9984\]|T\[181.0140\]/){
			$mz{"b_H3PO+"}{$n} = $mz{"b+"}{$n} - $Phos_loss if $b_pep =~ /[ST]_/;
			$mz{"y_H3PO4+"}{$n} = $mz{"y+"}{$n} - $Phos_loss if $y_pep =~ /[ST]_/;
		}
		#print join "\t",$b_pep,$y_pep,"\n";
		
	}
	#print Dumper \%mz;
	my %mzs;
	foreach my $ion_type(keys %mz){
		foreach my $n (keys %{$mz{$ion_type}}){
			$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type;
		}
	}
	return \%mzs;
}


sub fragmentation_mod {
	
	#print Dumper \@_;
	my $sequence = shift;
	my $precursor_charge = shift;
	#my $shift_pos = shift;
	my $mass_shift = shift;
	my $shift_pos = shift;
	my $mod_name = shift;
	my $neutral_loss = shift;
	my $pepLen = shift;

	my $lossadd = 1;
	
	#print "mass_shift: $mass_shift, neutral_loss:$neutral_loss","\n";
	
	if($shift_pos =~ "N-"){
		$shift_pos = 1;
	}elsif($shift_pos =~ "C-"){
		$shift_pos = length($pepLen);
	}
	#$pep,$assumed_charge,$mod_mass,$sites_candidate{$site_mod_test}
	#my $sequence = "GHRPLDK";
	my $phostag = 0;
	if($sequence =~ /S\[166.9984\]|T\[181.0140\]/){
		$phostag = 1;
	}
	#$sequence =~ s/M\[147.0354\]/M_/g;
	#$sequence =~ s/C\[160.0307\]/C_/g;
	#$sequence =~ s/S\[166.9984\]/S_/g;
	#$sequence =~ s/T\[181.0140\]/T_/g;
	#$sequence =~ s/Y\[243.0297\]/Y_/g;
	#$sequence =~ s/(?<=[^MC])\[[^][]+\]//g;
	
	while($sequence =~ /(.)\[([^][]+)\]/g){
		my $aa = $1;
		my $massdiff = $2;
		$mono_aa{${aa}."_"} = $massdiff;
		$sequence =~ s/$aa\[[^][]+\]/${aa}_/;
	}
	#print $sequence,"\n";
	my %mz;
	my $mz_b;
	my $mz_y;
	my @seq_arr = split /(?=[[:upper:]][_#]?)/,$sequence;
	#print Dumper \@seq_arr;
	my $n;
	my $y_pep;
	my $b_pep;
	foreach my $index (0..$#seq_arr){
		$b_pep .= $seq_arr[$index];
		
		$mz_b += $mono_aa{$seq_arr[$index]};
		if($index + 1 == $shift_pos){
			$mz_b += $mass_shift;
			#$b_pep .= "[mod]";
			$b_pep .= "[mod.".$seq_arr[$index]."]";
		}
		$mz{"b+"}{++$n} = $mz_b;
		
		$y_pep .= $seq_arr[$#seq_arr - $index];
		$mz_y += $mono_aa{$seq_arr[$#seq_arr - $index]};
		if($n + $shift_pos == $#seq_arr + 2){
			$mz_y += $mass_shift;
			#$y_pep .= "mod";
			$y_pep = "[mod.".$seq_arr[$#seq_arr + 1 - $n]."]" . $y_pep ;
		}
		$mz{"y+"}{$n} = $mz_y + $O + $H;
		
		#(fragM + (charge*proton_mass))/float(charge)
		$mz{"b+"}{$n} += $proton;
		$mz{"y+"}{$n} += $H + $proton;	
		
		
		# if ($precursor_charge >=4){
			# $mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
			# $mz{"b+++"}{$n} = ($mz{"b+"}{$n} + $proton )/3 ;
			# $mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
			# $mz{"y+++"}{$n} = ($mz{"y+"}{$n} + $proton )/3;
		# }
		# Isotopic ions 
		# mz <1010 M+1 M+2 M+3
		# mz >=1010 M+1 M+2 M+3 M+4
		# delta mz =~ 1.0	
		my $delta_mz = $massdiff_C12_C13;
		
		# if($mz{b}{$n} >= 1010){
		# $mz{b_m1}{$n} = $mz{b}{$n} + $delta_mz;
		# $mz{b_m2}{$n} = $mz{b}{$n} + 2*$delta_mz;
		# $mz{b_m3}{$n} = $mz{b}{$n} + 3*$delta_mz;
		# $mz{b_m4}{$n} = $mz{b}{$n} + 3*$delta_mz;
		# }else{
		# $mz{b_m1}{$n} = $mz{b}{$n} + 1;
		# }
		# if($mz{y}{$n} >= 1010){
		# $mz{y_m1}{$n} = $mz{y}{$n} + 1;
		# $mz{y_m2}{$n} = $mz{y}{$n} + 2;
		# }else{
		# $mz{y_m1}{$n} = $mz{y}{$n} + 1;
		# }
		
		# $mz{"b_H2O+"}{$n} = $mz{"b+"}{$n} - $H2O;
		# $mz{"b_NH3+"}{$n} = $mz{"b+"}{$n} - $NH3;
		# $mz{"y_H2O+"}{$n} = $mz{"y+"}{$n} - $H2O;
		# $mz{"y_NH3+"}{$n} = $mz{"y+"}{$n} - $NH3;
		
		# $mz{"b_H2O+"}{$n} = $mz{"b+"}{$n} - $H2O if $b_pep =~ /[DEST]/ and $b_pep !=~ /mod.[ST]/ and $lossadd;
		# $mz{"b_NH3+"}{$n} = $mz{"b+"}{$n} - $NH3 if $b_pep =~ /[KNQR]/ and $b_pep !=~ /mod.[ST]/ and $lossadd;
		# $mz{"y_H2O+"}{$n} = $mz{"y+"}{$n} - $H2O if $y_pep =~ /[DEST]/ and $y_pep !=~ /mod.[ST]/ and $lossadd;
		# $mz{"y_NH3+"}{$n} = $mz{"y+"}{$n} - $NH3 if $y_pep =~ /[KNQR]/ and $y_pep !=~ /mod.[ST]/ and $lossadd;
		
		$mz{"b_NL_+"}{$n} = $mz{"b+"}{$n} - $neutral_loss if $neutral_loss != 0 && $neutral_loss != $mass_shift && $b_pep =~ /mod/;
		$mz{"y_NL_+"}{$n} = $mz{"y+"}{$n} - $neutral_loss if $neutral_loss != 0 && $neutral_loss != $mass_shift && $y_pep =~ /mod/;
		
		# $mz{"b_PhosNL_+"}{$n} = $mz{"b+"}{$n} - $Phos_loss if $phostag && $b_pep =~ /[ST]_/; 
		# $mz{"y_PhosNL_+"}{$n} = $mz{"y+"}{$n} - $Phos_loss if $phostag && $y_pep =~ /[ST]_/;
		
		if ($precursor_charge >=3){
			$mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
			$mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
			$mz{"b_NL_++"}{$n} = ($mz{"b_NL_+"}{$n} + $proton )/2 if $neutral_loss != 0 && $neutral_loss != $mass_shift && $b_pep =~ /mod/;
			$mz{"y_NL_++"}{$n} = ($mz{"y_NL_+"}{$n} + $proton )/2 if $neutral_loss != 0 && $neutral_loss != $mass_shift && $y_pep =~ /mod/;
			
			# $mz{"b_PhosNL_++"}{$n} = ($mz{"b_PhosNL_+"}{$n} + $proton )/2 if $phostag && $b_pep =~ /[ST]_/; 
			# $mz{"y_PhosNL_++"}{$n} = ($mz{"y_PhosNL_+"}{$n} + $proton )/2 if $phostag && $y_pep =~ /[ST]_/;
		
		}
		if ($precursor_charge == 2 && ($y_pep =~ /H/)){
			$mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
		}
		#print join "\t",$sequence,$b_pep,$y_pep,"\n";
		
	}
	#print STDERR Dumper $sequence,\%mz;#exit;
	my %mzs;
	foreach my $ion_type(keys %mz){
		foreach my $n (keys %{$mz{$ion_type}}){
			if($ion_type =~ /y/){
				if($n >= $pepLen-$shift_pos + 1){
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type."[mod]";
				}else{
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type;
				}	
			}
			if($ion_type =~ /b/){
				if($n >= $shift_pos){
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type."[mod]";
				}else{
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type;
				}
			}
		}
	}
	return \%mzs;
}


sub fragmentation_phosSite {
	
	#print Dumper \@_;
	my $sequence = shift;
	my $precursor_charge = shift;
	#my $shift_pos = shift;
	my $mass_shift = 79.9663;
	my $shift_pos = shift;
	#my $mod_name = shift;
	my $neutral_loss = $Phos_loss;
	my $pepLen = shift;
	my $lossadd = 1;
	
	my %mod_pos;
	map{$mod_pos{$_} = 1}@{$shift_pos};
	#print "mass_shift: $mass_shift, neutral_loss:$neutral_loss","\n";

	while($sequence =~ /(.)\[([^][]+)\]/g){
		my $aa = $1;
		my $massdiff = $2;
		$mono_aa{${aa}."_"} = $massdiff;
		$sequence =~ s/$aa\[[^][]+\]/${aa}_/;
	}
	#print $sequence,"\n";
	my %mz;
	my $mz_b;
	my $mz_y;
	my @seq_arr = split /(?=[[:upper:]][_#]?)/,$sequence;
	#print Dumper \@seq_arr;
	my $n;
	my $y_pep;
	my $b_pep;
	foreach my $index (0..$#seq_arr){
		$b_pep .= $seq_arr[$index];
		
		$mz_b += $mono_aa{$seq_arr[$index]};
		#if($index + 1 == $shift_pos){
		if($mod_pos{$index + 1}){
			$mz_b += $mass_shift;
			$b_pep .= "[mod.".$seq_arr[$index]."]";
		}
		$mz{"b+"}{++$n} = $mz_b;
		
		$y_pep = $seq_arr[$#seq_arr - $index] . $y_pep;
		$mz_y += $mono_aa{$seq_arr[$#seq_arr - $index]};
		#if($n + $shift_pos == $#seq_arr + 2){
		if($mod_pos{ $#seq_arr + 2 - $n }){
			$mz_y += $mass_shift;
			$y_pep = "[mod.".$seq_arr[$#seq_arr + 1 - $n]."]" . $y_pep ;
		}
		$mz{"y+"}{$n} = $mz_y + $O + $H;
		
		#(fragM + (charge*proton_mass))/float(charge)
		$mz{"b+"}{$n} += $proton;
		$mz{"y+"}{$n} += $H + $proton;	
		
		
		# if ($precursor_charge >=4){
			# $mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
			# $mz{"b+++"}{$n} = ($mz{"b+"}{$n} + $proton )/3 ;
			# $mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
			# $mz{"y+++"}{$n} = ($mz{"y+"}{$n} + $proton )/3;
		# }
		# Isotopic ions 
		# mz <1010 M+1 M+2 M+3
		# mz >=1010 M+1 M+2 M+3 M+4
		# delta mz =~ 1.0	
		my $delta_mz = $massdiff_C12_C13;
		
		# if($mz{b}{$n} >= 1010){
		# $mz{b_m1}{$n} = $mz{b}{$n} + $delta_mz;
		# $mz{b_m2}{$n} = $mz{b}{$n} + 2*$delta_mz;
		# $mz{b_m3}{$n} = $mz{b}{$n} + 3*$delta_mz;
		# $mz{b_m4}{$n} = $mz{b}{$n} + 3*$delta_mz;
		# }else{
		# $mz{b_m1}{$n} = $mz{b}{$n} + 1;
		# }
		# if($mz{y}{$n} >= 1010){
		# $mz{y_m1}{$n} = $mz{y}{$n} + 1;
		# $mz{y_m2}{$n} = $mz{y}{$n} + 2;
		# }else{
		# $mz{y_m1}{$n} = $mz{y}{$n} + 1;
		# }
		
		# $mz{"b_H2O+"}{$n} = $mz{"b+"}{$n} - $H2O if $b_pep =~ /[DEST]/ and $b_pep !=~ /mod.[ST]/ and $lossadd;
		# $mz{"b_NH3+"}{$n} = $mz{"b+"}{$n} - $NH3 if $b_pep =~ /[KNQR]/ and $b_pep !=~ /mod.[ST]/ and $lossadd;
		# $mz{"y_H2O+"}{$n} = $mz{"y+"}{$n} - $H2O if $y_pep =~ /[DEST]/ and $y_pep !=~ /mod.[ST]/ and $lossadd;
		# $mz{"y_NH3+"}{$n} = $mz{"y+"}{$n} - $NH3 if $y_pep =~ /[KNQR]/ and $y_pep !=~ /mod.[ST]/ and $lossadd;
		
		#$mz{"b_NL_+"}{$n} = $mz{"b+"}{$n} - $neutral_loss if $b_pep =~ /mod.[ST]/;
		#$mz{"y_NL_+"}{$n} = $mz{"y+"}{$n} - $neutral_loss if $y_pep =~ /mod.[ST]/;
		
		if ($precursor_charge >=3){
			$mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
			$mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
			# $mz{"b_NL_++"}{$n} = ($mz{"b_NL_+"}{$n} + $proton )/2 if $b_pep =~ /mod.[ST]/;
			# $mz{"y_NL_++"}{$n} = ($mz{"y_NL_+"}{$n} + $proton )/2 if $y_pep =~ /mod.[ST]/;
		}
		
		if ($precursor_charge == 2 && ($y_pep =~ /H/)){
			$mz{"y++"}{$n} = ($mz{"y+"}{$n} + $proton )/2;
		}
		# if ($precursor_charge == 2 && ($b_pep =~ /[KR]/)){
			# $mz{"b++"}{$n} = ($mz{"b+"}{$n} + $proton )/2;
		# }
		
	}
	
	#print STDERR Dumper $sequence,\%mz;#exit;
	my %mzs;
	my $min_pos = min(@{$shift_pos});
	my $max_pos = max(@{$shift_pos});
	foreach my $ion_type(keys %mz){
		foreach my $n (keys %{$mz{$ion_type}}){
			if($ion_type =~ /y/){
				if($n >= $pepLen - $max_pos + 1){
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type."[mod]";
				}else{
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type;
				}	
			}
			if($ion_type =~ /b/){
				if($n >= $min_pos){
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type."[mod]";
				}else{
					$mzs{$mz{$ion_type}{$n}} = $n.".".$ion_type;
				}
			}
		}
	}
	#print STDERR Dumper $sequence,$shift_pos,\%mzs;exit;
	return \%mzs;
}


sub bin_mode_fragmentation{
	my $sequence  = shift;
	my $wd = shift;
	#my $wd = 0.2;
	my $step = shift;
	#my $step = 1;
	my $precursor_charge = shift;	
	my $scale = (10 ** $step);
	my $step_ = 1 / $scale;
	#my $step = 0.01;
	#$sequence = "GHRPLDK";
	my $fragments = fragmentation($sequence,$precursor_charge);
	
	#print Dumper $fragments;
	my $bin_data;
	foreach my $ion_type(keys %{$fragments}){
		foreach my $n (keys %{$fragments->{$ion_type}}){
			my $mz = $fragments->{$ion_type}->{$n};
			my ($start , $end) = (sprintf("%.${step}f",($mz - $wd)) , sprintf("%.${step}f",($mz + $wd)));
			my $i;
			for($i = $start; $i <= $end; $i += $step_){
				$bin_data->{sprintf("%.${step}f",$i)} = [$mz,$ion_type]; # 20180621
			}
		}
	}
	#print Dumper::Sortedkeys $bin_data;
	#print scalar (keys %{$bin_data}),"\n";
	return $bin_data;
}

sub calmz{
	my ($peptide,$charge) = @_;
	
	#print STDERR $peptide,"\n";
	
	while($peptide =~ /(.)\[([^][]+)\]/g){
		my $aa = $1;
		my $massdiff = $2;
		#print $aa,"\t",$massdiff,"\n";
		$mono_aa{${aa}."_"} = $massdiff;
		$peptide =~ s/S\[166.9984\]/S_/g;
		$peptide =~ s/T\[181.0140\]/T_/g;
		$peptide =~ s/Y\[243.0297\]/Y_/g;		
		$peptide =~ s/$aa\[[^][]+\]/${aa}_/;
	}
	
	my $monoMass;
	#print STDERR $peptide,"\n";
	
	my %AAs;
	my @AAs = split /(?=[[:upper:]]_?)/,$peptide;
	#print STDERR Dumper \@AAs;
	foreach my $aa(@AAs){
		$monoMass += $mono_aa{$aa};
	}
	$monoMass += $O + $H + $proton;
	
	my $calmz = ($monoMass + $charge * $proton)/$charge;
	return $calmz;
}


1;
