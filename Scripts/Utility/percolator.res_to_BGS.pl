#!/usr/bin/env perl

use strict;
use Data::Dumper;

my $filename = $ARGV[0]; # filename
my $mzml_dir = $ARGV[1]; 
my $ssl_file = $filename.".BGS.tsv";
my $mzmlfile = $mzml_dir."/".$filename.".mzML";
print $mzmlfile,"\n";

# file   scan    charge  sequence
# demo.ms2        8       3       VGAGAPVYLAAVLEYLAAEVLELAGNAAR
# demo.ms2        1806    2       LAESITIEQGK
# demo.ms2        2572    2       ELAEDGC[+57.0]SGVEVR
# demo.ms2        3088    2       TTAGAVEATSEITEGK
# demo.ms2        3266    2       DC[+57.0]EEVGADSNEGGEEEGEEC[+57.0]
# demo.ms2        9734    3       IWELEFPEEAADFQQQPVNAQ[-17.0]PQN
# demo.ms2        20919   3       VHINIVVIGHVDSGK
# ../elsewhere/spec.mzXML 00497   2       LKEPAQNTADNAK
# ../elsewhere/spec.mzXML 00680   2       ALEGPGPGEDAAHSENNPPR
# ../elsewhere/spec.mzXML 00965   2       FFSHEAEQK
# ../elsewhere/spec.mzXML 01114   2       C[+57.0]GPSQPLK
# ../elsewhere/spec.mzXML 01382   2       AVHVQVTDAEAGK

my $aa_mass = AA_mass();
my $rt_data;
open(PIN,$filename.".pin.target.pep.tsv") or die "no  $filename $!\n";
open(PINRAW,"$filename.pin") or die "no $filename.pin $!\n ";
while(<PINRAW>){
	chomp;
	s/\r//;
	next if $. == 1 && /PSMId/;
	my @line = split /\t/,$_;
	$rt_data->{$line[0]} = $line[25]/60;
}
close(PINRAW);

open(OUT,">".$ssl_file) or die "$ssl_file";
print OUT join "\t","Raw File","Stripped Sequence", "Precursor Charge", "Labeled Sequence", "Retention Time","Scan Number";
print OUT "\n";
while(<PIN>){
	chomp;
	s/\r//;
	next if $. == 1 && /PSMId/;
	my @line = split /\t/,$_;
	if($line[2] < 0.01){
		my @title = split /[.]/,$line[0];
		my $scan = $title[1];
		my $charge = $title[3];
		#my $file = $mzmlfile;
		my $sequence = $line[4];
		
		##if($filename eq "all"){
		my $file = $title[0];
		$file =~ s/raw_|res\d+_//g;
		$file = $file .".raw";
		#}
		my $ptm;
		if($sequence =~ /:\(/){
			next;
			#$sequence =~ s/....$//;
			$sequence =~ s/:\((.+)\)//;
			$ptm = $1;
		}
		
		$sequence =~ s/^..|..$//g;
		$sequence =~ s/^n|c$//g;
		my $stripped_seq = $sequence;
		$stripped_seq =~ s/\[[^][]+\]//g;
		my $mass;
		#$sequence =~ s/(.)(\[)([^][]+)(\])/${1}${2}($3>$aa_mass->{$1} ? sprintf("+%.01f",$3-$aa_mass->{$1}) : sprintf("-%.01f",$aa_mass->{$1}-$3))${4}/g;
		$sequence =~ s/(.)\[([^][]+)\](?{$mass = $2>$aa_mass->{$1} ? sprintf("+%.01f",$2-$aa_mass->{$1}) : sprintf("-%.01f",$aa_mass->{$1}-$2)})/$1\[$mass\]/g;
		$sequence =~ s/\[\+57.0\]/[CAM]/g;
		$sequence =~ s/\[\+16.0\]/[Oxi]/g;
		$sequence = "_".$sequence."_";
		#$sequence =~ s/(.)\[([^][]+)\]/$1\[$aa_mass->{$1}\]/g;
		
		# $sequence =~ s/C\[160.0307\]/C[+57.0]/g;
		# $sequence =~ s/M\[147.0354\]/M[+16.0]/g;
		# $sequence =~ s/S\[166.9984\]/S[+80.0]/g;
		# $sequence =~ s/T\[181.0140\]/T[+80.0]/g;
		# $sequence =~ s/Y\[243.0297\]/Y[+80.0]/g;

		
		
		print OUT join "\t",$file,$stripped_seq,$charge,$sequence,$rt_data->{$line[0]},$scan;
		print OUT "\n";
	}
}
close(PIN);
close(OUT);
print "Done: $ssl_file\n";


sub AA_mass{

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
	return \%mono_aa;
}
