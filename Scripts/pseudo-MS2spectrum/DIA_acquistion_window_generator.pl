#!/usr/bin/env perl

use strict;
use Data::Dumper;


=head file out format
Window #	Q1 Start (Da)	Q1 Stop (Da)	Collision Energy Spread (V)	Mass Window Center	Total Q1 Window Size (Da)	
1	399.5	408.2	10		403.85	8.7

=cut

my $mzMLfile = @ARGV[0];
my $outfile = "${mzMLfile}.DIA_acquisition_window.txt";
open(OUT,">".$outfile) or die "$! $outfile \n";
print OUT "Window #	Q1 Start (Da)	Q1 Stop (Da)	Collision Energy Spread (V)	Mass Window Center	Total Q1 Window Size (Da)\n";

my %isolationwindow;
my $mslevel;
my $scan;
my %tmp;
my $isolationwindow;
open(FH,$mzMLfile) or die "$!\n $mzMLfile";
while(<FH>){
	chomp;
	if(m{<spectrum index="[^"]+" id="controllerType=[^ ]+ controllerNumber=[^ ]+ scan=([^ "]+)}){
		$scan = $1;
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="(\d)"/>}){
		#my $tmp = $mslevel;
		($mslevel) = ($1);
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="([^"]+)" }){
		$isolationwindow->{target} = $1;
		last if (++$tmp{$1} == 2);
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="([^"]+)"}){
		$isolationwindow->{lower} = $1;
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="([^"]+)"}){
		$isolationwindow->{upper} = $1;
		if($mslevel == 2){
			print OUT join "\t",$scan-1,$isolationwindow->{target} - $isolationwindow->{lower},$isolationwindow->{target}+$isolationwindow->{upper},"","",$isolationwindow->{target},sprintf("%d",$isolationwindow->{upper} + $isolationwindow->{lower});
			print OUT "\n";
		}
	}
}
close(FH);
close(OUT);
