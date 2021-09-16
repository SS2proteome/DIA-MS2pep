#!/usr/bin/env perl

use strict;
use Data::Dumper;
use List::Util qw(max min sum);
use POSIX qw(floor ceil round);
use File::Basename;
use Parallel::ForkManager;
use lib dirname(__FILE__);
use peptide;
use SWATH_pseudo_ms2_InParallel;

$| = 1;
my $mgf_file = "$ARGV[0].mgf";
my $ms1_file = "$ARGV[0].ms1";
my $mzML_file = "$ARGV[0].mzML";
#my $window_file = $ARGV[1];
my $charge_dir = $ARGV[1];
my $ms2_precison = $ARGV[2];
my $max_processes = $ARGV[3];
my $num_of_data_points;
my $logfile = "pseudo.${mgf_file}.log";

open(LOGF,">$charge_dir/$logfile") or die "$charge_dir/$logfile $!\n";

my $window_file = "$ARGV[0].mzML.SWATH_acquisition_window.txt";

my $min_peaks_splited = 20;


my %ms2toms1;
my $start_premz ;
my $end_premzmz ;
my $windowsize ;
my $ms2_low_mz = 140;
my $ms2_min_intensity = 2;

my $ms2scan;
my %SWATH_window;
my %windows_index;

=head window_file
Window #	Q1 Start (Da)	Q1 Stop (Da)	Collision Energy Spread (V)	Mass Window Center	Total Q1 Window Size (Da)
1	399.5	408.2	10		403.85	8.7
2	407.2	415.8	10		411.5	8.6
3	414.8	422.7	10		418.75	7.9
4	421.7	429.7	10		425.7	8
5	428.7	437.3	10		433	8.6
6	436.3	444.8	10		440.55	8.5
7	443.8	451.7	10		447.75	7.9
=cut 
my %isolationlist;
my $number_of_windows;
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
	$windows_index{$line[0]} = $line[0];
	$number_of_windows = $line[0];
	$SWATH_window{$line[-2]}{window} = $line[0];
}
close(FH);
#print Dumper \%SWATH_window;exit;


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

# READ mzML file
print "*****Read mzML file: $mzML_file!\n";
open(FH,$mzML_file) or die "$! $mzML_file\n";
while(<FH>){
	if(m{<spectrum index="(\d+)" id="sample=[^ ]+ period=[^ ]+ cycle=(\d+) experiment=(\d+)"}){
		($index, $cycle,$experiment) = ($1,$2,$3-1);
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="([^"]+)"}){
		$scantime = $1;
		
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="([^"]+)" }){
		if($mslevel == 2){
			$basepeak_ms1->{$premz_target} = int($1);
		}
	}elsif(/name="ms level" value="(\d+)"/){
		$mslevel = $1;
		if($mslevel == 2){
			# $ms2scan = $experiment; 
			$ms2scan++; 
			#$index_ms2scan{$ms2scan} = $index_ms1;
			if($lastbasepeak_ms1->{$premz_target} == $basepeak_ms1->{$premz_target}){
				$basepeak_ms1_cnt->{$premz_target}++;
			}else{
				push @basepeak_ms1_hash,$basepeak_ms1_cnt->{$premz_target};
				#print STDERR  $scantime,"\t",$basepeak_ms1,"\t",$basepeak_ms1_cnt,"\n";
				$basepeak_ms1_cnt->{$premz_target} = 1;
			}
			$lastbasepeak_ms1->{$premz_target} = $basepeak_ms1->{$premz_target};
		}else{
			#$index_ms1++;
			$ms1scan++;
			$cycletime += $scantime - $lastscantime;
			$lastscantime = $scantime;
			
			$cyclecount++;
			#$ms1scan = $cycle;
			#$index_ms1scan{$index_ms1} = sprintf("%06.f",$ms1scan);
		}
	}elsif(m{accession="MS:1000827" name="isolation window target m/z" value="([\d.]+)" }){
		my $premz_ = $1;
		$ms2toms1{$ms2scan} = [sprintf("%06.f",$ms1scan),$premz_];
		#$ms2toms1{$index} = [sprintf("%06.f",$ms1scan),$premz_];
		#print join "#",$index,$ms1scan,$premz_,"\n";
		#$ms1toms2{sprintf("%06.f",$ms1scan)}{$ms2scan} = 1;
		printf "\r$ms1scan";
	}elsif(m{<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="([^"]+)"} and $mslevel == 2){
		$basepeak_intensity{$ms2scan} = $1;
	}
}
close(FH);
#print Dumper \%ms2toms1;exit;
print "\n";
my @basepeak_ms1_hash_2 = grep{$_ > 1}@basepeak_ms1_hash;

print LOGF join "\t",$ARGV[0],$number_of_scan,$number_of_windows,$cycletime/$cyclecount,$cyclecount, $scantime,ceil(sum(@basepeak_ms1_hash_2)/($#basepeak_ms1_hash_2 + 1))-1,"\n";# exit;
print join "\t",$ARGV[0],$number_of_scan,$number_of_windows,$cycletime/$cyclecount,$cyclecount, $scantime,ceil(sum(@basepeak_ms1_hash_2)/($#basepeak_ms1_hash_2 + 1))-1,"\n";# exit;

$num_of_data_points = ceil(sum(@basepeak_ms1_hash_2)/($#basepeak_ms1_hash_2 + 1))-1;
$num_of_data_points = $num_of_data_points < 3 ? 3 : $num_of_data_points;
$num_of_data_points = $num_of_data_points > 6 ? 6 : $num_of_data_points;
print LOGF "num_of_data_points: $num_of_data_points\n";


close(LOGF);


my %filehandle;
#my %filehandle_num;
mkdir "${charge_dir}/$ARGV[0]";
foreach my $key (1..$number_of_windows){
	open($filehandle{$key},">${charge_dir}/$ARGV[0]/index${key}_$mgf_file") or die "Can't write mgf file! ${charge_dir}/$ARGV[0]/index${key}_$mgf_file $!\n";
}


my $line_block;
my $premz;
open(FH,$mgf_file) or die "$! $mgf_file \n";
while(<FH>){
	
	chomp;
	my $Line = $_;
	#printf "\r$Line" if $Line =~ /TITLE/;
	$line_block .= $_."\n";
	if (/^PEPMASS=([\d.]+)/){
		$premz = $1;

	}
	if(/END/){
		my $fh =  $filehandle{$SWATH_window{$premz}{window}};
		print $fh $line_block;
		$line_block = "";
	}
	
}
foreach my $key (1..$number_of_windows){
		close($filehandle{$key});
}
close(FH);


print "\nMGF file spliting!\n";

my $pm = new Parallel::ForkManager($max_processes);  

my $i=1;

for($i=1;$i<=$number_of_windows;$i++){
	
	$pm->start and next;

	#print "number_of_windows: ",$i,"\n";
	my $index = $i;
	my $mgf_split_file = "index${index}_$mgf_file";
	my %thread_list;
	pseudo_ms2($mgf_split_file,"$charge_dir/@ARGV[0]",$ms2_precison,$num_of_data_points,@ARGV[0]);

	$pm->finish;
}

$pm->wait_all_children;

my @w_mgf_files = glob "$charge_dir/@ARGV[0]/w*@ARGV[0].mgf";
#print Dumper \@w_mgf_files;


my %filehandle;
#my %filehandle_num;
foreach my $key (keys %isolationlist){

	open($filehandle{$key},">${charge_dir}/w_${key}_$mgf_file") or die "Can't write mgf file! ${charge_dir}/w_${key}_$mgf_file $!\n";

}

foreach my $file (@w_mgf_files){
	my $windowsize;
	$file =~ m{.+/w_(\d+)_index\d+_.+};
	$windowsize = $1;
	#print join "\t",$file,$windowsize,"\n";
	my $fh_ = $filehandle{$windowsize};
	open(FH,$file);
	while(<FH>){
		print $fh_ $_;
	}
	close(FH);

}

foreach my $key (keys %isolationlist){
	#foreach my $z ($min_charge .. $max_charge){
		close($filehandle{$key});
		#close($filehandle{$key}{$z});
	#}
}

my @files = glob "$charge_dir/@ARGV[0]/*";
unlink $_ for @files;

rmdir "$charge_dir/@ARGV[0]";

#	print "****\n";

sub pseudo_ms2 {
	
	SWATH_pseudo_ms2_InParallel::demultiplexing(@_);
	#print "\nDemultiplexing Finished ..\n";
}






