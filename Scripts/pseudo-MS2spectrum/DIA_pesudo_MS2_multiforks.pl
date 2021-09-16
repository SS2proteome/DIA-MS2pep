#!/usr/bin/env perl

use strict;
use Data::Dumper;
use POSIX qw(floor);
use List::Util qw(max min sum);
use POSIX qw(floor ceil round);
use Parallel::ForkManager;
use File::Basename;
use lib dirname(__FILE__);
use pseudo_ms2_InParallel;

$| = 1;

my $mgf_file = "$ARGV[0].mgf";
my $mzML_file = "$ARGV[0].mzML";
my $charge_dir = $ARGV[1];
my $ms2_precison = $ARGV[2];# ppm;
my $max_processes = $ARGV[3];
my $num_of_data_points;# = $ARGV[3]; #@ARGV[5];
my $logfile = "pseudo.${mgf_file}.log";
open(LOGF,">$charge_dir/$logfile") or die "$charge_dir/$logfile $!\n";

my $window_file = "$ARGV[0].mzML.DIA_acquisition_window.txt";

my $ms2scan;
my %DIA_window;
my %isolationlist;
my $number_of_windows;
my %windows_index;
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
	$windows_index{$line[0]} = $line[0];
	$number_of_windows = $line[0];
	$DIA_window{$line[-2]}{window} = $line[0];
	#print $_,"\n";
}
close(FH);
#print $number_of_windows,"\n";
#print STDERR Dumper \%DIA_window; #exit;

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
			push @total_ms2scan,$scan;
			$ms2scan = $scan; 
			if($lastbasepeak_ms1->{$premz_target} == $basepeak_ms1->{$premz_target}){
				$basepeak_ms1_cnt->{$premz_target}++;
			}else{
				push @basepeak_ms1_hash,$basepeak_ms1_cnt->{$premz_target};
				#print STDERR  $scantime,"\t",$basepeak_ms1,"\t",$basepeak_ms1_cnt,"\n";
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
		$ms2toms1{$ms2scan} = [$ms1scan,$premz_target];
		printf "\r$ms1scan";
	}
	elsif(m{<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="([^"]+)"}){
		$basepeak_intensity_whole{$ms2scan} = $1;
	}
}
close(FH);

# my ($s,$n);
# foreach my $i (@basepeak_ms1_hash){
	
	# $s->{$i}++;
	# $n += $i;

# }
print "\n";
my @sort_basepeak_ms1_hash = sort {$a <=> $b}@basepeak_ms1_hash;
print LOGF join "\t",$ARGV[0],$number_of_scan,$number_of_windows,$cycletime/$cyclecount,$cyclecount, $scantime,	ceil(sum(@sort_basepeak_ms1_hash)/($#sort_basepeak_ms1_hash + 1)), @sort_basepeak_ms1_hash[int($#sort_basepeak_ms1_hash / 2) ],	"\n";
#print  join "\t",$ARGV[0],$number_of_scan,$number_of_windows,$cycletime/$cyclecount,$cyclecount, $scantime,	ceil(sum(@sort_basepeak_ms1_hash)/($#sort_basepeak_ms1_hash + 1)), @sort_basepeak_ms1_hash[int($#sort_basepeak_ms1_hash / 2) ],	"\n";

$num_of_data_points = ceil(sum(@sort_basepeak_ms1_hash)/($#sort_basepeak_ms1_hash + 1));
#$num_of_data_points = @sort_basepeak_ms1_hash[int($#sort_basepeak_ms1_hash / 2) ];
$num_of_data_points = $num_of_data_points < 3 ? 3 : $num_of_data_points;
$num_of_data_points = $num_of_data_points > 6 ? 6 : $num_of_data_points;
#print   "num_of_data_points: $num_of_data_points\n";
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
	printf "\r$Line" if $Line =~ /TITLE/;
	$line_block .= $_."\n";
	if (/^PEPMASS=([\d.]+)/){
		$premz = $1;

	}
	if(/END/){
		my $fh =  $filehandle{$DIA_window{$premz}{window}};
		#print STDERR join "\t",$premz,$DIA_window{$premz}{window},$fh,"\n";
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
	
	pseudo_ms2_InParallel::demultiplexing(@_);
	#print "\nDemultiplexing Finished ..\n";
}


