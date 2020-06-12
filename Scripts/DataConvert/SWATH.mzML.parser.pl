#!/usr/bin/env perl

use strict;
use Data::Dumper;
use MIME::Base64;
use MIME::Base32; 
use IO::Uncompress::Inflate qw(inflate $InflateError) ;

my $ms2_min_intensity = 1;
my $ms1_min_intensity = 3;

# usage perl script [mgf/ms1] [file]

if($ARGV[0] eq "mgf"){
	Wiff_mzMLtomgf(@ARGV[1]);
}
if($ARGV[0] eq "ms1"){
	Wiff_mzMLtoms1(@ARGV[1]);
}
#$mzMLfile = "Site10_DynRge_S1_Day_1.mzML";

#my $ref_mz_inten = base64decode(\$mz_binray);
#my $ref_count = base32decode(\$count_binary);
sub Wiff_mzMLtomgf{
	my $scan_count;
	my $run_id;
	my $spectrum_index;
	my $cycle;
	my $experiment;
	my $mslevel;
	my $retention_time;
	my $isolationwindow;
	my $decode_type;
	my $mz_binray;
	my $mz;
	my $count_binary;
	my $count;
	my $precusor_mz;
	my $title;
	my $mzMLfile = shift;
	#my $remove_scan_window_tag = shift;
	open(FH,$mzMLfile) or die "$!\n $mzMLfile";
	my $mgfname = $mzMLfile;
	$mgfname =~ s/mzML/mgf/;
	if($mgfname =~ /mgf/){
		open(OUT,">".$mgfname) or die "$!$mgfname";
	}
	while(<FH>){
		chomp;
		if(m{<run id="([^"]+)"}){
			$run_id = $1;
		}elsif(m{<spectrum index="(\d+)" id="sample=[^ ]+ period=[^ ]+ cycle=(\d+) experiment=(\d+)"}){
			$scan_count++;
			printf STDERR "\r$scan_count";
			($spectrum_index,$cycle,$experiment) = ($1,$2,$3-1);
			$title = join ".",$run_id,$spectrum_index,$cycle,$experiment;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="(\d+)"/>}){
			$mslevel = $1;
		} 
		if($mslevel == 1){
			next;
		}else{
		if(m{<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="([^"]+)" }){
			$retention_time = sprintf ("%.3f",$1 * 60); # seconds
		 }elsif(m{<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="([^"]+)" }){
			 $isolationwindow->{target} = $1;
			 $precusor_mz = $1;
		# }elsif(m{<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="([^"]+)"}){
			# $isolationwindow->{lower} = $1;
		# }elsif(m{<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="([^"]+)"}){
			# $isolationwindow->{upper} = $1;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000523" name="(\d+)-bit float"}){
			$decode_type = $1;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000514" name="m/z array"}){
			my $nextline = <FH>;
			$nextline =~ m{<binary>([^<>]+)</binary>};
			$mz_binray = $1;
			$mz = $decode_type == "64" ? base64decode(\$mz_binray) : base32decode(\$mz_binray);
			
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000521" name="(\d+)-bit float"}){
			$decode_type = $1;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000515" name="intensity array"}){
			my $nextline = <FH>;
			$nextline =~ m{<binary>([^<>]+)</binary>};
			$count_binary = $1;
			$count = $decode_type == "64" ? base64decode(\$count_binary) : base32decode(\$count_binary);
		}elsif(m{</binaryDataArrayList>}){
			if($mslevel == 2){
				mgf_print(	$title,
						$retention_time,
						$precusor_mz,
						$mz,
						$count
					);
				undef $mz;
				undef $count;
			}
		}
		}
	}
	close(FH);
	if($mgfname =~ /mgf/){close(OUT)};
	print "\nDone!\n";
}

sub mgf_print{
	my ($title,$retention_time,$precusor_mz,$mz,$count) = (@_);
	print OUT "BEGIN IONS\n";
	print OUT "TITLE=".$title,"\n";
	print OUT "RTINSECONDS=".$retention_time,"\n";
	print OUT "PEPMASS=".$precusor_mz,"\n";
	#print OUT "CHARGE=".$_."+","\n";
	my $len = $#$mz;
	foreach my $index (0..$len){
		print OUT $mz->[$index]," ",$count->[$index],"\n" if ($count->[$index] > $ms2_min_intensity);
	}
	print OUT "END IONS\n";	
}

sub Wiff_mzMLtoms1{
	
	my $scan_count;
	my $run_id;
	my $spectrum_index;
	my $cycle;
	my $experiment;
	my $mslevel;
	my $retention_time;
	my $isolationwindow;
	my $decode_type;
	my $mz_binray;
	my $mz;
	my $count_binary;
	my $count;
	my $precusor_mz;
	my $title;
	my $mzMLfile = shift;
	#my $remove_scan_window_tag = shift;
	open(FH,$mzMLfile) or die "$!\n $mzMLfile";
	my $ms1name = $mzMLfile;
	$ms1name =~ s/mzML/ms1/;
	if($ms1name =~ /ms1/){
		open(OUT,">".$ms1name) or die "$!$ms1name";
	}
	while(<FH>){
		chomp;
		if(m{<run id="([^"]+)"}){
			$run_id = $1;
		}elsif(m{<spectrum index="(\d+)" id="sample=[^ ]+ period=[^ ]+ cycle=(\d+) experiment=(\d+)"}){
			$scan_count++;
			printf STDERR "\r$scan_count";
			($spectrum_index,$cycle,$experiment) = ($1,$2,$3-1);
			$title = join ".",$run_id,$spectrum_index,$cycle,$experiment;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="(\d+)"/>}){
			$mslevel = $1;
		}
		if($mslevel == 2){
			next;
		}else{
		if(m{<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="([^"]+)" }){
			$retention_time = sprintf ("%.3f",$1 * 60); # seconds
		# }elsif(m{<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="([^"]+)" }){
			# $isolationwindow->{target} = $1;
			# $precusor_mz = $1;
		# }elsif(m{<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="([^"]+)"}){
			# $isolationwindow->{lower} = $1;
		# }elsif(m{<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="([^"]+)"}){
			# $isolationwindow->{upper} = $1;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000523" name="(\d+)-bit float"}){
			$decode_type = $1;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000514" name="m/z array"}){
			my $nextline = <FH>;
			$nextline =~ m{<binary>([^<>]+)</binary>};
			$mz_binray = $1;
			$mz = $decode_type == "64" ? base64decode(\$mz_binray) : base32decode(\$mz_binray);
			
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000521" name="(\d+)-bit float"}){
			$decode_type = $1;
		}elsif(m{<cvParam cvRef="MS" accession="MS:1000515" name="intensity array"}){
			my $nextline = <FH>;
			$nextline =~ m{<binary>([^<>]+)</binary>};
			$count_binary = $1;
			$count = $decode_type == "64" ? base64decode(\$count_binary) : base32decode(\$count_binary);
		}elsif(m{</binaryDataArrayList>}){
			if($mslevel == 1){
				ms1_print($cycle,
						$retention_time,
						#$precusor_mz,
						$mz,
						$count
					);
				undef $mz;
				undef $count;
			}
		}
		}
	}
	close(FH);
	if($ms1name =~ /ms1/){close(OUT)};
	print "\nDone!\n";
}

sub ms1_print{
	my ($cycle,$retention_time,$mz,$count) = (@_);
	print OUT sprintf("%s\t%06.f\t%06.f\n","S",$cycle,$cycle);
	print OUT sprintf("%s\t%s\t%.3f\n","I","RTime",$retention_time);
	my $len = $#$mz;
	foreach my $index (0..$len){
		print OUT $mz->[$index]," ",$count->[$index],"\n" if $count->[$index] > $ms1_min_intensity ;
	}
}


sub base64decode{
	my $peaks = shift;
	my $base64_decoded = decode_base64(${$peaks});
	#print Dumper "base64_decoded",$base64_decoded;
	my $uncompress_;
	my $input = \$base64_decoded;
	inflate $input => \$uncompress_;
	#print Dumper "uncompress_",$uncompress_;
	my @mzs_intensities = unpack "d*",$uncompress_;
	#print Dumper "mzs_intensities",\@mzs_intensities;
	return \@mzs_intensities;
}

sub base32decode{
	my $peaks = shift;
	my $base64_decoded = decode_base64(${$peaks});
	my $uncompress_;
	my $input = \$base64_decoded;
	inflate $input => \$uncompress_;
	my @peaks = unpack("V*",$uncompress_);
	my @mzs_intensities;
	@mzs_intensities = map{
		sprintf("%.5f",unpack("f",pack("I",$_)));
	}(@peaks);
	#print Dumper "mzs_intensities",\@mzs_intensities;
	return \@mzs_intensities;
}

