#!/usr/bin/env perl

use strict;
use Data::Dumper;
$|=1;

my %pep_info;
my %pep_proteins;
my $peptidefile = shift;
my $outfile = "Protein_Group_$peptidefile";
my $fastafile = shift;
open(fh,$peptidefile) or die "$peptidefile$!\n";
while(<fh>){
	chomp;
	s/\r//;
	s/\[[^][]+\]//g;
	my @line = split /\t/,$_; # peptide in column 1
	$pep_info{$line[0]}++ if length($line[0]);
}


print "Read protein database...\n";
#my $fastafile = @ARGV[-1];
my $id;
my $pro_disc;
my $gene;
my $seq;
my %pep_pro;  # {pep => [pro1,pro2,...]}
my %database; # {UNIPROT-AC => [PROTEIN NAME, GENE NAME]}

print STDERR "ProteinRetrive\n";

ProteinRetrive($fastafile,\%pep_info);

#print Dumper \%database;

sub ProteinRetrive(){

	my $fastafile = shift;
	my $pep_info = shift;
	my $pep_proteins = shift;
	my $id;
	my $pro_disc;
	my $gene;
	my $seq;
	my @pep_info = keys %{$pep_info};  # {pep => [pro1,pro2,...]}
	#my %database; # {UNIPROT-AC => [PROTEIN NAME, GENE NAME]}
	my $line_;

	# use Parallel::ForkManager;
	# my $MAX_processes = 4;
	# my $pm = Parallel::ForkManager->new($MAX_processes);	
	my @sequences;
	
	open(fasta,"$fastafile") or die "no $fastafile $!\n";
	while(<fasta>){
	chomp;
	s/\r$//;
	my $line = $_;
	if(/^>/){
		#if($seq and $id !~ /[Cc]ontaminan/ and length $id and $pro_disc !~ /\(Fragment\)/ and length $gene) {
		if($seq and $id !~ /REV_/) {
			$database{$id} = [$pro_disc,$gene,length $seq];
			push @sequences, [$id,$seq];
		}
		$seq = "";
		if($line =~ /^>(REV_[^ ]+) Randomized Protein Sequence/){
			($id,$pro_disc) = ($1,"Randomized Protein Sequence");
			$gene = "REV_";
			
		}else{
			($id,$pro_disc) = (/[^|]+\|([^|]+)\|([^=]+) OS/);
			($id,$pro_disc) = ($line,$line) if ! defined $id;
			$gene = $_ =~ /GN=([^ ]+)/ ? " ".$1 : "";
		}
	}else{
		$seq .= $_;
	}
	}
	print "\n";
	close(fasta);
	$database{$id} = [$pro_disc,$gene,length $seq];
	push @sequences, [$id,$seq];
	
	# my %tmpfile;
	# my $step  = $#sequences%($MAX_processes) ? $#sequences/($MAX_processes) : int($#sequences/($MAX_processes))+1 ;
	# foreach my $fork (0..$MAX_processes-1){
		# $pm->start and next;
		# open($tmpfile{$fork},">".$fork.".tmp") or die "Can't write to ".$fork.".tmp $!\n";
		# my @subsequences = @sequences[$fork*$step..($fork == $MAX_processes-1 ? $#sequences : $fork*$step+$step-1)];
		foreach my $subseq (@sequences){
			print "\r".++$line_;
			map{
				if(index($subseq->[1],$_)>=0){
					push @{$pep_pro{$_}},$subseq->[0];
					# $tmpfile{$fork}->print  (join "\t", $_,$subseq->[0]);
					# $tmpfile{$fork}->print  ("\n");
				}
				
			}@pep_info;
		}
		#close($tmpfile{$fork});
		#$pm->finish;
	#}
	#$pm->wait_all_children;
	print "\n";

	# my @tempfiles = glob "*tmp";
	# #my %pep_pro;
	# foreach my $file (@tempfiles){#print $file,"\n";
		# open(FH,"$file") or die "Can't open $file $!\n";
		# while(<FH>){
			# chomp;
			# my @line = split;
			# #print $_,"\n";
			
		# }	
		# close(FH);
		# #unlink($file);
	# }
	#return (\%pep_pro,\%database);
}


my %pro_peps; # pro => [pep1,pep2,...]
my %peps_pros; # pep1#pep2#... => [pro1,pro2,...];
my %pros_cnt; # {pro1;pro2;...}++
my @cnt_pros_pep; # (pros_cnt,$pro1;$pro2;...,$pep)

print "\nOrganize dataset...\n";
foreach my $pep (keys %pep_pro){
	map{#print $pep,$_,"\n";
		push @{$pro_peps{$_}},$pep;
	}@{$pep_pro{$pep}};
}


foreach my $pro (keys %pro_peps){
	my $peps = join "#",@{$pro_peps{$pro}};
	push @{$peps_pros{$peps}},$pro;
}


foreach my $peps (keys %peps_pros){
	my $pros = join ";", @{$peps_pros{$peps}};
	map{
		$pros_cnt{$pros}++;
	}split /#/,$peps;
	map{
		#print join "\t", $pros_cnt{$pros},$pros,$_;
		#print "\n";
		push @cnt_pros_pep,[$pros_cnt{$pros},$pros,$_];
	}split /#/,$peps;
}

my %line_cnt; # {$pros_cnt{$pros} \t pro1;pro2;... \t pep}++ 
my %pep_cnt; # {pep}++;
		
map{
	my $line = join "\t",@{$_};
	$line_cnt{$line} = ++$pep_cnt{$_->[-1]};
}sort{$b->[0] <=> $a->[0] || $b->[1] cmp $a->[1]}@cnt_pros_pep;

my %label1_pros; # {pros => 1/na}
foreach my $line (keys %line_cnt){
	if($line_cnt{$line} == 1){
		$label1_pros{(split /\t/,$line)[1]} = 1;
	}
}

my %label2_pros;  # {pro1;pro2;... \t pep }
my %pep_cnts;  # {pep}++
foreach my $line (keys %line_cnt){
	my ($pros,$pep) = ((split /\t/,$line)[1,2]);
	if($label1_pros{$pros} == 1){
		$label2_pros{$pros."\t".$pep} = ++$pep_cnts{$pep};  
	}
}

print "Calculate Spc...\n";

my %spc_uniq_pep; # {pro1;pro2;... => spc of unique peptides}
my @data1;
foreach my $line (keys %label2_pros){
	my ($pros,$pep) = (split /\t/,$line);
	push @data1,[(split /\t/,$line),$pep_cnts{$pep},($pep_cnts{$pep} == 1 ? "Y" : "N")];
	if($pep_cnts{$pep} == 1){
		$spc_uniq_pep{$pros} += $pep_info{$pep};
	}
}

my %spc_share_pep; # {pro1;pro2;... => spc of shared peptides}
foreach my $data (@data1){
	my ($pros,$pep,$n,$u_s) = @{$data};
	if($u_s eq "N"){
		$spc_share_pep{$pep} += $spc_uniq_pep{$pros};
	}
}

my @data2;
foreach my $data(@data1){
	my ($pros,$pep,$n,$u_s) = @{$data};
	if($u_s eq "Y"){
		push @data2,[@{$data},$pep_info{$pep},$spc_uniq_pep{$pros},$spc_uniq_pep{$pros}];
	}else{
		push @data2,[@{$data},$pep_info{$pep},$spc_uniq_pep{$pros},$spc_share_pep{$pep}];		
	}
}

my %spc_pros; #
my %spc_pros_peps;  
foreach my $data(@data2){
	if($data->[-1] > 0 and length $data->[0]){ # remove the share protein(B) like: A(pep1,pep2),B(pep2,pep3),C(pep3,pep4);
		$spc_pros{$data->[0]} += sprintf ("%0.2f",$data->[4]*$data->[5]/$data->[-1]);
		if($data->[3] eq "Y"){
			$spc_pros_peps{$data->[0]}->[0] .= $data->[1]."; ";
			$spc_pros_peps{$data->[0]}->[1]++;
		}else{
			$spc_pros_peps{$data->[0]}->[2] .= $data->[1]."; ";
			$spc_pros_peps{$data->[0]}->[3]++;		
		}
	}
}

print "Calculate NSAF...\n";
my %SAF_pros;
my $total_SAF;
foreach my $pros (keys %spc_pros){
	if($spc_pros{$pros} > 0){
		$SAF_pros{$pros} = $spc_pros{$pros}/$database{max_len_pro($pros)}->[-1]; # SAF/PROTEIN LENGTH
		$total_SAF += $SAF_pros{$pros};
	}
}


open (group, ">$outfile") or die "Can't write the data to $outfile$!\n";
print group join"\t", qw/Proteins Genes Protein_with_max_len Protein_disc Gene Protein_len SpC SAF NSAF Ln(NSAF) uniq_evidence #uniq share_evidence #share/;
print group "\n";

foreach my $pros(keys %SAF_pros){
	my $max_pro = max_len_pro($pros);
	print group join "\t",$pros,genes($pros),$max_pro,@{$database{$max_pro}},$spc_pros{$pros},$SAF_pros{$pros},$SAF_pros{$pros}/$total_SAF,log($SAF_pros{$pros}/$total_SAF),@{$spc_pros_peps{$pros}};
	print group "\n";
}


sub max_len_pro{
	my $pros = shift;
	my $pro;
	my $max_len;
	map{
		if($database{$_}->[-1] > $max_len){
			$max_len = $database{$_}->[-1];
			$pro = $_;
		}
	}split /;/,$pros;
	return $pro;
}

sub genes{
	my $pros = shift;
	my $genes;
	map{$genes .= $database{$_}->[1].";"}split /;/,$pros;
	chop $genes;
	return $genes;
}
