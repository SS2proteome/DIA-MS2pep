package IsotopeDistribution;

use strict;
use Data::Dumper;

# mono_mass_list
my %atom;
$atom{'O'} = 15.99491463;
$atom{'H'} = 1.007825035;
$atom{'N'} = 14.003074;
$atom{'C'} = 12.00000;
$atom{'I'} = 126.904473;
$atom{'S'} = 31.9720707;
$atom{'P'} = 30.973762;
$atom{'hex'} = 162.052824;
$atom{'hexnac'} = 203.079373;

$atom{'C13'} = 13.0033554;
$atom{'N15'} = 15.000108;

my $proton = 1.00727646677;
my $H2O = $atom{'H'} * 2 + $atom{'O'};
my $NH3 = $atom{'H'} *3 + $atom{'N'};
my $CO = $atom{'C'} + $atom{'O'}; 
my $CAM = 2 * $atom{'C'} + 3 * $atom{'H'} + $atom{'O'} + $atom{'N'};

#<aminoacid_modification aminoacid="S" massdiff="365.1322" mass="452.1642" variable="Y"/>
#<aminoacid_modification aminoacid="T" massdiff="365.1322" mass="466.1799" variable="Y"/>

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
#$mono_aa{'M_'} = 131.040485;
# $mono_aa{'M_'} = 131.040485 + $atom{'O'}; # Oxidation
$mono_aa{'F'} = 147.068414;
$mono_aa{'P'} = 97.052764;
$mono_aa{'S'} = 87.032028;
$mono_aa{'T'} = 101.047679;
$mono_aa{'U'} = 150.95363;
$mono_aa{'W'} = 186.079313;
$mono_aa{'Y'} = 163.06332;
$mono_aa{'V'} = 99.068414;
# $mono_aa{'R_'} = $mono_aa{'R'} + 10.008269;  # heavy labeled R
# $mono_aa{'K_'} = $mono_aa{'K'} + 8.014199; # heavy labeled K
$mono_aa{'S_'} = $mono_aa{'S'} + 79.9663; # PhosSer
$mono_aa{'T_'} = $mono_aa{'T'} + 79.9663; # PhosThr
$mono_aa{'Y_'} = $mono_aa{'Y'} + 79.9663; # PhosTyr


my %AA_elements;

$AA_elements{'A'} = {'H',5,  'C',3,  'O', 1, 'N', 1};
$AA_elements{'C_'} = {'H', 8,  'C', 5,  'O', 2, 'N', 2, 'S', 1};  ## with cam iodoacetamide
$AA_elements{'C'} = {'H', 5,  'C', 3,  'O', 1, 'N', 1, 'S', 1};  # no modification!
$AA_elements{'D'} = {'H', 5,  'C', 4,  'O', 3, 'N', 1};
$AA_elements{'E'} = {'H', 7,  'C', 5,  'O', 3, 'N', 1};
$AA_elements{'F'} = {'H', 9,  'C', 9,  'O', 1, 'N', 1};
$AA_elements{'G'} = {'H', 3,  'C', 2,  'O', 1, 'N', 1};
$AA_elements{'H'} = {'H', 7,  'C', 6,  'O', 1, 'N', 3};
$AA_elements{'I'} = {'H', 11, 'C', 6,  'O', 1, 'N', 1};
$AA_elements{'K'} = {'H', 12, 'C', 6,  'O', 1, 'N', 2};
$AA_elements{'L'} = {'H', 11, 'C', 6,  'O', 1, 'N', 1};
$AA_elements{'M'} = {'H', 9,  'C', 5,  'O', 1, 'N', 1, 'S', 1};
$AA_elements{'M_'} = {'H', 9,  'C', 5,  'O', 2, 'N', 1, 'S', 1}; # oxidation
$AA_elements{'N'} = {'H', 6,  'C', 4,  'O', 2, 'N', 2};
$AA_elements{'P'} = {'H', 7,  'C', 5,  'O', 1, 'N', 1};
$AA_elements{'Q'} = {'H', 8,  'C', 5,  'O', 2, 'N', 2};
$AA_elements{'R'} = {'H', 12, 'C', 6,  'O', 1, 'N', 4};
$AA_elements{'S'} = {'H', 5,  'C', 3,  'O', 2, 'N', 1};
$AA_elements{'T'} = {'H', 7,  'C', 4,  'O', 2, 'N', 1};
$AA_elements{'U'} = {'H', 5,  'C', 3,  'O', 1, 'N', 1};  # and 'Se',1
$AA_elements{'V'} = {'H', 9,  'C', 5,  'O', 1, 'N', 1};
$AA_elements{'W'} = {'H', 10, 'C', 11, 'O', 1, 'N', 2};
$AA_elements{'Y'} = {'H', 9,  'C', 9,  'O', 2, 'N', 1};

#$AA_elements{'X'} = {'H', 11, 'C', 6,  'O', 1, 'N', 1}   # L or I
#$AA_elements{'Z'} = {'H',999, 'C',999, 'O',999,'N',999}  # E or Q, differ by 1 Da, always ignored
# heavy labeled c-terminal R, K
$AA_elements{'R_'} = {'H', 12, 'C13', 6,  'O', 1, 'N15', 4};
$AA_elements{'K_'} = {'H', 12, 'C13', 6,  'O', 1, 'N15', 2};
$AA_elements{'S_phos'} = {'H', 5,  'C', 3,  'O', 5, 'N', 1, 'P', 1};
$AA_elements{'T_phos'} = {'H', 7,  'C', 4,  'O', 5, 'N', 1, 'P', 1};
$AA_elements{'Y_phos'} = {'H', 9,  'C', 9,  'O', 5, 'N', 1, 'P', 1};

my %isotope;
# Molecular Isotopic Distribution Analysis (MIDAs) with Adjustable Mass Accuracy 
# Journal of The American Society for Mass Spectrometry
# January 2014, Volume 25, Issue 1, pp 57–70

$isotope{'H'} = {       
					1 => 0.999885,  #  1.0078250321     
					2 => 0.000115   #  2.0141017780
				};
$isotope{'C'} = {
					1 => 0.9893, #	 12.0000000000
					2 => 0.0107  #   13.0033548378
				};
$isotope{'O'} = {
					1 => 0.99757, #  15.9949146
					2 => 0.00038, #  16.9991312
					3 => 0.002050 #  17.9991603
				};
$isotope{'N'} = {
					1 => 0.99632, #  14.0030740052
					2 => 0.003680 #  15.0001088984
				};
$isotope{'S'} = {	
					1 => 0.949300, #  31.97207070
					2 => 0.007600, #  32.97145843
					3 => 0.04290,  #  33.96786665
					4 => 0.00020   #  35.96708062
				};
$isotope{'P'} = {
					1 => 1.0000
				};
# heavy 13C and 15N
$isotope{'C13'} = { 
					1 => 0,
					2 => 1.000  # 13.0033548378
				};
$isotope{'N15'} = { 
					1 => 0,
					2 => 1.000  # 15.0001088984
				};

#my $peptide = "LENITTGTYTIHAQR[166.109]";

#CAPTT($peptide);

sub CAPTT{
	my ($peptide,$mods) = @_;
	
	#print $peptide,"\n";
	
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
	#print $peptide,"##\n";
	
	#print Dumper $peptide;
	# $peptide =~ s/M\[147.0354\]/M_/g;
	# $peptide =~ s/C\[160.0307\]/C_/g;

	
	#my %uniq_;
	my ($P, $Q) = (1,1);
	my $PREC  = 0.000001;
	my $distribution;
	my $Captt;
	$Captt->{1} = 1;
	my $monoMass;
	
	
	my %AAs;
	my @AAs = split /(?=[[:upper:]]_?)/,$peptide;
	#my @seq_arr = split /(?=[[:upper:]]_?)/,$sequence;

	#print Dumper "\@AAs",\@AAs;
	my %natom;
	foreach my $aa(@AAs){
		$AAs{$aa}++;
		$monoMass += $mono_aa{$aa};
		#print Dumper $aa,$AA_elements{$aa};
		foreach my $atom (keys %{$AA_elements{$aa}}){
			$natom{$atom} += $AA_elements{$aa}{$atom};
		}
	}
	$monoMass += $atom{'O'} + $atom{'H'} + $proton;
	
	#print "\$monoMass = $monoMass\n";
	#print Dumper "\%atom",\%natom;
	#print Dumper \%mono_aa;
	#%natom = ('C', 254, 'H', 377, 'N', 65, 'O', 75, 'S', 6);
	#%natom = ('C', 1185, 'H', 1850, 'N', 282, 'O', 339, 'S', 18);
	#%natom = ('C', 28, 'H', 49, 'N', 7, 'O', 8, 'S', 1);
	#print Dumper "\%atom",\%natom;
	
	
	foreach my $J (keys %natom){ # [H,O,C,N,S]
		#print "\$J = $J\n";
		foreach my $I (1..$natom{$J}){ # number of atoms
			#print Dumper "\$distribution",$distribution;
			undef $distribution;
			my $Npeak;
			$Npeak = scalar (keys %{$isotope{$J}}) ;
			#print "\$Npeak, $Npeak\n";
			foreach my $K ($P .. $Q){ # Isotopic distribution
				foreach my $L (1 .. $Npeak){
					$distribution->{$K + $L - 1} = $distribution->{$K + $L - 1} + $Captt->{$K} * $isotope{$J}{$L};
					
					#print "\$distribution->\{$K + $L - 1\} = \$distribution->\{$K + $L - 1\} + \$Captt->\{$K\} * \$isotope\{$J\}\{$L\};\n";
					#print Dumper "\$distribution",$distribution,$isotope{$J}{$K};
					
				}
			}
			$Q = $Q + $Npeak - 1;
			#print "\$J: $J; \$I: $I; \$Q = $Q\n";
			
			my $max = 0 ;
			foreach my $K ($P .. $Q){
				$max = $distribution->{$K} if $max < $distribution->{$K};
			}
			
			#print "\$max, $max\n";
			#print Dumper "\$--distribution",$distribution;
			foreach my $K ($P .. $Q){
				$distribution->{$K} /= $max;
			}
			
			#print Dumper "\$%%distribution",$distribution;
			
			#my %rm_list;
			#foreach my $K ($P .. $Q){
			my $K = $P;
			while($K <= $Q){
				#if D(K) > PREC then P = K : K = Q
				#print "$distribution->{$K} > $PREC\n";
				#print "\$K: $K; \$P: $P; \$Q: $Q\n";
				if ($distribution->{$K} > $PREC){
					$P = $K;
					$K = $Q;
				}
				$K++;
				#print "---\$K: $K; \$P: $P; \$Q: $Q\n";
			}
			
			my $K = $Q;
			while($K >= $P){
			#foreach my $K ($Q .. $P){
				if ($distribution->{$K} > $PREC){
					$Q = $K;
					$K = $P;
				}
				$K--;
			}
			
			#print Dumper "\$oldCaptt",$Captt;
			undef $Captt;
			undef $K;
			#print "---\$Q: $Q; \$P: $P;\n";
			foreach my $K ($P .. $Q){
				$Captt->{$K} = $distribution->{$K};
			}
			#print Dumper "\$newCaptt",$Captt;
		}
	}
	
	#print Dumper "\$newCaptt",$Captt;
	#print Dumper "\$distribution",$distribution;
	my @iso_dis;
	foreach my $K (sort{$a <=> $b}(keys %{$Captt})){
		push @iso_dis,$Captt->{$K};
		#print join "\t", $K, $Captt->{$K};
		#print "\n";
	}
	return [@iso_dis];
}


=head 
Analytica Chimica Acta
Volume 247, Issue 1, 14 June 1991, Pages 107-119
Calculation of isotope distributions in mass spectrometry. A trivial solution for a non-trivial problem

REM --------------------------------------------------------
REM 			calculation of isotope distributions
REM --------------------------------------------------------

REM list of atoms = H, He, Li, Be, B, C, N, O, ...

P = 1: Q = 1					:REM 	initialize P (=low), Q (=high)
PREC = 0.000001					:REM 	define pruning treshold factor PREC
CAPTT(1) = 1					:REM 	CAPTT = calculated pattern
for J = 1 to natom				:REM 	natom = number of atoms in list
	for I = 1 to C(J)			:REM 	C(J) = number of atoms J in the molecules
		ERASE D					:REM 	D = dummy variable
		for K = P TO Q			:REM 	calculation of new pattern
			for L = 1 to NPEAK(J) :REM 	NPEAK(J) = number of heaviest isotope
				D(K+L-1) = D(K+L-1) + CAPTT(K) * ABUND(J,L)
			next L-1			:REM 	ABUND(J,L) = abundances of iostopes of the atom type J	
		next K									
		Q = Q + NPEAK(J) - 1	:REM 	Q = new number of heaviest iostope of the pattern
		max = 0 				:REM 	set variable
		for K = P TO Q			:REM 	find largest value of D(K)
			if(D(K) > max) then max = D(K)
		next
		for K = P to Q			:REM 	normalize D(K) to D(max) = 1
			D(K) = D(K)/max 
		next
		for K = P to Q 			:REM 	eliminate small peaks at left side
			if D(K) > PREC then P = K : K = Q
		next K
		for K = Q TO P step-1		:REM 	eliminate small peaks at right side
			if D(K) > PREC then Q = K : K = P
		next K
		ERASE CAPTT				:REM 	erase old isotope peak pattern
		for K = P TO Q 
			CAPTT(K) = D(K)		:REM 	create new isotope peak pattern
		next
	next						:REM 	next atom of same atom type
next							:REM 	next atom type
end
=cut



1;
