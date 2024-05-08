#!/usr/bin/perl

# here is a Perl script that takes in sequences of yeast, converts to AA and reports AA statistics info

mkdir 'ProcessedAASequences';
$dir = "YeastGenes"; opendir(DIR,$dir) or die "can't open directory $dir:$!"; print"\n";
open(OUTFILE2, ">myAAsequence.txt") or die "can't open file\n";

while ($filename = readdir DIR){ 
  $ORFname = substr($filename, 0, 7);
  print "\nmy ORF name is "."$ORFname\n";
  $filelocation = "./YeastGenes/"."$filename";
  if (length $ORFname == 7){
    open(INFILE, $filelocation) or die "Cannot open file";
  }else {next;}
    open (OUTFILE, ">"."./ProcessedAASequences/"."$filename") || die " could not open output file\n";

#loop for each file
  while(<INFILE>){
    chomp;
    #Calculate GC Content
    my $totalCount = 0;
    my $CGCount = 0;
    my $DNA = <INFILE>;
    print($DNA . "\n");
    print OUTFILE ($DNA . "\n");
    my $position = 0;
    my $DNAsize = length($DNA);
    my $counter = 0;
    while($counter !=  $DNAsize){
      my $base = substr($DNA, $position, 1);
      if($base eq "A" || $base eq "T"  || $base eq "C"  ||  $base eq "G" ){
        $totalCount++;
      }
      if($base eq "G" || $base eq "C"){
        $CGCount++;
      }
      $position++;
      $counter++;
    }
    #print GC Content Info;
print "countGC "."$ORFname\t"."$CGCount\n";
print "countTTL "."$ORFname\t"."$totalCount\n";
$freqGC = $CGCount/$totalCount;
print OUTFILE "countGC "."$ORFname\t"."$CGCount\n";
print OUTFILE "countTTL "."$ORFname\t"."$totalCount\n";
print OUTFILE "freqGC "."$ORFname\t"."$freqGC\n";

#Convert each Yeast String into AA
# Define the codon table
my %codon_table = (
    'TTT' => 'F', 'TTC' => 'F', 'TTA' => 'L', 'TTG' => 'L',
    'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
    'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'ATG' => 'M',
    'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
    'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
    'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
    'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
    'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
    'TAT' => 'Y', 'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*',
    'CAT' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
    'AAT' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
    'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
    'TGT' => 'C', 'TGC' => 'C', 'TGA' => '*', 'TGG' => 'W',
    'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
    'AGT' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
    'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',
);
# Split the DNA sequence into codons
my @codons = $DNA =~ /.../g;
# Translate each codon to amino acid
my $protein = '';
foreach my $codon (@codons) {
    if (exists $codon_table{$codon}) {
        $protein .= $codon_table{$codon};
    } else {
        die "Error: Invalid codon encountered: $codon\n";
    }
}
# Print the resulting protein sequence
print OUTFILE "Protein sequence: $protein\n";

# Calculate the percentage of polar and non-polar amino acids
my $polar_count = ($protein =~ tr/DEHKNQRSTY/DEHKNQRSTY/);  # Aspartic Acid, Glutamic Acid, Histidine, Lysine, Asparagine, Glutamine, Serine, Threonine
my $non_polar_count = length($protein) - $polar_count;
my $total_amino_acids = length($protein);
my $polar_percent = ($polar_count / $total_amino_acids) * 100;
my $non_polar_percent = ($non_polar_count / $total_amino_acids) * 100;
#Print polar non polar calculations
print OUTFILE "Percentage of Polar Amino Acids: $polar_percent%\n";
print OUTFILE "Percentage of Non-Polar Amino Acids: $non_polar_percent%\n";

# Calculate the percentage of positively and negatively charged amino acids
my $positively_charged_count = ($protein =~ tr/KRH/KRH/);
my $negatively_charged_count = ($protein =~ tr/DE/DE/);
my $total_amino_acid = length($protein);
my $positively_charged_percent = ($positively_charged_count / $total_amino_acid) * 100;
my $negatively_charged_percent = ($negatively_charged_count / $total_amino_acid) * 100;
#Print +/- info
print OUTFILE "Number of Positively Charged Amino Acids: $positively_charged_count\n";
print OUTFILE "Percentage of Positively Charged Amino Acids: $positively_charged_percent%\n";
print OUTFILE "Number of Negatively Charged Amino Acids: $negatively_charged_count\n";
print OUTFILE "Percentage of Negatively Charged Amino Acids: $negatively_charged_percent%\n";

#sleep(1);
print "\n\n";
close OUTFILE;
close INFILE;
  }
} # end while loop
print "end program\n";
close OUTFILE2;
exit;