#!/usr/bin/env perl

use strict;
use warnings;
use utf8;
use Getopt::Long;
use REST::Client;
use JSON;
use Pod::Usage;

my $VERSION = 0.01;

=head1 NAME

Targeted Mutation Identification : TaMi

=head1 SYNOPSIS

./TaMi.pl -g GENE_NAME [-k KMER_LENGTH] FASTQ_FILES

=head1 DESCRIPTION

TaMi can get a DNA Sequence from database 'ensembl' or from a fasta file, and giving a read length and a FASTQ file, will generate a VCF file with every SNPs that have been found. The sequence uses by TaMi is also stored as "Sequence.fa".

=head1 VERSION

0.01

=head1 AUTHORS

J.Audoux / A.Soriano

=head1 OPTIONS

  -man                            Print the manual.
  -help                           Print this help dialog.
  -v,--verbose                    Verbose mode.

input:

  -g,--gene STR                   Name of the Gene to look for in the database (hum...)
  -b,--bed  FILE                  Use a bed file as input instead of a gene name.
  -f,--FASTA FILE                 Use a FASTA file as input instead of a gene name 
  -s,--Species STR                The species you are working with (def : homo_sapiens)
  --grch37                        Use Ensembl API with the human genome in
                                  version GRCh37 (alias of Hg19)

output:

  -o,--output FILE                Name of the output vcf file. (default: STDOUT)

algorithmic features:

  -k,--kmer_length INT            K-mer length (default : 30)
  --disable-canonical             A k-mer and its reverse complement are
                                  considered as different sequences.the
                                  specified sequence. (default : false)

input filter:

  -F,--min-alternate-fraction N   Require at least this fraction of observations
                                  supporting an alternate allele in the in order
                                  to evaluate the position. (default: 0.01)
  -C,--min-alternate-count N      Require at least this count of observations
                                  supporting an alternate allele to evaluate the
                                  position. (default : 2)
  --min-coverage N                Require at least this coverage to process a
                                  site. (default: 0)
  --max-coverage N                Do not process sites with greater than this
                                  coverage. (default: no limit)

=cut

my ($help, $man, $verbose);
my $output_fileName;
my $geneName='';
my $inputFASTA='';
my $k=30;
my $specie='homo_sapiens';
my $disable_RC;
my $min_alternate_fraction = 0.01;
my $min_alternate_count = 2;
my $min_coverage = 0;
my $max_coverage = undef;
my $use_grch37;
my $ensembl_api_url = 'rest.ensembl.org';
my $nameBed = '';
my $commandLine = join " ", $0, @ARGV;
my $debug;

GetOptions( "v|verbose"           => \$verbose,
  "man"                         => \$man,
  "help"                        => \$help,
  "debug"                       => \$debug,

  # input
  "g|gene_name=s"               => \$geneName,
  "b|bed=s"                     => \$nameBed,
  # FIXME -f option is broken right now...
  #"f|FASTA-file=s"              => \$inputFASTA,
  "s|specie=s"                  => \$specie,
  "grch37"                      => \$use_grch37,

  # output
  "o|output-filename=s"         => \$output_fileName,

  # algorithmic features
  "k|kmer_length=i"             => \$k,
  "disable-canonical"           => \$disable_RC,

  # input filters
  "F|min-alternate-fraction=f"  => \$min_alternate_fraction,
  "C|min-alternate-count=i"     => \$min_alternate_count,
  "--min-coverage=i"            => \$min_coverage,
  "--max-coverage=i"            => \$max_coverage,

) or pod2usage (-verbose => 1);

#Now some test to check if everything's okay.
#No FASTQ
# FIXME No gene and no FASTA file, or both of them.

pod2usage(-verbose => 1)  if ($help);
pod2usage(-verbose => 2)  if ($man);

pod2usage(
  -message => "Only one input genome can be specified",
  -verbose => 1,
) unless ($geneName xor $nameBed);

my @FASTQ_files = @ARGV;

pod2usage(
  -message => "Mandatory argument 'FASTQ_file' is missing",
  -verbose => 1,
) if scalar @FASTQ_files == 0;

if ($use_grch37){# A dot is needed after the name of the genome.
  $ensembl_api_url = "grch37".$ensembl_api_url;
}

print STDERR "\n\n\t -- TaMi : Targeted Mutation Identification --\n\n\n";

my $client = REST::Client->new();

my %listingKmer=();

my @genomic_intervals;

# #############################################################################
# STEP 1 : GENOMIC INTERVALS
# 
# Retrieve genomic intervals that will be used to build the mutated k-mer hash
#
# Case 1: If the Name of a Gene has been specified
if ($geneName){
  my ($chromosome, $limInf, $limSup);
  $client->GET("http://".$ensembl_api_url."/xrefs/symbol/$specie/$geneName?content-type=application/json");

  my $xrefs = decode_json $client->responseContent();

  my $nbRefs=1;

  #This portion of code is not very nice... but it works !
  foreach my $ref (@{$xrefs}) {
    $client->GET("http://".$ensembl_api_url."/lookup/id/".$ref->{'id'}."?content-type=application/json");
    my $gene = decode_json $client->responseContent();
    #If only one key is present, we can store the information. If not, two different cases. The first case is when more than one entry are present, but that only one is usefull. The second case is when two different entry may be usefull for the user.
    if (scalar @{$xrefs} ==1){ #Easyest case, we store the values and then exit the loop
      $chromosome=$gene->{'seq_region_name'};
      $limInf = $gene->{'start'};
      $limSup = $gene->{'end'};
      print STDERR "\nName : ".$gene->{'display_name'}."\n";
      print STDERR "Description : ".$gene->{'description'}."\n";
      last;
    }
    else{ #This part will able us to see how many entrys are realy usefull if the doccument contains more than one. the condition can be translated by "if the entry is usefull, then nbRefs++. The informations are each time because they will be good if only one entry is present.
      if ((index($gene->{'source'}, 'havana') != -1) && ($gene->{'object_type'} eq 'Gene')){
        print STDERR "\n$nbRefs - Name : ".$gene->{'display_name'}."\n";
        print STDERR "Description : ".$gene->{'description'}."\n";
        $chromosome=$gene->{'seq_region_name'};
        $limInf = $gene->{'start'};
        $limSup = $gene->{'end'};
        $nbRefs++;
        next;
      }
      else{
        next;
      }
    }
  }

  if ($nbRefs>2){ #Only usefull if two or more usefull entrys are present. In this case we ask the user to choose the gene he want to use by typing his name. The user can see the name of each gene and the short description provided.
    print STDERR "\n".($nbRefs-1)." Genes have been found.\nType the name of the gene you wan't to work with : ";
    my $name = <STDIN>;
    my $count=0;
    chomp($name);
    foreach my $ref2 (@{$xrefs}){    
      $client->GET("http://".$ensembl_api_url."/lookup/id/".$ref2->{'id'}."?content-type=application/json");
      my $gene2 = decode_json $client->responseContent();
      if ((index($gene2->{'source'}, 'havana') != -1) && ($gene2->{'object_type'} eq 'Gene')){
        $count++;
        if ($count==$name){
          $chromosome=$gene2->{'seq_region_name'};
          $limInf = $gene2->{'start'};
          $limSup = $gene2->{'end'};
          $geneName = $gene2->{'display_name'};
          last;
        }
      }
    }
  }

  print STDERR "\nThe gene $geneName is located on chromosome $chromosome between position $limInf and $limSup.\n\n";

  push @genomic_intervals, { chr => $chromosome, start => $limInf, end => $limSup };

# case 2: If a bed file has been specified by the user
} elsif ($nameBed){
  print STDERR "Building the Kmer list using a Bed file...\n";
  open(my $inputBed, '<', $nameBed) or die("open $!");
  my @ligne;
  while (<$inputBed>){
    my ($chr,$start,$end) = split /\s+/, $_;
    # remove possible "chr" prefix
    $chr =~ s/^chr//i;
    # FIXME BED are half-open intervals, we should do $end-1
    push @genomic_intervals, { chr => $chr, start => $start + 1, end => $end };
  }
}

# Remove overlapping intervals
my @sorted_intervals = sort { $a->{chr} cmp $b->{chr} || $a->{start} <=> $b->{start} } @genomic_intervals;
my @merged_intervals;
my $curr_interval = shift @sorted_intervals;
for(my $i = 0; $i < scalar @sorted_intervals; $i++) {
  # If the two intervals are overlapping we merge them
  if($curr_interval->{chr} eq $sorted_intervals[$i]->{chr} && 
    $sorted_intervals[$i]->{start} <= $curr_interval->{end}) {
    # Update the "end" position
    $curr_interval->{end} = $sorted_intervals[$i]->{end} unless $sorted_intervals[$i]->{end} < $curr_interval->{end};
  # Otherwise we push the current interval to the merged intervals
  } else {
    push @merged_intervals, $curr_interval;
    $curr_interval = $sorted_intervals[$i];
  }
}
# Push the last interval
push @merged_intervals, $curr_interval;

# #############################################################################
# STEP 2 : MUTATED K-MER DICTIONNARY
#
# Create the mutated k-mer dictionnary from the merged intervals
my $kDiv2=int($k/2);
my $beenReverse=0;
my $refNucReverse;
my $kIsOdd=$k%2;
my $kmer;
foreach my $interval (@merged_intervals) {
  
  my $chromosome = $interval->{chr};
  my $limInf     = $interval->{start};
  my $limSup     = $interval->{end};

  print STDERR "Check interval : $chromosome:$limInf-$limSup\n" if $debug;

  $client->GET("http://".$ensembl_api_url."/sequence/region/$specie/$chromosome:$limInf..$limSup:1?content-type=text/plain");

  $inputFASTA = $client->responseContent();

  #Build a hash with every kmer with all the possible mutations on the midle base.
  #print STDERR ("Building the Kmer list...");

  my $seq_length = length $inputFASTA;

  for (my $i=0;$i<length($inputFASTA);$i++) #Build every kmer with a mutation on the middle base
  {
    $beenReverse  = 0;
    my $ref_kmer;
    my $mut_pos;
    if($i < $kDiv2) {
      $ref_kmer = substr($inputFASTA, 0, $k);
      $mut_pos  = $i;
    } elsif($i >= $seq_length - $kDiv2) {
      $ref_kmer = substr($inputFASTA, -$k);
      $mut_pos  = $k + $i - $seq_length;
    } else {
      $ref_kmer = substr($inputFASTA, $i - $kDiv2, $k);
      $mut_pos  = $kDiv2;
    }

    my $refNuc    = substr($ref_kmer, $mut_pos, 1);
    my $canonical_ref_kmer = canonicalKmer($ref_kmer);

    foreach my $nuc ("A", "G", "T", "C") {
      if($nuc ne $refNuc){ 

        # Create a mutation on the middle of the k-mer
        $kmer = mutationSimple($ref_kmer, $mut_pos, $nuc);

        # We check whether the mutated k-mer or its revcomp is smaller
        # even if this have very low chance to happen.
        $kmer = canonicalKmer($kmer);

        # FIXME we should check that the mutated k-mer is not already in the hash, otherwise
        # their is a collision that need to be handled
        if(exists $listingKmer{$kmer}) {
          print STDERR "Dupliacted entry, next...\n";
          next;
        }
        
        # FIXME we should store all these informations in a file to avoid having a huge amount
        # of memory used for nothing
        $listingKmer{$kmer}{'count'}    = 0;
        $listingKmer{$kmer}{'chromo'}   = $chromosome;
        $listingKmer{$kmer}{'ref_kmer'} = $canonical_ref_kmer;
        $listingKmer{$kmer}{'mut'}      = $nuc;
        $listingKmer{$kmer}{'position'} = $limInf + $i;

        
      #Will store the total number of kmer mapped derived from the ref.
      } else { 
        $listingKmer{$canonical_ref_kmer}{'count'}    = 0;
        $listingKmer{$canonical_ref_kmer}{'ref_nuc'}  = $refNuc;
      }
    }
  }
}

# Add all derived k-mer with one more mutation
foreach my $kmer (keys %listingKmer) {
  next if !defined $listingKmer{$kmer}{'ref_kmer'};
  for(my $j = 0; $j < $k; $j++) {
    my $refNuc = substr($kmer, $j, 1);
    foreach my $nuc ("A", "G", "T", "C") {
      if($nuc ne $refNuc) { 
        my $derived_kmer = canonicalKmer(mutationSimple($kmer, $j, $nuc));
        if(!exists $listingKmer{$derived_kmer}) {
          $listingKmer{$derived_kmer} = $listingKmer{$kmer};
        }
      }
    }
  }
}

# #############################################################################
# STEP 3 : QUANTIFY ABUNDANCE OF MUTATED K-MER
#
print STDERR ("\n\nReading FASTQ file...\n");

my $kmerRead;
my $nbRead=0;

foreach my $file (@FASTQ_files) {
  open(my $inputFASTQ, '<', $file) or die("open $!");
  while (<$inputFASTQ>){
    my $ligneQ = $_;

    # select only sequences lines
    if ($.%4 == 2){
      $nbRead++;
      chomp($ligneQ);
      for (my $i=0;$i<=length($ligneQ)-$k;$i++){
        $kmerRead = substr($ligneQ, $i, $k);
        $kmerRead = canonicalKmer($kmerRead) unless $disable_RC;

        # If this k-mer is defiened on the mutated k-mer dictionnary
        if(defined $listingKmer{$kmerRead}) {
          $listingKmer{$kmerRead}{'count'}++;

          # Update the reference k-mer, is this is a mutated k-mer
          if(defined $listingKmer{$kmerRead}{'ref_kmer'}) {
            $listingKmer{$listingKmer{$kmerRead}{'ref_kmer'}}{'count'}++;
          }
        }
      }
      if ($nbRead%50000==0){
        print STDERR "*";
      }
    }
  }
  close ($inputFASTQ);
}

print STDERR "\n$nbRead reads were present.\n\n";



# #############################################################################
# STEP 4 : FILTER AND PRINT OUTPUT
#
print STDERR ("Writing the output file as $output_fileName\n");

my $outputVCF;
if(defined $output_fileName) {
  open($outputVCF, '>', $output_fileName) or die ("open $!");
} else {
  $outputVCF = \*STDOUT;
}

# Print VCF headers
print $outputVCF "##fileformat=VCFv4.1\n";
print $outputVCF "##source=TaMI v$VERSION\n";
print $outputVCF '##commandline="'.$commandLine.'"',"\n";
print $outputVCF '##INFO=<ID=DP,Number=1,Type=Integer,',
                 'Description="Total read depth at the locus">',"\n";
print $outputVCF '##INFO=<ID=AC,Number=A,Type=Integer,',
                 'Description="Total number of alternate alleles in called genotypes">',"\n";
print $outputVCF '##INFO=<ID=AF,Number=A,Type=Float,',
                 'Description="Estimated allele frequency in the range (0,1]">',"\n";
print $outputVCF "#".join("\t",qw(CHROM POS ID REF ALT QUAL FILTER INFO)),"\n";

# Get kmers sorted by chr and positions
my @sorted_kmers = sort { $listingKmer{$a}{'chromo'} cmp $listingKmer{$b}{'chromo'} || 
                          $listingKmer{$a}->{'position'} <=> $listingKmer{$b}->{'position'} ||
                          $listingKmer{$a}->{'mut'} cmp $listingKmer{$b}{'mut'}
                  } grep { defined $listingKmer{$_}{'ref_kmer'} &&
                           $listingKmer{$_}{'count'} >= $min_alternate_count
                  } keys %listingKmer;

# Print VCF records
my $prev_mutation = "";
foreach my $key (@sorted_kmers){
  my $mutation = join("@",$listingKmer{$key}{'chromo'},
    $listingKmer{$key}{'position'},
    $listingKmer{$key}{'mut'},
  );
  next if $mutation eq $prev_mutation;
  $prev_mutation = $mutation;
  # Compute stats
  my $refNuc = $listingKmer{$listingKmer{$key}{'ref_kmer'}}{'ref_nuc'};
  my $DP = $listingKmer{$listingKmer{$key}{'ref_kmer'}}{'count'};
  my $AF = (($listingKmer{$key}{'count'})/($DP));

  # Apply filters
  next if $AF < $min_alternate_fraction;
  next if $DP < $min_coverage;
  next if defined $max_coverage && $DP > $max_coverage;

  # Print VCF line
  print $outputVCF join("\t",
    $listingKmer{$key}{'chromo'},
    $listingKmer{$key}{'position'},
    $key,
    $refNuc,
    $listingKmer{$key}{'mut'},
    join(";","DP=$DP","AF=$AF","AC=$listingKmer{$key}{'count'}")
  ),"\n";
}
print STDERR "\n\n --- END --- \n\n";


####################################################################
###########################***FONCTIONS***##########################
####################################################################

sub mutationSimple #Mute un nucléotide à une position donnée en un nucléotide de notre choix
{
  my ($sequence, $position, $nucleotide)=@_;
  substr($sequence, $position, 1, $nucleotide);
  return $sequence;
}

sub reverseComplement #Take a DNA sequence as input and output his reverse complement.
{
  my ($seq) = @_;
  chomp($seq);
  $seq =~ tr /atcgATCG/tagcTAGC/;
  $seq = reverse($seq);
  return $seq 
}

sub canonicalKmer {
  my $kmer = shift;
  my $beenReverse = 0;
  my  $reverseKmer = reverseComplement($kmer);
  if ($kmer gt $reverseKmer){
    $kmer         = $reverseKmer;
    $beenReverse  = 1;
  }
  return ($beenReverse, $kmer);
}
