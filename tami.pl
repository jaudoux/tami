#!/usr/bin/env perl
#===============================================================================
#
#         FILE: compteKmer.pl
#
#        USAGE: ./compteKmer.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 06/06/2016 17:37:38
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;
use REST::Client;
use Data::Dumper;
use JSON;
use Pod::Usage;

=pod

=head1 NAME

Targeted Mutation Identification : TaMi

=head1 SYNOPSIS

./TaMi.pl -g geneName -k Kmer_Length -q FASTQ_File

=head1 DESCRIPTION

TaMi can get a DNA Sequence from database 'ensembl' or from a fasta file, and giving a read length and a FASTQ file, will generate a VCF file with every SNPs that have been found.

=head1 VERSION

0.01

=head1 AUTHORS

J.Audoux / A.Soriano

=head1 OPTIONS

  -man                    Print the manual
  -help                   Print the... help !
  -v,--verbose            Verbose...
  -o,--output-dir         Output directory 
  -k,--kmer_length        Kmer Length (def : 29)
  -g,--Gene               Name of the Gene to look for in the database (hum...)
  -f,--FASTA              Use a FASTA file as input instead of a gene name 
  -q,--FASTQ              FASTQ File to work with
  -s,--Species            The species you are working with (def : human)
  -r,--Reverse_complement Specified if tami must use the reverse complement of the specified sequence. Set to 0 to disable. (def : true)
  -c,--cut                Set the AF above which a mutation will be selected. All mutations with a AF value lower than this one will not be selected : value must be between 0 and 1, where 0 means disabled. (def 0.01)
  -n,--nbCut              Specified the number of map above which a kmer will be chosen (def : 2)
  -a,--another_genome     Worst name ever. Specified if another version of the human genome must be used.

=cut

my ($help, $man, $verbose);
my $output_dir;
my $geneName='';
my $inputFASTA='';
my $refFASTQ='';
my $k=29;
my $specie='human';
my $RC = 1;
my $cut = 0.01;
my $nbCut = 2;
my $genome = 'rest.ensembl.org';

GetOptions( "v|verbose"           => \$verbose,
            "man"                 => \$man,
            "help"                => \$help,
            "o|output-dir=s"      => \$output_dir,
            "k|kmer_length=i"     => \$k,
            "g|gene_name=s"       => \$geneName,
            "f|FASTA_file=s"      => \$inputFASTA,
            "q|FASTQ_file=s"      => \$refFASTQ,
            "s|specie=s"          => \$specie,
            "r|reverse_c=i"       => \$RC,
            "c|cut=f"             => \$cut,
            "n|nbCut=i"           => \$nbCut,
            "a|another_genome=s"  => \$genome,
        ) or pod2usage (-verbose => 1);

#Now some test to check if everything's okay.
        #No FASTQ
        #No gene and no FASTA file, or both of them.
 
pod2usage(-verbose => 1)  if ($help);
pod2usage(-verbose => 2)  if ($man);

pod2usage(
    -message => "Mandatory argument 'FASTQ_file' is missing",
    -verbose => 1,
) unless defined $refFASTQ;

pod2usage(
    -message => "Only one input genome can be specified",
    -verbose => 1,
) unless ($geneName xor $inputFASTA);

open(my $inputFASTQ, '<', $refFASTQ) or die("open $!");
open(my $outputVCF, '>', 'Output.vcf') or die ("open $!");

if ($genome ne 'rest.ensembl.org' )# A dot is needed after the name of the genome.
{
    $genome =  $genome.".rest.ensembl.org";
}

my $client = REST::Client->new();



$client->GET("http://".$genome."/xrefs/symbol/homo_sapiens/$geneName?content-type=application/json");

#print STDERR Dumper($client->responseContent());

my $xrefs = decode_json $client->responseContent();

my $chromosome;
my $limInf;
my $limSup;
my $nbRefs=0;

#This portion of code is not very nice... but it works !
foreach my $ref (@{$xrefs}) {
    $client->GET("http://".$genome."/lookup/id/".$ref->{'id'}."?content-type=application/json");
	my $gene = decode_json $client->responseContent();
    #If only one key is present, we can store the information. If not, two different cases. The first case is when more than one entry are present, but that only one is usefull. The second case is when two different entry may be usefull for the user.
    if (scalar keys $xrefs ==1){ #Easyest case, we store the values and then exit the loop
        $chromosome=$gene->{'seq_region_name'};
        $limInf = $gene->{'start'};
   	    $limSup = $gene->{'end'};
        print STDERR "\nName : ".$gene->{'display_name'}."\n";
        print STDERR "Description : ".$gene->{'description'}."\n";
        last;
    }
    else{ #This part will able us to see how many entrys are realy usefull if the doccument contains more than one. the condition can be translated by "if the entry is usefull, then nbRefs++. The informations are each time because they will be good if only one entry is present.
	    if ((index($gene->{'source'}, 'havana') != -1) && ($gene->{'object_type'} eq 'Gene') ){
            print STDERR "\nName : ".$gene->{'display_name'}."\n";
            print STDERR "Description : ".$gene->{'description'}."\n";
            $chromosome=$gene->{'seq_region_name'};
            $limInf = $gene->{'start'};
   	        $limSup = $gene->{'end'};
            $nbRefs++;
            next;
        }
        else
        {
            next;
        }
    }
}

if ($nbRefs>1) #Only usefull if two or more usefull entrys are present. In this case we ask the user to choose the gene he want to use by typing his name. The user can see the name of each gene and the short description provided.
{
    print STDERR "\n$nbRefs Genes have been found.\nType the name of the gene you wan't to work with : ";
    my $name = <STDIN>;
    chomp($name);
    foreach my $ref2 (@{$xrefs})
    {    
        $client->GET("http://".$genome."/lookup/id/".$ref2->{'id'}."?content-type=application/json");
        my $gene2 = decode_json $client->responseContent();
	    if ((index($gene2->{'source'}, 'havana') != -1) && ($gene2->{'object_type'} eq 'Gene') )
        {
            if ($gene2->{'display_name'} eq $name)
            {
    	        $chromosome=$gene2->{'seq_region_name'};
    	        $limInf = $gene2->{'start'};
    	        $limSup = $gene2->{'end'};
             }
         }
    }
}

print STDERR "\nThe gene $geneName is located on chromosome $chromosome between position $limInf and $limSup.\n\n";

$client->GET("http://".$genome."/sequence/region/human/$chromosome:$limInf..$limSup:1?content-type=text/plain");

#Now the genome will be analysed.

$inputFASTA = $client->responseContent();
open(my $FastaGenome, '>', 'Sequence.fa') or die ("open $!"); #Récupère la séquence, pour moi.
print $FastaGenome $inputFASTA;
close $FastaGenome;
my $kmer;
my $nbKmer;
my $name; #Permet de stocker le nom du read
my $position=0; #Permet de stocker la position du match...
my $kDiv2=$k/2;
my $beenReverse=0;
my $refNucReverse;

#Construction d'une Hash devant stocker tous les Kmer non mutés. Il faudra gérer les bords aussi.

my %listingKmer=();

print STDERR ("Building the Kmer list...\n");

for (my $i=0;$i<length($inputFASTA)-$k+1;$i++) #Build every kmer with a mutation on the middle base
{
    my $refNuc;
	my $ref_kmer = substr ($inputFASTA, $i, $k);
    if ($k%2!=0){ #If k is an even number, then we must do this... If we don't the wrong base will be changed.
    	$refNuc= substr($ref_kmer, $kDiv2, 1);
    }
    else{
        $refNuc=substr($ref_kmer, $kDiv2-1, 1);
    }
    $beenReverse=0;

    if ($RC != 0) #Make the RC of the RefKMer, and use it if he is lower than the ref as a string.
    {
        my  $reverseKmer = reverseComplement($ref_kmer);
        if ($ref_kmer gt $reverseKmer)
        {
            $ref_kmer = $reverseKmer;
            $beenReverse=1;
        }
    }
    if (!$beenReverse && $k%2!=0){ #Theses 3 conditions will choose which base will be the reference one depending on k and if ref_kmer has been reversed or not.
        $refNucReverse=substr($ref_kmer, $kDiv2, 1);
    }
    elsif(!$beenReverse && $k%2 == 0){
        $refNucReverse=substr($ref_kmer, $kDiv2-1, 1);
    }
    else{
        $refNucReverse=substr($ref_kmer, $kDiv2, 1); #Same as the if... 
    }
    if (!exists($listingKmer{$ref_kmer})){
    	foreach my $nuc ("A", "G", "T", "C")#I wonder if making two distinct loop with and without the reverse won't be faster... in the code above the $RC is tested each time...
	    {
		    if($nuc ne $refNucReverse) #We must be carefull here. If the wrong nucleotide is used as a reference, then everything will be wrong... 
	   	    {
                if (!$beenReverse && $k%2!=0){
			        $kmer = mutationSimple($ref_kmer,$kDiv2, $nuc);
                }
                elsif($k%2==0 && !$beenReverse){
                    $kmer = mutationSimple($ref_kmer, $kDiv2-1, $nuc);
                }
                else{
                    $kmer = mutationSimple($ref_kmer, $kDiv2, $nuc);
                }
			    $listingKmer{$kmer}{'count'}=0;
			    $listingKmer{$kmer}{'ref_kmer'}=$ref_kmer;
                if ($beenReverse){
		    	    $listingKmer{$kmer}{'mut'}=reverseComplement($nuc);
                }
                else{
                    $listingKmer{$kmer}{'mut'}=$nuc;
                }
			    $listingKmer{$kmer}{'position'}=int($i+$kDiv2+$limInf); #Idiot ?
		    }
		    else #Will store the total number of kmer mapped derived from the ref.
		    {
			$listingKmer{$ref_kmer}{'count'}=0;
            $listingKmer{$ref_kmer}{'ref_nuc'}=$refNuc;
		    }
        }
    }
}

my $nbRead=0;

print STDERR ("Reading FASTQ file...\n");

my $kmerRead;

while (<$inputFASTQ>) #Reading the fastQ file.
{
	my $ligneQ = $_;
	if ($.%4 == 2) #Read selection.
	{
		$nbRead++;
		chomp($ligneQ);
		for (my $i=0;$i<=length($ligneQ)-$k;$i++) #Build the Kmer of the selected read.
		{
			$kmerRead = substr($ligneQ, $i, $k);
            if ($RC!=0)#Reverse Complement
            {
                my $kmerReverseRead = reverseComplement($kmerRead);
                if ($kmerRead gt $kmerReverseRead)
                {
                    $kmerRead = $kmerReverseRead;
                }
            }
			if (defined($listingKmer{$kmerRead})) #If this part of the read can be found somwhere in the hash.
			{
				$listingKmer{$kmerRead}{'count'}++;
                # print Dumper($listingKmer{$kmerRead}{'count'});
                if (defined($listingKmer{$kmerRead}{'ref_kmer'})) #If the Kmer has a reference kmer or a position, it means that it's a mutated one. So we will increment the ref in order to have an access to the DP.
                {
#                    print "$listingKmer{$listingKmer{$kmerRead}{'ref_kmer'}}{'count'}\n";
                    $listingKmer{$listingKmer{$kmerRead}{'ref_kmer'}}{'count'}++;
                }
			}
		}
        if ($nbRead%50000==0){
            print STDERR "*";
	    }
    }
}

print STDERR "\n$nbRead reads were present.\n\n";

close ($inputFASTQ);

my $compteur=0; 

print STDERR ("Writing the output file as Output.vcf...\n");

print $outputVCF "Chrom\tPos\tID\tRef\tAlt\tInfo\n";

my $refNuc;
my $DP; #Tout
my $AF; #Ref/somme (une moyenne quoi...)

foreach my $key ( sort {$listingKmer{$a}->{'position'} <=> $listingKmer{$b}->{'position'}} grep { defined $listingKmer{$_}{'position'}} keys %listingKmer)
{
		if ($listingKmer{$key}{'count'}>0 && defined($listingKmer{$key}{'ref_kmer'})) #On peut rajouter cette condition dans le grep. 
		{
            $refNuc = $listingKmer{$listingKmer{$key}{'ref_kmer'}}{'ref_nuc'};
            $DP = $listingKmer{$listingKmer{$key}{'ref_kmer'}}{'count'};#$listingKmer{$key}{'ref_kmer'}{'count'}
  			$AF = (($listingKmer{$key}{'count'})/($DP));
            if ($cut <= $AF && $listingKmer{$key}{'count'}>=$nbCut && defined($listingKmer{$key}{'ref_kmer'}))
            {
			    print $outputVCF "$chromosome\t$listingKmer{$key}{'position'}\t$key\t$refNuc\t$listingKmer{$key}{'mut'}\tDP=$DP;AF=$AF\n";
            }
		}
}

print STDERR "\n\n --- END --- \n\n";


####################################################################
###########################***FONCTIONS***##########################
####################################################################

sub mutation
{
	my ($sequence, $position)=@_;
	my $seqLength = length($sequence);
	my @nucleotides = ('A', 'G', 'C', 'T');
	my $nucleotideRand=$nucleotides[rand(@nucleotides)];
	substr($sequence, $position, 1, $nucleotideRand);
	return $sequence;
}

sub mutationCentrale #Prend une chaine de caractère, et en modifie le nucléotide central, se rapelle elle même si le résultat est le même.
{
	my ($maSequence)=@_;
	my $taille = (length($maSequence)/2);
	return mutation($maSequence, $taille);
}

sub mutationSimple #Mute un nucléotide à une position donnée en un nucléotide de notre choix
{
	my ($sequence, $position, $nucleotide)=@_;
	substr($sequence, $position, 1, $nucleotide);
	return $sequence;
}

sub importFastaFile #UNTESTED ! Take the name of a FASTA file as input, and output the sequenc without any header or useless char.
{
    open(my $inputFASTA, '<', @_) or die("open $!");
    my $Fasta;

    while (<$inputFASTA>) #Lecture du fichier Fasta et Stockage de son contenu dans $Fasta 
    {
	    if ($_ !~ /^>/) #Dégage les lignes qui correspondent à une en-tête... Mais rajoute un test à chaque fois donc pas glorieux, dans notre cas faudrat juste sauter la première lignecar on a qu'un gène avec les données de Benois. =~ /^>/
	    {
	        chomp($_);
	        $Fasta=$Fasta.$_;	
	    }
	    else
	    {
		    chomp($_);
	    }
    }
    return $Fasta;
}

sub reverseComplement #Take a DNA sequence as input and output his reverse complement.
{
    my ($seq) = @_;
    chomp($seq);
    $seq =~ tr /atcgATCG/tagcTAGC/;
    $seq = reverse($seq);
    return $seq 
}
