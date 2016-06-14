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

./TaMi.pl Gene FastQ KmerLength (for the moment)

=head1 DESCRIPTION

ToDo

=head1 VERSION

0.0

=head1 AUTHORS

J.Audoux / A.Soriano

=head1 OPTIONS

  -o,--output-dir         Output directory (DEFAULT: 'simCT_simulation') 
  -k                      Kmer Length
  -g,--Gene               Name of the Gene to look for in the database (hum...)
  -F,--FASTA              Use a FASTA file as input instead of a gene name 
  -Q,--FASTQ              FASTQ File to work with
  -s,--Species            The species you are working with

=cut


my $geneName=$ARGV[0];

#GetOptions ('inputFASTA' => \$inputFASTA, 'inputFASTQ' => \$all, 'k' => \$k);
#Il faudra rendre l'utilisation d'un fichier facultative.
#my $organism = $ARGV[1]; Different organism from human can also be used.

my $client = REST::Client->new();

$client->GET("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/$geneName?content-type=application/json");

#print STDERR Dumper($client->responseContent());

my $xrefs = decode_json $client->responseContent();

my $chromosome;
my $limInf;
my $limSup;
my $nbRefs=0;

#This portion of code is not very nice... but it works !
foreach my $ref (@{$xrefs}) {
    $client->GET("https://rest.ensembl.org/lookup/id/".$ref->{'id'}."?content-type=application/json");
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
    print STDERR "\n$nbRefs Genes have been found.\nType the name of the gene you wan't to work with : \n";
    my $name = <STDIN>;
    chomp($name);
    foreach my $ref2 (@{$xrefs})
    {    
        $client->GET("https://rest.ensembl.org/lookup/id/".$ref2->{'id'}."?content-type=application/json");
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

$client->GET("https://rest.ensembl.org/sequence/region/human/$chromosome:$limInf..$limSup:1?content-type=text/plain");

my $inputFASTA = $client->responseContent();

open(my $inputFASTQ, '<', $ARGV[1]) or die("open $!");
open(my $outputVCF, '>', 'Output.vcf') or die ("open $!");
my $k=$ARGV[2];

my $kmer;
my $nbKmer;
my $name; #Permet de stocker le nom du read
my $position=0; #Permet de stocker la position du match...

#Construction d'une Hash devant stocké tous les Kmer non mutés. Il faudra gérer les bords aussi.

my %listingKmer=();

print STDERR ("Building the Kmer list...\n");

for (my $i=0;$i<length($inputFASTA)-$k+1;$i++) #Construction de tous les kmer mutés au centre.
{
	my $ref_kmer = substr ($inputFASTA, $i, $k);
	my $refNuc= substr($ref_kmer, (length($ref_kmer)/2), 1);
	foreach my $nuc ("A", "G", "T", "C")
	{
		if($nuc ne $refNuc)
	   	{
			$kmer = mutationSimple($ref_kmer,(length($ref_kmer)/2) , $nuc);
			$listingKmer{$kmer}{'count'}=0; #Construire ici les Kmer mutés au centre. Technique de rat d'égout où on code tout en dur !
			$listingKmer{$kmer}{'ref_kmer'}=$ref_kmer;
			$listingKmer{$kmer}{'mut'}=$nuc;
			$listingKmer{$kmer}{'position'}=int($i+($k/2)+1+$limInf); #Idiot ?
		}
		else
		{
			$listingKmer{$ref_kmer}{'count'}=0;
		}
	}
}

my $nbRead=0;

print STDERR ("Reading FASTQ file...\n");

my $kmerRead;

while (<$inputFASTQ>) #On lit le FastQ
{
	my $ligneQ = $_;
	if ($.%4 == 2) #Selection des reads
	{
		$nbRead++;
		chomp($ligneQ);
		for (my $i=0;$i<length($ligneQ);$i++)
		{
			$kmerRead = substr($ligneQ, $i, $k);
			if (defined($listingKmer{$kmerRead}))
			{
				$listingKmer{$kmerRead}{'count'}++;
			}
		}
	}
}

print STDERR "$nbRead reads were present.\n\n";

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
			$refNuc = substr($listingKmer{$key}{'ref_kmer'}, $k/2, 1);
			$DP = $listingKmer{$listingKmer{$key}{'ref_kmer'}}{'count'}+$listingKmer{$key}{'count'};
			$AF = $listingKmer{$listingKmer{$key}{'ref_kmer'}}{'count'}/$DP;
			print $outputVCF "$chromosome\t$listingKmer{$key}{'position'}\t$key\t$refNuc\t$listingKmer{$key}{'mut'}\tDP=$DP;AF=$AF\n";
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
	my @nucleotides = ('a', 'g', 'c', 't');
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