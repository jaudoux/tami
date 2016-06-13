#!/usr/bin/env perl
#===============================================================================
#
#         FILE: FastaMere.pl
#
#        USAGE: ./FastaMere.pl  
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

#my $inputFASTA;
#my $inputFASTQ;
#my $k;

#GetOptions ('inputFASTA' => \$inputFASTA, 'inputFASTQ' => \$all, 'k' => \$k);
#Il faudra rendre l'utilisation d'un fichier facultative.

open(my $inputFASTA, '<', $ARGV[0]) or die("open $!");
open(my $inputFASTQ, '<', $ARGV[1]) or die("open $!");
open(my $outputVCF, '>', 'Output.vcf') or die ("open $!");
my $k=$ARGV[2];

my $Fasta;
my $kmer;
my $nbKmer;
my $name; #Permet de stocker le nom du read
my $position=0; #Permet de stocker la position du match...
my $chromosomeName;

print STDERR ("\n\nStockage du Fasta...\n");

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
		$chromosomeName = $_;
	}
}
close($inputFASTA);

#Construction d'une Hash devant stocké tous les Kmer non mutés. Il faudra gérer les bords aussi.

my %listingKmer=();

print STDERR ("Construction de la liste des Kmer mutes\n");

for (my $i=0;$i<length($Fasta)-$k+1;$i++) #Construction de tous les kmer mutés au centre.
{
	my $ref_kmer = substr ($Fasta, $i, $k);
	my $refNuc= substr($ref_kmer, (length($ref_kmer)/2), 1);
	foreach my $nuc ("A", "G", "T", "C")
	{
		if($nuc ne $refNuc)
	   	{
			$kmer = mutationSimple($ref_kmer,(length($ref_kmer)/2) , $nuc);
			$listingKmer{$kmer}{'count'}=0; #Construire ici les Kmer mutés au centre. Technique de rat d'égout où on code tout en dur !
			$listingKmer{$kmer}{'ref_kmer'}=$ref_kmer;
			$listingKmer{$kmer}{'mut'}=$nuc;
			$listingKmer{$kmer}{'position'}=int($i+($k/2)+1); #Idiot ?
		}
		else
		{
			$listingKmer{$ref_kmer}{'count'}=0;
		}
	}
}

my $nbRead=0;

print STDERR ("Lecture du FastQ...\n");

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

print STDERR "Il y avait $nbRead reads.\n";

close ($inputFASTQ);

my $compteur=0; 

print STDERR ("Ecriture...\n");

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
			print $outputVCF "$chromosomeName\t$listingKmer{$key}{'position'}\t$key\t$refNuc\t$listingKmer{$key}{'mut'}\tDP=$DP;AF=$AF\n";
		}
}

print STDERR "\n\n --- Fini --- \n\n";

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
