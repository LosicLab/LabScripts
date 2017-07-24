#!/usr/bin/perl
use warnings;
use strict;

my $source = "custom";
my $nstring = "n"x10000;
my $chrmpos = 1 ; 
open CONTIGS, "<$ARGV[0]" or die $!;
my $contigCount = 0; 
my $gtffile = "SuperGenome.gtf" ;
my $fastafile = "SuperGenome.fasta" ;  
my $gtflineGene ;
my $gtflineExon ; 
open GTF, ">$gtffile" or die $!;
open FASTA, ">$fastafile" or die $!; 
print FASTA ">SuperGenome\n";

my $fastaID ;
my $sequence ; 
my $ExonID ;  
while (my $y = <CONTIGS>) {
	chomp $y; 
	#if it's a >fasta header, add the previous fasta
	if ($y =~ m/^>/ ) {
		if ($contigCount == 0 ) {
			$fastaID = $y;
			$ExonID = $fastaID; 
			$gtflineGene="SuperGenome\t$source\tgene\t$chrmpos\t";
			$gtflineExon="SuperGenome\t$source\texon\t$chrmpos\t";
		}		
		else {
			$fastaID = $y ; 
			$ExonID = $fastaID;
			#position of the end of this 'gene'
			$chrmpos = length($sequence); 
			$gtflineGene .= "$chrmpos\t.\t+\t.\tgene_id \"$fastaID\";\n";
			$gtflineExon .= "$chrmpos\t.\t+\t.\tgene_id \"$fastaID\"; exon_id \"$ExonID\";\n";
			print GTF "$gtflineGene"; 
			print GTF "$gtflineExon";
			$sequence .= $nstring ; 
			#position of the start of the next gene
			$chrmpos = 1 + length($sequence); 
			$gtflineGene = "SuperGenome\t$source\tgene\t$chrmpos\t";
			$gtflineExon = "SuperGenome\t$source\texon\t$chrmpos\t";
		}
		$contigCount++;
	}
	else {
		#add the sequence to the seq variable
		$sequence .= $y ; 
	}
}
$chrmpos = length($sequence);
$gtflineGene .= "$chrmpos\t.\t+\t.\tgene_id \"$fastaID\";\n";
$gtflineExon .=	"$chrmpos\t.\t+\t.\tgene_id \"$fastaID\"; exon_id \"$ExonID\";\n";
print GTF "$gtflineGene"; 
print GTF "$gtflineExon";

while (my $chunk = substr($sequence, 0, 80, "")) {
	print FASTA "$chunk\n";
}
