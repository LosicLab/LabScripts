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
my $ChromosomeIndex = 1; 
my $gtflineGene ;
my $gtflineExon ; 
open GTF, ">$gtffile" or die $!;
open FASTA, ">$fastafile" or die $!; 
print FASTA ">SuperGenome_$ChromosomeIndex\n";

my $fastaID ;
my $sequence ; 
my $Geneus ; 
my $Species ; 
my %datahash ; # format: datahash{genus}[index][0=sequence, 1=species, 2=size]
my %geneuscount ; 
while (my $y = <CONTIGS>) {
	chomp $y; 
	#if it's a >fasta header, add the previous fasta
	if ($y =~ m/^>/ ) {
		if ($contigCount == 0 ) { #initiate
			$fastaID = $y;
			$fastaID =~ s/\'//g; 
			my $templine = $fastaID ;
			# Working with Fasta headers that look like:
			#>gi|56899872|ref|NC_006578.1| Bacillus thuringiensis serovar konkukian str. 97-27 plasmid pBT9727, complete sequence 
			$templine =~ s/.*\| //; #remove everything before the |
 			$Species = $templine; 
			my @templinesplit = split (/\s/, $templine);
			$Geneus = $templinesplit[0];
		}		
		else {
			#position of the end of this 'gene'
			$geneuscount{$Geneus}++; 
			$datahash{$Geneus}[$geneuscount{$Geneus}][0] = $sequence ; 
			$datahash{$Geneus}[$geneuscount{$Geneus}][1] = $Species ;  
			$datahash{$Geneus}[$geneuscount{$Geneus}][2] = length($sequence) ;  
			#Now start the next sequence
			$fastaID = $y ; 
			$fastaID =~ s/\'//g;
			my $templine = $fastaID ;
			$templine =~ s/.*\| //; #remove everything before the |
                        $Species = $templine; 
			my @templinesplit = split (/\s/, $templine);
                        $Geneus = $templinesplit[0];
			$sequence = ""; 
			#position of the start of the next gene
		}
		$contigCount++;
	}
	else {
		#add the sequence to the seq variable
		$sequence .= $y ; 
	}
}
$geneuscount{$Geneus}++;
$datahash{$Geneus}[$geneuscount{$Geneus}][0] = $sequence ;
$datahash{$Geneus}[$geneuscount{$Geneus}][1] = $Species ;
$datahash{$Geneus}[$geneuscount{$Geneus}][2] = length($sequence) ;

my %geneSize ;
my $chrmpos1 ;  
my $fullsequence; 
foreach my $geneus (sort keys %datahash) {
	$chrmpos1 = $chrmpos ; 
	#Print the Gene Line
	print GTF "SuperGenome_$ChromosomeIndex\t$source\tgene\t$chrmpos\t";
	foreach my $index ( 1..$geneuscount{$geneus}) { #go through all species in gene, add up size/sequence.
		$geneSize{$geneus} += $datahash{$geneus}[$index][2] ; 
		$geneSize{$geneus} += length($nstring) ; 
		$fullsequence .= $datahash{$geneus}[$index][0];
		$fullsequence .= $nstring ;
	}
	$chrmpos += $geneSize{$geneus} - 1 ; 
	print GTF "$chrmpos\t.\t+\t.\tgene_id \"$geneus\";\n"; 
	$chrmpos += 1; 
	my $chrmpos2 ; 
	#Print the Exon Lines
	foreach my $index ( 1..$geneuscount{$geneus}) {
		$chrmpos2 = $chrmpos1 + length($datahash{$geneus}[$index][0]) -1; 
		print GTF "SuperGenome_$ChromosomeIndex\t$source\texon\t$chrmpos1\t$chrmpos2\t.\t+\t.\tgene_id \"$geneus\"; exon_id \"$datahash{$geneus}[$index][1]\";\n"; 
		$chrmpos1 = $chrmpos2 + length($nstring) + 1 ;
	}
	if (length($fullsequence) > 200000000 ) { #approaching the size limit for indexing.
		while (my $chunk = substr($fullsequence, 0, 80, "")) {
	        	print FASTA "$chunk\n";
		}
		$chrmpos = 1; 
		$ChromosomeIndex++;
		$fullsequence="";
		print FASTA ">SuperGenome_$ChromosomeIndex\n";
	}
}
while (my $chunk = substr($fullsequence, 0, 80, "")) {
	print FASTA "$chunk\n";
}
#system('rm SuperGenome.fata.temp');
#print out the fasta 80 bases per line
#while (my $chunk = substr($fullsequence, 0, 80, "")) {
#	print FASTA "$chunk\n";
#}
