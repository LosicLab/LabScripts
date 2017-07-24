#!/usr/bin/perl
use warnings;
use strict;

#usage: vcf_parser.pl myfile.vcf > output
### my current usage involves prefiltering VCFs for lines with a LoF tag. 

#read through VCF
open VCF, "<$ARGV[0]" or die $!;
while (my $x = <VCF>) {
	chomp $x; 
	if ($x ~~ /^#/) { print "$x\n"; next ; }
	my @tabsplit= split (/\s+/, $x);# 0:chr 1pos 2id 3gen1 4gen2 5score 6pass/np 7annotation 8genoinfo 9genovalues
	##plan: 3 reads, quality over 25, PASS)
	if ($tabsplit[6] eq "PASS" ) {
		if ($tabsplit[5] >=25 ) {
			my $allele1depth = 0;
			my $allele2depth = 0; 
			foreach my $counts (9..scalar(@tabsplit)-1) {
				my $readcount = $tabsplit[$counts]; 
				if ($readcount ~~ /:/) {
					my ($genotype, $alleleDepth) = split(/:/, $readcount); #count values ie 0/1:7,103:110:99:4305,0,175
					my ($allele1depth_temp, $allele2depth_temp) = split( /,/, $alleleDepth);
					$allele1depth += $allele1depth_temp;
					$allele2depth += $allele2depth_temp;
				}	
			}

			#my $readcounts = $tabsplit[-1]; #count values ie 0/1:7,103:110:99:4305,0,175
			#my ($genotype, $alleleDepth) = split(/:/, $readcounts);
			#my ($allele1depth, $allele2depth) = split( /,/, $alleleDepth);
			if ($allele1depth >=3 && $allele2depth >=3 ) {
				print "$x\n";	
			}
		}
	}
}

