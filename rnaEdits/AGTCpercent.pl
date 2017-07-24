#!/usr/bin/perl
use strict;
use warnings;
#input (via stdin) should be the result of something like :
# cut -f1,2,4,5 myfiles*.vcf | sort | uniq | cut -f3,4 |sort |uniq -c 

my $sum =0;
my $agtcsum = 0;
my $ag = 0;
my $ct = 0; 
while (my $x = <STDIN>) {
	chomp $x;
	next if ($x =~ /#/) ; #should remove headers 
	my ($nothing, $count, $g1, $g2) = split(/\s+/, $x);
	#print "G1:$g1\tG2:$g2\tCount:$count\n";
	$sum +=$count;
	if ($g1 eq "A") {
		if ($g2 eq "G") {
			$agtcsum+=$count;
			$ag+=$count;
		}
	}
	if ($g1 eq "T")	{
                if ($g2 eq "C")	{
                        $agtcsum+=$count;
			$ct+=$count;
                }
        }
}
my $fraction = $agtcsum/$sum; 
print "Total:$sum Edits:$agtcsum A>G:$ag C>T:$ct Percent:$fraction\n";
