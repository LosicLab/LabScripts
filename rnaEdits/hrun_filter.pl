#!/usr/bin/perl
use strict;
use warnings;

my $limit = 4;

while (my $line = <STDIN>) {
	chomp $line;
	if ($line =~ m/^#/ ) { print "$line\n"; next;}
	my ($trash, $keep ) = split (/HRun=/, $line); 
	my ($hrun, $trash2 ) = split (/;/, $keep);
	if ($hrun <= $limit) {
		print "$line\n";
	}
}
