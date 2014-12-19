#!/usr/bin/perl
use strict;

my ( $file ) = @ARGV;
my %dmr;

open I, $file or die "Cannot open file : $file\n";
while ( <I> ) {
	chomp;
	my @tmp = split;
	my $key = join ":",( @tmp[0..2] );
	if ( ! defined $dmr{ $key } ) {
		@{$dmr{ $key }} = @tmp;
	} else {
		@{$dmr{ $key }} = @tmp if $dmr{ $key }[-3] < $tmp[-3];
	}
}
close ( I );

for my $key ( keys %dmr ) {

	print join "\t", ( @{$dmr{ $key }} ); print "\n";
}
