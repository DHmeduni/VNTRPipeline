#!/usr/bin/env perl
use strict;
use warnings;
use Text::Fuzzy;

my ($start, $end, $max_mismatches) = @ARGV;
my $query = do { local $/; <STDIN> };   # read sequence
$query =~ s/\s+//g;                     # strip newlines

my $tf_start = Text::Fuzzy->new($start);
my $tf_end   = Text::Fuzzy->new($end);

for my $i (0 .. length($query) - length($start)) {
    my $sub = substr($query, $i, length($start));
    if ($tf_start->distance($sub) <= $max_mismatches) {
        # found approximate start
        for my $j ($i + length($start) .. length($query) - length($end)) {
            my $sub_end = substr($query, $j, length($end));
            if ($tf_end->distance($sub_end) <= $max_mismatches) {
                my $match = substr($query, $i, $j + length($end) - $i);
                print "$match\n";
                last;
            }
        }
    }
}

