#!/usr/bin/perl

use 5.010;
use strict;
use warnings;
use utf8;

use Test::More;

use Math::Polygon::Tree qw/ polygon_centroid /;

our @TESTS = (
    [ [[0,0],[2,4]],             [1,2], 'single segment' ],
    [ [[0,0],[3,0],[0,3]],       [1,1], 'triangle' ],
    [ [[0,0],[4,0],[4,4],[0,4]], [2,2], 'square' ],
);


for my $test ( @TESTS ) {
    my ($in, $expected, $name) = @$test;
    my $got = polygon_centroid($in);
    is_deeply( $got, $expected, $name );
}


done_testing();



