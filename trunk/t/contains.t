#!/usr/bin/perl

use 5.010;
use strict;
use warnings;
use utf8;

use Test::More;

use Math::Polygon::Tree;

my $t = Math::Polygon::Tree->new([[0,0],[2,0],[0,2]]);

our @TESTS = (
    [ [0,0], -1, 'vertex point' ],
    [ [0,1], -1, 'edge point 1' ],
    [ [1,1], -1, 'edge point 2' ],
    [ [0.5,0.5], 1, 'inside' ],
    [ [1.5,1.5], 0, 'outside' ],
    [ [8,8], 0, 'far outside' ],
);


for my $test ( @TESTS ) {
    my ($in, $expected, $name) = @$test;
    my $got = $t->contains($in);
    is( $got, $expected, $name );
}


done_testing();



