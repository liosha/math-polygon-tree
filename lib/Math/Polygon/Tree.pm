package Math::Polygon::Tree;

# ABSTRACT: fast check if point is inside polygon

# $Id$

=head1 DESCRIPTION

Math::Polygon::Tree creates a B-tree of polygon parts for fast check if object is inside this polygon.
This method is effective if polygon has hundreds or more segments.

=head1 SYNOPSIS

    use Math::Polygon::Tree;

    my $poly  = [ [0,0], [0,2], [2,2], ... ];
    my $bound = Math::Polygon::Tree->new( $poly );

    if ( $bound->contains( [1,1] ) )  { ... }

=cut

use 5.010;
use strict;
use warnings;
use utf8;
use Carp;

use base qw{ Exporter };

use List::Util qw{ sum min max };
use List::MoreUtils qw{ uniq };

# FIXME: remove and use simple bbox clip?
use Math::Geometry::Planar::GPC::Polygon qw{ new_gpc };



our @EXPORT_OK = qw{
    polygon_bbox
    polygon_contains_point
};



my $MAX_LEAF_POINTS = 16;       # minimum 6


=method new

Takes [at least one] contour and creates a tree structure. All polygons are outer, inners in not implemented.

Contour is an arrayref of points:

    my $poly1 = [ [0,0], [0,2], [2,2], ... ];   
    ...
    my $bound = Math::Polygon::Tree->new( $poly1, $poly2, ... );

or a .poly file

    my $bound = Math::Polygon::Tree->new( 'boundary.poly' );

=cut

sub new {
    my $class = shift;
    my $self  = {};

    ##  load and close polys, calc bbox
    while ( my $chain_ref = shift ) {

        if ( ref $chain_ref ) {
            croak "Polygon should be a reference to array of points" 
                unless ref $chain_ref eq 'ARRAY';
        
            my @epoint = ();
            push @epoint, $chain_ref->[0]
                unless  $chain_ref->[0]->[0] == $chain_ref->[-1]->[0]
                    &&  $chain_ref->[0]->[1] == $chain_ref->[-1]->[1];

            my $poly = [ @$chain_ref, @epoint ];
            push @{$self->{poly}}, $poly;

            my ($xmin, $ymin, $xmax, $ymax) = polygon_bbox( @$poly );

            $self->{xmin} = $xmin       if  !(exists $self->{xmin})  ||  $xmin < $self->{xmin};
            $self->{xmax} = $xmax       if  !(exists $self->{xmax})  ||  $xmax > $self->{xmax};
            $self->{ymin} = $ymin       if  !(exists $self->{ymin})  ||  $ymin < $self->{ymin};
            $self->{ymax} = $ymax       if  !(exists $self->{ymax})  ||  $ymax > $self->{ymax};
        }
        else {

            open my $in, '<', $chain_ref
                or croak "Couldn't open $chain_ref: $!";

            my @bound;
            my $pid;

            while ( my $line = readline $in ) {
                if ( $line =~ /^(-?\d+)/ ) {
                    $pid = $1;
                }
                elsif ( $line =~ /^\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)/ ) {
                    push @bound, [ $1+0, $2+0 ];
                }
                elsif ( $line =~ /^END/  &&  $pid < 0 ) {
                    @bound = ();
                }
                elsif ( $line =~ /^END/  &&  @bound ) {

                    push @bound, $bound[0] 
                        unless  $bound[0]->[0] == $bound[-1]->[0]
                            &&  $bound[0]->[1] == $bound[-1]->[1];
                    push @{$self->{poly}}, [ @bound ];

                    my ($xmin, $ymin, $xmax, $ymax) = polygon_bbox( @bound );

                    $self->{xmin} = $xmin       if  !(exists $self->{xmin})  ||  $xmin < $self->{xmin};
                    $self->{xmax} = $xmax       if  !(exists $self->{xmax})  ||  $xmax > $self->{xmax};
                    $self->{ymin} = $ymin       if  !(exists $self->{ymin})  ||  $ymin < $self->{ymin};
                    $self->{ymax} = $ymax       if  !(exists $self->{ymax})  ||  $ymax > $self->{ymax};

                    @bound = ();
                }
            }

            close $in;
        }
    }

    my $nrpoints = sum map { scalar @$_ } @{$self->{poly}};

    ##  full square?
    if ( $nrpoints == 5 ) {
        my $poly = $self->{poly}->[0];
        my @xs = uniq map { $_->[0] } @$poly;
        my @ys = uniq map { $_->[1] } @$poly;

        if ( @xs == 2  &&  @ys == 2 ) {
            $self->{full} = 1;
        }
    }

    ##  branch if big poly
    if ( $nrpoints > $MAX_LEAF_POINTS ) {
        # 0 - horisontal split, 1 - vertical
        $self->{hv}  =  my $hv  =  ($self->{xmax}-$self->{xmin}) < ($self->{ymax}-$self->{ymin});
        $self->{avg} =  my $avg =  $hv  ?  ($self->{ymax}+$self->{ymin})/2  :  ($self->{xmax}+$self->{xmin})/2;

        my $gpc = new_gpc();
        for my $poly ( @{$self->{poly}} ) {
            $gpc->add_polygon( $poly, 0 );
        }

        my $part1_gpc = new_gpc();
        $part1_gpc->add_polygon( [
                [ $self->{xmin}, $self->{ymin} ],
                $hv  ?  [ $self->{xmin}, $avg          ]  :  [ $self->{xmin}, $self->{ymax} ],
                $hv  ?  [ $self->{xmax}, $avg          ]  :  [ $avg         , $self->{ymax} ],
                $hv  ?  [ $self->{xmax}, $self->{ymin} ]  :  [ $avg         , $self->{ymin} ],
                [ $self->{xmin}, $self->{ymin} ],
            ], 0 );
        $part1_gpc = $gpc->clip_to( $part1_gpc, 'INTERSECT' );
        $self->{part1} = Math::Polygon::Tree->new( $part1_gpc->get_polygons() );

        my $part2_gpc = new_gpc();
        $part2_gpc->add_polygon( [
                [ $self->{xmax}, $self->{ymax} ],
                $hv  ?  [ $self->{xmax}, $avg          ]  :  [ $self->{xmax}, $self->{ymin} ],
                $hv  ?  [ $self->{xmin}, $avg          ]  :  [ $avg         , $self->{ymin} ],
                $hv  ?  [ $self->{xmin}, $self->{ymax} ]  :  [ $avg         , $self->{ymax} ],
                [ $self->{xmax}, $self->{ymax} ],
            ], 0 );
        $part2_gpc = $gpc->clip_to( $part2_gpc, 'INTERSECT' );
        $self->{part2} = Math::Polygon::Tree->new( $part2_gpc->get_polygons() );

        delete $self->{poly};
    }

    bless ($self, $class);
    return $self;
}


=method contains

    if ( $bound->contains( [1,1] ) )  { ... }

Checks if point is inside bound polygon.

Returns 1 if point is inside polygon, -1 if it lays on polygon boundary (dirty), or 0 otherwise.

=cut

sub contains {
    my $self  = shift;
    my $point = shift;
    croak "Point should be a reference" 
        unless ref $point;

    my ($px, $py) = @$point;

    return 0
        if      $px < $self->{xmin}  ||  $px > $self->{xmax}
            ||  $py < $self->{ymin}  ||  $py > $self->{ymax};

    return $self->{full}    if  exists $self->{full};

    if ( exists $self->{hv} ) {
        if ( $point->[$self->{hv}] < $self->{avg} ) {
            return $self->{part1}->contains( $point );
        }
        else {
            return $self->{part2}->contains( $point );
        }
    }

    if ( exists $self->{poly} ) {
        for my $poly ( @{$self->{poly}} ) {
            return polygon_contains_point( $point, @$poly );
        }
    }

    return 0;
}


=method contains_points

Checks if points are inside bound polygon.

Returns 1 if all points are inside polygon, 0 if all outside, or B<undef>.

    if ( $bound->contains_points( [1,1], [2,2] ... ) )  { ...

=cut

sub contains_points {
    my $self  = shift;
    my $result = undef;
    
    while ( my $point = shift ) {
        next unless ref $point;

        my $isin = abs $self->contains( $point );
        if ( defined $result ) {
            return undef  unless  $isin == $result;
        }
        else {
            $result = $isin;
        }
    }

    return $result;
}


=method contains_bbox_rough

Checks if box is inside bound polygon.

Returns 1 if box is inside polygon, 0 if box is outside polygon or B<undef> if it 'doubts'. 

    my ($xmin, $ymin, $xmax, $ymax) = ( 1, 1, 2, 2 );
    if ( $bound->contains_bbox_rough( $xmin, $ymin, $xmax, $ymax ) )  { ... }

=cut

sub contains_bbox_rough {
    my $self  = shift;
    croak "Box should be 4 values xmin, ymin, xmax, ymax"
        unless @_ == 4;

    my ($xmin, $ymin, $xmax, $ymax) = @_;

    return 0
        if   $xmax < $self->{xmin}  ||  $xmin > $self->{xmax}
         ||  $ymax < $self->{ymin}  ||  $ymin > $self->{ymax};

    if (  $xmin > $self->{xmin}  &&  $xmax < $self->{xmax}
      &&  $ymin > $self->{ymin}  &&  $ymax < $self->{ymax} ) {

        return $self->{full}    if  exists $self->{full};

        if ( exists $self->{hv} ) {
            if ( $self->{hv} ) {
                return $self->{part1}->contains_bbox_rough( @_ )
                    if  $ymax < $self->{avg}; 
                return $self->{part2}->contains_bbox_rough( @_ )
                    if  $ymin > $self->{avg}; 
            }
            else {
                return $self->{part1}->contains_bbox_rough( @_ )
                    if  $xmax < $self->{avg}; 
                return $self->{part2}->contains_bbox_rough( @_ )
                    if  $xmin > $self->{avg}; 
            }
        }
    }

    return undef;     
}


=method contains_polygon_rough

Checks if polygon is inside bound polygon.

Returns 1 if inside, 0 if outside or B<undef> if 'doubts'. 

    if ( $bound->contains_polygon_rough( [ [1,1], [1,2], [2,2], ... ] ) )  { ... }

=cut

sub contains_polygon_rough {
    my $self = shift;
    my $poly = shift; 

    croak "Polygon should be a reference to array of points" 
        unless ref $poly;

    my @bbox = polygon_bbox( @$poly );
    return $self->contains_bbox_rough( @bbox );
}



=method bbox

Returns polygon's bounding box. 

    my ( $xmin, $ymin, $xmax, $ymax ) = $bound->bbox();

=cut

sub bbox {
    my $self  = shift;
    return ( $self->{xmin}, $self->{ymin}, $self->{xmax}, $self->{ymax} );
}




=func polygon_bbox

Function that returns polygon's bbox.

    my ( $xmin, $ymin, $xmax, $ymax ) = polygon_bbox( [1,1], [1,2], [2,2], ... );

=cut

sub polygon_bbox (@) {

    return (
        ( min map { $_->[0] } @_ ),
        ( min map { $_->[1] } @_ ),
        ( max map { $_->[0] } @_ ),
        ( max map { $_->[1] } @_ ),
    );
}



=func polygon_contains_point

Function that tests if polygon contains point (modified one from Math::Polygon::Calc).

Returns -1 if point lays on polygon's boundary

=cut

sub polygon_contains_point ($@) {

    my $point = shift;

    my ( $x,  $y)  =  @$point;
    my ($px, $py)  =  @{ (shift) };
    my ($nx, $ny);

    my $inside = 0;

    while( @_ ) {
        ($nx, $ny) =  @{ (shift) };
        
        return -1
            if  $y == $py  &&  $py == $ny
                && ( $x >= $px  ||  $x >= $nx )
                && ( $x <= $px  ||  $x <= $nx );

        next    if  $py == $ny;
        next    if  $y < $py  &&  $y < $ny;
        next    if  $y > $py  &&  $y > $ny;
        next    if  $x > $px  &&  $x > $nx;

        my $xx = ($y-$py)*($nx-$px)/($ny-$py)+$px;
        return -1   if  $x == $xx;

        next    if  $y <= $py  &&  $y <= $ny;

        $inside = 1 - $inside
            if  $px == $nx  ||  $x < $xx;
    }
    continue { ($px, $py) = ($nx, $ny); }

    return $inside;
}


1;

