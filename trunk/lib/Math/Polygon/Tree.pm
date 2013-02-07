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

use List::Util qw{ reduce first sum min max };
use List::MoreUtils qw{ all uniq };
use POSIX qw/ floor ceil /;

# FIXME: remove and use simple bbox clip?
use Math::Geometry::Planar::GPC::Polygon qw{ new_gpc };



our @EXPORT_OK = qw{
    polygon_bbox
    polygon_centroid
    polygon_contains_point
};



our $MAX_LEAF_POINTS = 16;
our $SLICE_COEF = 2;


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
    my ($class, @in_contours) = @_;
    my $self = bless {}, $class;

    ##  load and close polys, calc bbox
    my @contours;
    while ( @in_contours ) {
        my $contour = shift @in_contours;

        if ( ref $contour ne 'ARRAY' ) {
            unshift @in_contours, read_poly_file($contour);
            next;
        }

        my @points = @$contour;
        push @points, $points[0]  if !( $points[0] ~~ $points[-1] );

        push @contours, \@points;

        my $bbox = polygon_bbox(\@points);
        $self->{bbox} = bbox_union($bbox, $self->{bbox});
    }

    croak "No contours"  if !@contours;

    my $nrpoints = sum map { scalar @$_ } @contours;


    # small polygon - no need to slice
    if ( $nrpoints <= $MAX_LEAF_POINTS ) {
        $self->{poly} = \@contours;
        return $self;
    }


    # calc number of pieces (need to tune!)
    my ($xmin, $ymin, $xmax, $ymax) = @{$self->{bbox}};
    my $xy_ratio = ($xmax-$xmin) / ($ymax-$ymin);
    my $nparts = $SLICE_COEF * log( exp(1) * ($nrpoints/$MAX_LEAF_POINTS) );

    my $x_parts = $self->{x_parts} = ceil( sqrt($nparts * $xy_ratio) );
    my $y_parts = $self->{y_parts} = ceil( sqrt($nparts / $xy_ratio) );
    my $x_size  = $self->{x_size}  = ($xmax-$xmin) / $x_parts;
    my $y_size  = $self->{y_size}  = ($ymax-$ymin) / $y_parts;


    # slice
    my $subparts = $self->{subparts} = [];
    
    my $gpc_poly = new_gpc();
    $gpc_poly->add_polygon( $_, 0 )  for @contours;
    
    for my $j ( 0 .. $y_parts-1 ) {
        for my $i ( 0 .. $x_parts ) {

            my $x0 = $xmin + ($i-0.0001)*$x_size;
            my $y0 = $ymin + ($j-0.0001)*$y_size;
            my $x1 = $xmin + ($i+1.0001)*$x_size;
            my $y1 = $ymin + ($j+1.0001)*$y_size;

            my $gpc_slice = new_gpc();
            $gpc_slice->add_polygon([ [$x0,$y0],  [$x0,$y1], [$x1,$y1], [$x1,$y0], [$x0,$y0] ], 0);

            my @slice_parts = $gpc_poly->clip_to($gpc_slice, 'INTERSECT')->get_polygons();

            # empty part
            if ( !@slice_parts ) {
                $subparts->[$i]->[$j] = 0;
                next;
            }

            # filled part
            if ( @slice_parts == 1 && @{$slice_parts[0]} == 4
                && all { $_->[0] ~~ [$x0,$x1] && $_->[1] ~~ [$y0,$y1] } @{$slice_parts[0]}
            ) {
                $subparts->[$i]->[$j] = 1;
                next;
            }

            # complex subpart
            $subparts->[$i]->[$j] = Math::Polygon::Tree->new(@slice_parts);
        }
    }

    return $self;
}


=func read_poly_file

    my @contours = read_poly_file( \*STDIN )

Reads content of .poly-file

=cut

sub read_poly_file {
    my ($file) = @_;

    my $need_to_open = !ref $file || ref $file eq 'SCALAR';
    my $fh = $need_to_open
        ? do { open my $in, '<', $file  or croak "Couldn't open $file: $@"; $in }
        : $file;

    my @contours;
    my $pid;
    my @cur_points;
    while ( my $line = readline $fh ) {
        # new contour
        if ( $line =~ /^(-?\d+)/ ) {
            $pid = $1;
            next;
        }

        # point
        if ( $line =~ /^\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)/ ) {
            push @cur_points, [ $1+0, $2+0 ];
            next;
        }

        # !!! inner contour - skipping
        if ( $line =~ /^END/  &&  $pid < 0 ) {
            @cur_points = ();
            next;
        }

        # outer contour
        if ( $line =~ /^END/  &&  @cur_points ) {
            push @contours, [ @cur_points ];
            @cur_points = ();
            next;
        }
    }

    close $fh  if $need_to_open;
    return @contours;
}


=method contains

    if ( $bound->contains( [1,1] ) )  { ... }

Checks if point is inside bound polygon.

Returns 1 if point is inside polygon, -1 if it lays on polygon boundary, or 0 otherwise.

=cut

sub contains {
    my ($self, $point) = @_;

    croak "Point should be a reference"  if ref $point ne 'ARRAY';

    # check bbox
    my ($px, $py) = @$point;
    my ($xmin, $ymin, $xmax, $ymax) = @{ $self->{bbox} };
    return 0  if $px < $xmin  ||  $px > $xmax  ||  $py < $ymin  ||  $py > $ymax;

    # leaf
    if ( exists $self->{poly} ) {
        my $result = first {$_} map {polygon_contains_point($point, @$_)} @{$self->{poly}};
        return $result // 0;
    }

    # branch
    my $i = floor( ($px-$xmin) / $self->{x_size} );
    my $j = floor( ($py-$ymin) / $self->{y_size} );

    my $subpart = $self->{subparts}->[$i]->[$j];
    return $subpart  if !ref $subpart;
    return $subpart->contains($point);
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

    my $bbox = polygon_bbox( $poly );
    return $self->contains_bbox_rough( @$bbox );
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

    my $bbox = polygon_bbox( [[1,1], [1,2], [2,2], ... ] );
    my ($xmin, $ymin, $xmax, $ymax) = @$bbox;

Returns polygon's bounding box.

=cut

sub polygon_bbox {
    my ($contour) = @_;

    no warnings 'once';
    return bbox_union(@$contour) if @$contour <= 2;
    return reduce { bbox_union($a, $b) } @$contour;
}


=func bbox_union

    my $united_bbox = bbox_union($bbox1, $bbox2);

Returns united bbox for two bboxes/points.

=cut

sub bbox_union {
    my ($bbox1, $bbox2) = @_;

    $bbox2 //= $bbox1;

    my @bbox = (
        min( $bbox1->[0], $bbox2->[0] ),
        min( $bbox1->[1], $bbox2->[1] ),
        max( $bbox1->[2] // $bbox1->[0], $bbox2->[2] // $bbox2->[0] ),
        max( $bbox1->[3] // $bbox1->[1], $bbox2->[3] // $bbox2->[1] ),
    );

    return \@bbox;
}


=func polygon_centroid

Function that returns polygon's weightened center.

    my ( $x, $y ) = polygon_centroid( [1,1], [1,2], [2,2], ... );

=cut

sub polygon_centroid {

    my $slat = 0;
    my $slon = 0;
    my $ssq  = 0;

    for my $i ( 1 .. $#_-1 ) {
        my $tlon = ( $_[0]->[0] + $_[$i]->[0] + $_[$i+1]->[0] ) / 3;
        my $tlat = ( $_[0]->[1] + $_[$i]->[1] + $_[$i+1]->[1] ) / 3;

        my $tsq = ( ( $_[$i]  ->[0] - $_[0]->[0] ) * ( $_[$i+1]->[1] - $_[0]->[1] )
                  - ( $_[$i+1]->[0] - $_[0]->[0] ) * ( $_[$i]  ->[1] - $_[0]->[1] ) );

        $slat += $tlat * $tsq;
        $slon += $tlon * $tsq;
        $ssq  += $tsq;
    }

    if ( $ssq == 0 ) {
        return (
            ((min map { $_->[0] } @_) + (max map { $_->[0] } @_)) / 2,
            ((min map { $_->[1] } @_) + (max map { $_->[1] } @_)) / 2 );
    }

    return ( $slon/$ssq , $slat/$ssq );
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

