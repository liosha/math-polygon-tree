package Math::Polygon::Tree;

# ABSTRACT: fast check if point is inside polygon

# $Id$

=head1 DESCRIPTION

Math::Polygon::Tree creates a tree of polygon parts for fast check if object is inside this polygon.
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
use List::MoreUtils qw{ all };
use POSIX qw/ floor ceil /;

# todo: remove gpc and use simple bbox clip
use Math::Geometry::Planar::GPC::Polygon qw{ new_gpc };



our @EXPORT_OK = qw{
    polygon_bbox
    polygon_centroid
    polygon_contains_point
};

our %EXPORT_TAGS = ( all => \@EXPORT_OK );


# tree options

our $MAX_LEAF_POINTS = 16;
our $SLICE_FIELD = 0.0001;
our $SLICE_NUM_COEF = 2;
our $SLICE_NUM_SPEED_COEF = 1;


=method new

Takes contours and creates a tree structure.
All polygons are outers, inners are not implemented.

Contour is an arrayref of points:

    my $poly1 = [ [0,0], [0,2], [2,2], ... ];   
    ...
    my $bound = Math::Polygon::Tree->new( $poly1, $poly2, ... );

or a .poly file

    my $bound1 = Math::Polygon::Tree->new( \*STDIN );
    my $bound2 = Math::Polygon::Tree->new( 'boundary.poly' );

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
    my $nparts = $SLICE_NUM_COEF * log( exp(1) * ($nrpoints/$MAX_LEAF_POINTS)**$SLICE_NUM_SPEED_COEF );

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

            my $x0 = $xmin + ($i  -$SLICE_FIELD)*$x_size;
            my $y0 = $ymin + ($j  -$SLICE_FIELD)*$y_size;
            my $x1 = $xmin + ($i+1+$SLICE_FIELD)*$x_size;
            my $y1 = $ymin + ($j+1+$SLICE_FIELD)*$y_size;

            my $gpc_slice = new_gpc();
            $gpc_slice->add_polygon([ [$x0,$y0],  [$x0,$y1], [$x1,$y1], [$x1,$y0], [$x0,$y0] ], 0);

            my @slice_parts = $gpc_poly->clip_to($gpc_slice, 'INTERSECT')->get_polygons();

            # empty part
            if ( !@slice_parts ) {
                $subparts->[$i]->[$j] = 0;
                next;
            }

            # filled part
            if (
                @slice_parts == 1 && @{$slice_parts[0]} == 4
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

    my @contours = read_poly_file( \*STDIN );
    my @contours = read_poly_file( 'bound.poly' )

Reads content of .poly-file. See http://wiki.openstreetmap.org/wiki/.poly

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
        if ( $line =~ /^([\-\!]?) (\d+)/x ) {
            $pid = $1 ? -$2 : $2;
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

    my $is_inside = $bound->contains( [1,1] );
    if ( $is_inside ) { ... }

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
        my $result = first {$_} map {polygon_contains_point($point, $_)} @{$self->{poly}};
        return $result // 0;
    }

    # branch
    my $i = min( floor( ($px-$xmin) / $self->{x_size} ), $self->{x_parts}-1 );
    my $j = min( floor( ($py-$ymin) / $self->{y_size} ), $self->{y_parts}-1 );

    my $subpart = $self->{subparts}->[$i]->[$j];
    return $subpart  if !ref $subpart;
    return $subpart->contains($point);
}


=method contains_points

    # list of points
    if ( $bound->contains_points( [1,1], [2,2] ... ) )  { ... }

    # arrayref of points
    if ( $bound->contains_points( [[1,1], [2,2] ...] ) )  { ... }

Checks if all points are inside or outside polygon.

Returns 1 if all points are inside polygon, 0 if all outside, or B<undef> otherwise.

=cut

sub contains_points {
    my ($self, @points) = @_;

    my $iter_list = @points==1 && ref $points[0]->[0]  ? $points[0]  : \@points;

    my $result;
    for my $point ( @$iter_list ) {
        my $point_result = 0 + !!$self->contains($point);
        return undef  if defined $result && $point_result != $result;
        $result = $point_result;
    }

    return $result;
}


=method contains_bbox_rough

Rough check if box is inside bound polygon.

Returns 1 if box is inside polygon, 0 if box is outside polygon or B<undef> if it 'doubts'. 

    my $bbox = [ 1, 1, 2, 2 ];
    if ( $bound->contains_bbox_rough( $bbox ) )  { ... }

=cut

sub contains_bbox_rough {
    my ($self, @bbox)  = @_;
    my $bbox = ref $bbox[0] ? $bbox[0] : \@bbox;

    croak "Box should be 4 values array: xmin, ymin, xmax, ymax" if @$bbox != 4;

    my ($x0, $y0, $x1, $y1) = @$bbox;
    my ($xmin, $ymin, $xmax, $ymax) = @{$self->{bbox}};

    # completely outside bbox
    return 0       if    $x1 < $xmin  ||  $x0 > $xmax  ||  $y0 < $ymin  ||  $y1 > $ymax;

    # partly inside
    return undef   if !( $x0 > $xmin  &&  $x1 < $xmax  &&  $y0 > $ymin  &&  $y1 < $ymax );

    return undef   if !$self->{subparts};

    # lays in defferent subparts 
    my $i0 = min( floor( ($x0-$xmin) / $self->{x_size} ), $self->{x_parts}-1 );
    my $i1 = min( floor( ($x1-$xmin) / $self->{x_size} ), $self->{x_parts}-1 );
    return undef if $i0 != $i1;
 
    my $j0 = min( floor( ($y0-$ymin) / $self->{y_size} ), $self->{y_parts}-1 );
    my $j1 = min( floor( ($y1-$ymin) / $self->{y_size} ), $self->{y_parts}-1 );
    return undef if $j0 != $j1;

    my $subpart = $self->{subparts}->[$i0]->[$j0];
    return $subpart  if !ref $subpart;
    return $subpart->contains_bbox_rough($bbox);
}


=method contains_polygon_rough

Checks if polygon is inside bound polygon.

Returns 1 if inside, 0 if outside or B<undef> if 'doubts'. 

    if ( $bound->contains_polygon_rough( [ [1,1], [1,2], [2,2], ... ] ) )  { ... }

=cut

sub contains_polygon_rough {
    my ($self, $poly) = @_; 
    croak "Polygon should be a reference to array of points" if ref $poly ne 'ARRAY';

    return $self->contains_bbox_rough( polygon_bbox($poly) );
}


=method bbox

    my $bbox = $bound->bbox();
    my ($xmin, $ymin, $xmax, $ymax) = @$bbox;

Returns polygon's bounding box.

=cut

sub bbox {
    return shift()->{bbox};
}



=func polygon_bbox

    my $bbox = polygon_bbox( [[1,1], [1,2], [2,2], ... ] );
    my ($xmin, $ymin, $xmax, $ymax) = @$bbox;

Returns polygon's bounding box.

=cut

sub polygon_bbox {
    my ($contour) = @_;

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

    my $center_point = polygon_centroid( [ [1,1], [1,2], [2,2], ... ] );

Returns polygon's weightened center.

Math::Polygon 1.02+ has the same function, but it is very inaccurate.

=cut

sub polygon_centroid {
    my (@poly) = @_;
    my $contour = ref $poly[0]->[0] ? $poly[0] : \@poly;

    return $contour->[0]  if @$contour < 2;

    my $sx = 0;
    my $sy = 0;
    my $sq = 0;

    my $p0 = $contour->[0];
    for my $i ( 1 .. $#$contour-1 ) {
        my $p  = $contour->[$i];
        my $p1 = $contour->[$i+1];

        my $tsq = ( ( $p->[0]  - $p0->[0] ) * ( $p1->[1] - $p0->[1] )
                  - ( $p1->[0] - $p0->[0] ) * ( $p->[1]  - $p0->[1] ) );
        next if $tsq == 0;
        
        my $tx = ( $p0->[0] + $p->[0] + $p1->[0] ) / 3;
        my $ty = ( $p0->[1] + $p->[1] + $p1->[1] ) / 3;

        $sx += $tx * $tsq;
        $sy += $ty * $tsq;
        $sq += $tsq;
    }

    if ( $sq == 0 ) {
        my $bbox = polygon_bbox($contour);
        return [ ($bbox->[0]+$bbox->[2])/2, ($bbox->[1]+$bbox->[3])/2 ];
    }

    return [$sx/$sq, $sy/$sq];
}


=func polygon_contains_point

    my $is_inside = polygon_contains_point($point, $polygon);

Function that tests if polygon contains point (modified one from Math::Polygon::Calc).

Returns -1 if point lays on polygon's boundary

=cut

sub polygon_contains_point {
    my ($point, @poly) = @_;
    my $contour = ref $poly[0]->[0] ? $poly[0] : \@poly;

    my ($x, $y) = @$point;
    my ($px, $py) = @{ $contour->[0] };
    my ($nx, $ny);

    my $inside = 0;

    for my $i ( 1 .. scalar @$contour ) { 
        ($nx, $ny) =  @{ $contour->[ $i % scalar @$contour ] };

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

