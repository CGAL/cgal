// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_KERNEL_PREDICATE_OBJECTS_2_H
#define CGAL_KERNEL_PREDICATE_OBJECTS_2_H

#include <CGAL/predicate_objects_on_points_2.h>
#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

template < class R >
struct Kernel_predicate_objects_2
{
typedef typename R::FT                      FT;
typedef typename R::RT                      RT;
typedef typename R::Point_2                 Point_2;
typedef typename R::Vector_2                Vector_2;
typedef typename R::Direction_2             Direction_2;
typedef typename R::Segment_2               Segment_2;
typedef typename R::Line_2                  Line_2;
typedef typename R::Ray_2                   Ray_2;
typedef typename R::Triangle_2              Triangle_2;
typedef typename R::Circle_2                Circle_2;
typedef typename R::Iso_rectangle_2         Iso_rectangle_2;
typedef typename R::Aff_transformation_2    Aff_transformation_2;

/*
typedef CGALi::std::equal_to                                  Equal_2;
Equal_2
equal_2_object() const { return Equal_2(); }

typedef CGALi::p_Equal_xy                              Equal_xy_2;
Equal_xy_2
equal_xy_2_object() const { return Equal_xy_2(); }
*/

typedef p_Less_xy<Point_2>                       Less_xy_2;
Less_xy_2
less_xy_2_object() const { return Less_xy_2(); }

/*
typedef CGALi::p_Less_dist_to_point<Point_2>           Less_distance_to_point_2;
Less_distance_to_point_2
less_distance_to_point_2_object() const { return Less_distance_to_point_2(); }
*/

typedef p_Less_dist_to_line_2p<Point_2>         Less_signed_distance_to_line_2;
Less_signed_distance_to_line_2
less_signed_distance_to_line_2_object() const { return Less_signed_distance_to_line_2(); }

typedef p_Leftturn<Point_2>                     Leftturn_2;
Leftturn_2
leftturn_2_object() const { return Leftturn_2(); }

typedef p_Left_of_line_2p<Point_2>              Left_of_line_2;
Left_of_line_2
left_of_line_2_object() const { return Left_of_line_2(); }

typedef CGALi::Collinear                               Collinear_2;
Collinear_2
collinear_2_object() const { return Collinear_2(); }

typedef p_Orientation<Point_2>                  Orientation_2;
Orientation_2
orientation_2_object() const { return Orientation_2(); }

typedef CGALi::Side_of_oriented_circle                 Side_of_oriented_circle_2;
Side_of_oriented_circle_2
side_of_oriented_circle_2_object() const { return Side_of_oriented_circle_2(); }

typedef CGALi::Call_is_horizontal                      Is_horizontal_2;
Is_horizontal_2
is_horizontal_2_object() const { return Is_horizontal_2(); }

typedef CGALi::Call_is_vertical                        Is_vertical_2;
Is_vertical_2
is_vertical_2_object() const { return Is_vertical_2(); }

typedef CGALi::Call_is_degenerate                      Is_degenerate_2;
Is_degenerate_2
is_degenerate_2_object() const { return Is_degenerate_2(); }

typedef CGALi::Call_has_on_bounded_side                Has_on_bounded_side_2;
Has_on_bounded_side_2
has_on_bounded_side_2_object() const { return Has_on_bounded_side_2(); }

typedef CGALi::Call_has_on_unbounded_side              Has_on_unbounded_side_2;
Has_on_unbounded_side_2
has_on_unbounded_side_2_object() const { return Has_on_unbounded_side_2(); }

typedef CGALi::Call_has_on_boundary                    Has_on_boundary_2;
Has_on_boundary_2
has_on_boundary_2_object() const { return Has_on_boundary_2(); }

typedef CGALi::Call_has_on_positive_side               Has_on_positive_side_2;
Has_on_positive_side_2
has_on_positive_side_2_object() const { return Has_on_positive_side_2(); }

typedef CGALi::Call_has_on_negative_side               Has_on_negative_side_2;
Has_on_negative_side_2
has_on_negative_side_2_object() const { return Has_on_negative_side_2(); }

typedef CGALi::Call_oriented_side                      Oriented_side_2;
Oriented_side_2
oriented_side_2_object() const { return Oriented_side_2(); }

typedef CGALi::Compare_x                               Compare_x_2;
Compare_x_2
compare_x_2_object() const { return Compare_x_2(); }

typedef CGALi::Compare_y                               Compare_y_2;
Compare_y_2
compare_y_2_object() const { return Compare_y_2(); }

typedef CGALi::Compare_y_at_x                          Compare_y_at_x_2;
Compare_y_at_x_2
compare_y_at_x_2_object() const { return Compare_y_at_x_2(); }

typedef CGALi::Are_ordered_along_line                  Are_ordered_along_line_2 ;
Are_ordered_along_line_2
are_ordered_along_line_2_object() const { return Are_ordered_along_line_2(); }

typedef CGALi::Are_strictly_ordered_along_line         Are_strictly_ordered_along_line_2;
Are_strictly_ordered_along_line_2
are_strictly_ordered_along_line_2_object() const { return Are_strictly_ordered_along_line_2(); }

};

CGAL_END_NAMESPACE

// This macro is provided for convenience in defining the Kernel
// function objects inside a new representation class.
// See Cartesian_2.h and Cartesian.h

#define CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS_2(PR) \
 \
typedef PR::Less_xy_2    Less_xy_2; \
Less_xy_2 \
less_xy_2_object() const { return Less_xy_2(); } \
 \
typedef PR::Less_signed_distance_to_line_2     Less_signed_distance_to_line_2; \
Less_signed_distance_to_line_2 \
less_signed_distance_to_line_2_object() const { return Less_signed_distance_to_line_2(); } \
 \
typedef PR::Leftturn_2    Leftturn_2; \
Leftturn_2 \
leftturn_2_object() const { return Leftturn_2(); } \
 \
typedef PR::Left_of_line_2    Left_of_line_2; \
Left_of_line_2 \
left_of_line_2_object() const { return Left_of_line_2(); } \
 \
typedef PR::Collinear_2    Collinear_2; \
Collinear_2 \
collinear_2_object() const { return Collinear_2(); } \
 \
typedef PR::Orientation_2   Orientation_2; \
Orientation_2 \
orientation_2_object() const { return Orientation_2(); } \
 \
typedef PR::Side_of_oriented_circle_2    Side_of_oriented_circle_2; \
Side_of_oriented_circle_2 \
side_of_oriented_circle_2_object() const { return Side_of_oriented_circle_2(); } \
 \
typedef PR::Is_horizontal_2       Is_horizontal_2; \
Is_horizontal_2 \
is_horizontal_2_object() const { return Is_horizontal_2(); } \
 \
typedef PR::Is_vertical_2    Is_vertical_2; \
Is_vertical_2 \
is_vertical_2_object() const { return Is_vertical_2(); } \
 \
typedef PR::Is_degenerate_2     Is_degenerate_2; \
Is_degenerate_2 \
is_degenerate_2_object() const { return Is_degenerate_2(); } \
 \
typedef PR::Has_on_bounded_side_2     Has_on_bounded_side_2; \
Has_on_bounded_side_2 \
has_on_bounded_side_2_object() const { return Has_on_bounded_side_2(); } \
 \
typedef PR::Has_on_unbounded_side_2        Has_on_unbounded_side_2; \
Has_on_unbounded_side_2 \
has_on_unbounded_side_2_object() const { return Has_on_unbounded_side_2(); } \
 \
typedef PR::Has_on_boundary_2   Has_on_boundary_2; \
Has_on_boundary_2 \
has_on_boundary_2_object() const { return Has_on_boundary_2(); } \
 \
typedef PR::Has_on_positive_side_2      Has_on_positive_side_2; \
Has_on_positive_side_2 \
has_on_positive_side_2_object() const { return Has_on_positive_side_2(); } \
 \
typedef PR::Has_on_negative_side_2    Has_on_negative_side_2; \
Has_on_negative_side_2 \
has_on_negative_side_2_object() const { return Has_on_negative_side_2(); } \
 \
typedef PR::Oriented_side_2   Oriented_side_2; \
Oriented_side_2 \
oriented_side_2_object() const { return Oriented_side_2(); } \
 \
typedef PR::Compare_x_2      Compare_x_2; \
Compare_x_2 \
compare_x_2_object() const { return Compare_x_2(); } \
 \
typedef PR::Compare_y_2    Compare_y_2; \
Compare_y_2 \
compare_y_2_object() const { return Compare_y_2(); } \
 \
typedef PR::Compare_y_at_x_2    Compare_y_at_x_2; \
Compare_y_at_x_2 \
compare_y_at_x_2_object() const { return Compare_y_at_x_2(); } \
 \
typedef PR::Are_ordered_along_line_2    Are_ordered_along_line_2 ; \
Are_ordered_along_line_2 \
are_ordered_along_line_2_object() const { return Are_ordered_along_line_2(); } \
 \
typedef PR::Are_strictly_ordered_along_line_2     Are_strictly_ordered_along_line_2; \
Are_strictly_ordered_along_line_2 \
are_strictly_ordered_along_line_2_object() const { return Are_strictly_ordered_along_line_2(); }

#endif // CGAL_KERNEL_PREDICATE_OBJECTS_2_H
