// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_KERNEL_PREDICATE_OBJECTS_3_H
#define CGAL_KERNEL_PREDICATE_OBJECTS_3_H

#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Kernel_predicate_objects_3
{
public:
typedef typename R::FT                         FT;
typedef typename R::RT                         RT;
typedef typename R::Point_3                    Point_3;
typedef typename R::Vector_3                   Vector_3;
typedef typename R::Direction_3                Direction_3;
typedef typename R::Line_3                     Line_3;
typedef typename R::Plane_3                    Plane_3;
typedef typename R::Ray_3                      Ray_3;
typedef typename R::Segment_3                  Segment_3;
typedef typename R::Triangle_3                 Triangle_3;
typedef typename R::Tetrahedron_3              Tetrahedron_3;
typedef typename R::Aff_transformation_3       Aff_transformation_3;

typedef CGALi::Equal                           Equal_3;
Equal_3
equal_3_object() const 
{ return Equal_3(); }

typedef CGALi::Equal_x                         Equal_x_3;
Equal_x_3
equal_x_3_object() const 
{ return Equal_x_3(); }

typedef CGALi::Equal_y                         Equal_y_3;
Equal_y_3
equal_y_3_object() const 
{ return Equal_y_3(); }

typedef CGALi::Equal_z                         Equal_z_3;
Equal_z_3
equal_z_3_object() const 
{ return Equal_z_3(); }

typedef CGALi::Equal_xy                        Equal_xy_3;
Equal_xy_3
equal_xy_3_object() const 
{ return Equal_xy_3(); }

typedef CGALi::Equal_xyz                       Equal_xyz_3;
Equal_xyz_3
equal_xyz_3_object() const 
{ return Equal_xyz_3(); }

typedef CGALi::Less_x                          Less_x_3;
Less_x_3
less_x_3_object() const 
{ return Less_x_3(); }

typedef CGALi::Less_y                          Less_y_3;
Less_y_3
less_y_3_object() const 
{ return Less_y_3(); }

typedef CGALi::Less_z                          Less_z_3;
Less_z_3
less_z_3_object() const 
{ return Less_z_3(); }

typedef CGAL::p_Less_xy<Point_3>               Less_xy_3;
Less_xy_3
less_xy_3_object() const 
{ return Less_xy_3(); }

typedef CGALi::Less_xyz                        Less_xyz_3;
Less_xyz_3
less_xyz_3_object() const 
{ return Less_xyz_3(); }

typedef CGALi::Compare_x                       Compare_x_3;
Compare_x_3
compare_x_3_object() const 
{ return Compare_x_3(); }

typedef CGALi::Compare_y                       Compare_y_3;
Compare_y_3
compare_y_3_object() const 
{ return Compare_y_3(); }

typedef CGALi::Compare_z                       Compare_z_3;
Compare_z_3
compare_z_3_object() const 
{ return Compare_z_3(); }

typedef CGALi::Compare_xy                      Compare_xy_3;
Compare_xy_3
compare_xy_3_object() const 
{ return Compare_xyz_3(); }

typedef CGALi::Compare_xyz                     Compare_xyz_3;
Compare_xyz_3
compare_xyz_3_object() const 
{ return Compare_xyz_3(); }

typedef CGAL ::p_Less_dist_to_point<Point_3>   Less_distance_to_point_3;
Less_distance_to_point_3
less_distance_to_point_3_object() const 
{ return Less_distance_to_point_3(); }

typedef CGALi::Collinear                       Collinear_3;
Collinear_3
collinear_3_object() const 
{ return Collinear_3(); }

typedef CGALi::Coplanar                        Coplanar_3 ;
Coplanar_3
coplanar_3_object() const 
{ return Coplanar_3(); }

typedef CGAL ::p_Orientation<Point_3>          Orientation_3;
Orientation_3
orientation_3_object() const 
{ return Orientation_3(); }

typedef CGALi::Call_is_degenerate              Is_degenerate_3;
Is_degenerate_3
is_degenerate_3_object() const 
{ return Is_degenerate_3(); }

typedef CGALi::Call_has_on                     Has_on_3;
Has_on_3
has_on_3_object() const 
{ return Has_on_3(); }

typedef CGALi::Call_has_on_bounded_side        Has_on_bounded_side_3;
Has_on_bounded_side_3
has_on_bounded_side_3_object() const 
{ return Has_on_bounded_side_3(); }

typedef CGALi::Call_has_on_unbounded_side      Has_on_unbounded_side_3;
Has_on_unbounded_side_3
has_on_unbounded_side_3_object() const 
{ return Has_on_unbounded_side_3(); }

typedef CGALi::Call_has_on_boundary            Has_on_boundary_3;
Has_on_boundary_3
has_on_boundary_3_object() const 
{ return Has_on_boundary_3(); }

typedef CGALi::Call_has_on_positive_side       Has_on_positive_side_3;
Has_on_positive_side_3
has_on_positive_side_3_object() const 
{ return Has_on_positive_side_3(); }

typedef CGALi::Call_has_on_negative_side       Has_on_negative_side_3;
Has_on_negative_side_3
has_on_negative_side_3_object() const 
{ return Has_on_negative_side_3(); }

typedef CGALi::Call_oriented_side              Oriented_side_3;
Oriented_side_3
oriented_side_3_object() const 
{ return Oriented_side_3(); }

typedef CGALi::Are_ordered_along_line          Are_ordered_along_line_3 ;
Are_ordered_along_line_3
are_ordered_along_line_3_object() const 
{ return Are_ordered_along_line_3(); }

typedef CGALi::Are_strictly_ordered_along_line Are_strictly_ordered_along_line_3;
Are_strictly_ordered_along_line_3
are_strictly_ordered_along_line_3_object() const 
{ return Are_strictly_ordered_along_line_3(); }

typedef CGALi::Collinear_are_ordered_along_line Collinear_are_ordered_along_line_3;
Collinear_are_ordered_along_line_3
collinear_are_ordered_along_line_3_object() const 
{ return Collinear_are_ordered_along_line_3(); }

typedef CGALi::Collinear_are_strictly_ordered_along_line Collinear_are_strictly_ordered_along_line_3;
Collinear_are_strictly_ordered_along_line_3
collinear_are_strictly_ordered_along_line_3_object() const 
{ return Collinear_are_strictly_ordered_along_line_3(); }

typedef CGALi::Side_of_oriented_sphere         Side_of_oriented_sphere_3;
Side_of_oriented_sphere_3
side_of_oriented_sphere_3_object() const 
{ return Side_of_oriented_sphere_3(); }

typedef CGALi::Side_of_bounded_sphere          Side_of_bounded_sphere_3;
Side_of_bounded_sphere_3
side_of_bounded_sphere_3_object() const 
{ return Side_of_bounded_sphere_3(); }

};

CGAL_END_NAMESPACE

// This macro is provided for convenience in defining the Kernel
// function objects inside a new representation class.
// See Cartesian_3.h and Cartesian.h

#define CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS_3(PR) 

#endif // CGAL_KERNEL_PREDICATE_OBJECTS_3_H
