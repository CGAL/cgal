// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Kernel_checker.h
// revision      : $Revision$
// revision_date : $Date$
// package       : ???
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_KERNEL_CHECKER_H
#define CGAL_KERNEL_CHECKER_H

// This file contains the definition of a kernel traits checker.
//
// TODO:
// - have a look at the PM_checker from Tel-Aviv, and Stefan's NT checker.
// - At the moment, only predicates are checked.  To handle constructions as
//   well, the best approach is probably to have objects be pairs.
//   So the template parameter will be a comparator, not a converter.

#include <CGAL/basic.h>
#include <utility>

CGAL_BEGIN_NAMESPACE

// Class used by Kernel_checker.
// This works only for predicates.
// Something else must be done to handle the constructions.
// Maybe the solution is to store pairs of objects ?
template <class O1, class O2, class Conv>
class Predicate_checker
{
    O1 o1;
    O2 o2;
    Conv c;

public:

    Predicate_checker() {}

    typedef typename O1::result_type result_type;

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2) const
    {
	typename O1::result_type res1 = o1(a1, a2);
	typename O2::result_type res2 = o2(c(a1), c(a2));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << std::endl;
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3);
	typename O2::result_type res2 = o2(c(a1), c(a2), c(a3));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << std::endl;
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3, a4);
	typename O2::result_type res2 = o2(c(a1), c(a2), c(a3), c(a4));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << std::endl;
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4, class A5>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3, a4, a5);
	typename O2::result_type res2 = o2(c(a1), c(a2), c(a3), c(a4), c(a5));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << ", " << a5 << std::endl;
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    // Same thing with more arguments...
};

// For now, we inherit all geometric objects and constructions from K1, and
// just overload the predicates.
template <class K1, class K2, class Conv>
class Kernel_checker
  : public K1
{
    typedef K1     Kernel1;
    typedef K2     Kernel2;
    // typedef Comp   Comparator;
    typedef Conv   c;

    // typedef std::pair<K1::Point_2, K2::Point_2>  Point_2;
    // ...  Same thing for all objects.

#define CGAL_check_pred(X, Y) \
    typedef Predicate_checker<typename K1::X, typename K2::X, Conv> X; \
    X Y() const { return X(); }

public:
    CGAL_check_pred(Equal_2, equal_2_object)
    CGAL_check_pred(Equal_x_2, equal_x_2_object)
    CGAL_check_pred(Equal_y_2, equal_y_2_object)
    CGAL_check_pred(Equal_xy_2, equal_xy_2_object)
    CGAL_check_pred(Less_x_2, less_x_2_object)
    CGAL_check_pred(Less_y_2, less_y_2_object)
    CGAL_check_pred(Less_xy_2, less_xy_2_object)
    CGAL_check_pred(Less_yx_2, less_yx_2_object)
    CGAL_check_pred(Compare_x_2, compare_x_2_object)
    CGAL_check_pred(Compare_y_2, compare_y_2_object)
    CGAL_check_pred(Compare_xy_2, compare_xy_2_object)
    CGAL_check_pred(Compare_y_at_x_2, compare_y_at_x_2_object)
    CGAL_check_pred(Compare_distance_2, compare_distance_2_object)
    CGAL_check_pred(Counterclockwise_in_between_2,
	counterclockwise_in_between_2_object)
    CGAL_check_pred(Leftturn_2, leftturn_2_object)
    CGAL_check_pred(Collinear_2, collinear_2_object)
    CGAL_check_pred(Orientation_2, orientation_2_object)
    CGAL_check_pred(Side_of_oriented_circle_2,
	    side_of_oriented_circle_2_object)
    CGAL_check_pred(Side_of_bounded_circle_2, side_of_bounded_circle_2_object)
    CGAL_check_pred(Is_horizontal_2, is_horizontal_2_object)
    CGAL_check_pred(Is_vertical_2, is_vertical_2_object)
    CGAL_check_pred(Is_degenerate_2, is_degenerate_2_object)
    CGAL_check_pred(Has_on_2, has_on_2_object)
    CGAL_check_pred(Collinear_has_on_2, collinear_has_on_2_object)
    CGAL_check_pred(Has_on_bounded_side_2, has_on_bounded_side_2_object)
    CGAL_check_pred(Has_on_unbounded_side_2, has_on_unbounded_side_2_object)
    CGAL_check_pred(Has_on_boundary_2, has_on_boundary_2_object)
    CGAL_check_pred(Has_on_positive_side_2, has_on_positive_side_2_object)
    CGAL_check_pred(Has_on_negative_side_2, has_on_negative_side_2_object)
    CGAL_check_pred(Oriented_side_2, oriented_side_2_object)
    CGAL_check_pred(Are_ordered_along_line_2, are_ordered_along_line_2_object)
    CGAL_check_pred(Are_strictly_ordered_along_line_2,
	are_strictly_ordered_along_line_2_object)
    CGAL_check_pred(Collinear_are_ordered_along_line_2,
	collinear_are_ordered_along_line_2_object)
    CGAL_check_pred(Collinear_are_strictly_ordered_along_line_2,
	collinear_are_strictly_ordered_along_line_2_object)

    CGAL_check_pred(Equal_3, equal_3_object)
    CGAL_check_pred(Equal_x_3, equal_x_3_object)
    CGAL_check_pred(Equal_y_3, equal_y_3_object)
    CGAL_check_pred(Equal_z_3, equal_z_3_object)
    CGAL_check_pred(Equal_xy_3, equal_xy_3_object)
    CGAL_check_pred(Equal_xyz_3, equal_xyz_3_object)
    CGAL_check_pred(Less_x_3, less_x_3_object)
    CGAL_check_pred(Less_y_3, less_y_3_object)
    CGAL_check_pred(Less_z_3, less_z_3_object)
    CGAL_check_pred(Less_xy_3, less_xy_3_object)
    CGAL_check_pred(Less_xyz_3, less_xyz_3_object)
    CGAL_check_pred(Compare_x_3, compare_x_3_object)
    CGAL_check_pred(Compare_y_3, compare_y_3_object)
    CGAL_check_pred(Compare_z_3, compare_z_3_object)
    CGAL_check_pred(Compare_xy_3, compare_xy_3_object)
    CGAL_check_pred(Compare_xyz_3, compare_xyz_3_object)
    CGAL_check_pred(Compare_distance_3, compare_distance_3_object)
    CGAL_check_pred(Collinear_3, collinear_3_object)
    CGAL_check_pred(Coplanar_3, coplanar_3_object)
    CGAL_check_pred(Coplanar_orientation_3, coplanar_orientation_3_object)
    CGAL_check_pred(Coplanar_side_of_bounded_circle_3,
	    coplanar_side_of_bounded_circle_3_object)
    CGAL_check_pred(Orientation_3, orientation_3_object)
    CGAL_check_pred(Is_degenerate_3, is_degenerate_3_object)
    CGAL_check_pred(Has_on_3, has_on_3_object)
    CGAL_check_pred(Has_on_bounded_side_3, has_on_bounded_side_3_object)
    CGAL_check_pred(Has_on_unbounded_side_3, has_on_unbounded_side_3_object)
    CGAL_check_pred(Has_on_boundary_3, has_on_boundary_3_object)
    CGAL_check_pred(Has_on_positive_side_3, has_on_positive_side_3_object)
    CGAL_check_pred(Has_on_negative_side_3, has_on_negative_side_3_object)
    CGAL_check_pred(Oriented_side_3, oriented_side_3_object)
    CGAL_check_pred(Are_ordered_along_line_3, are_ordered_along_line_3_object)
    CGAL_check_pred(Are_strictly_ordered_along_line_3,
	    are_strictly_ordered_along_line_3_object)
    CGAL_check_pred(Collinear_are_ordered_along_line_3,
	    collinear_are_ordered_along_line_3_object)
    CGAL_check_pred(Collinear_are_strictly_ordered_along_line_3,
	    collinear_are_strictly_ordered_along_line_3)
    CGAL_check_pred(Side_of_oriented_sphere_3,
	    side_of_oriented_sphere_3_object)
    CGAL_check_pred(Side_of_bounded_sphere_3, side_of_bounded_sphere_3_object)
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_CHECKER_H
