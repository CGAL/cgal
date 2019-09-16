// Copyright (c) 2000-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Herve Bronnimann, Sylvain Pion, Susan Hert

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred, CGAL_Kernel_cons and CGAL_Kernel_obj.
// And they are #undefed at the end of this file.

#ifndef CGAL_Kernel_pred
#  define CGAL_Kernel_pred(X, Y)
#endif

#ifndef CGAL_Kernel_cons
#  define CGAL_Kernel_cons(X, Y)
#endif

#ifndef CGAL_Kernel_obj
#  define CGAL_Kernel_obj(X)
#endif

CGAL_Kernel_obj(Point_d)
CGAL_Kernel_obj(Vector_d)
CGAL_Kernel_obj(Direction_d)
CGAL_Kernel_obj(Hyperplane_d)
CGAL_Kernel_obj(Sphere_d)
CGAL_Kernel_obj(Iso_box_d)
CGAL_Kernel_obj(Segment_d)
CGAL_Kernel_obj(Ray_d)
CGAL_Kernel_obj(Line_d)


CGAL_Kernel_pred(Affinely_independent_d,
                 affinely_independent_d_object)
CGAL_Kernel_pred(Affine_rank_d,
                 affine_rank_d_object)
CGAL_Kernel_pred(Compare_lexicographically_d,
                 compare_lexicographically_d_object)
CGAL_Kernel_pred(Contained_in_affine_hull_d,
                 contained_in_affine_hull_d_object)
CGAL_Kernel_pred(Contained_in_linear_hull_d,
                 contained_in_linear_hull_d_object)
CGAL_Kernel_pred(Contained_in_simplex_d,
                 contained_in_simplex_d_object)
// TODO: create a Do_intersect_d functor
//CGAL_Kernel_pred(Do_intersect_d,
//                 do_intersect_d_object)
CGAL_Kernel_pred(Less_lexicographically_d,
                 less_lexicographically_d_object)
CGAL_Kernel_pred(Less_or_equal_lexicographically_d,
                 less_or_equal_lexicographically_d_object)
CGAL_Kernel_pred(Linearly_independent_d,
                 linearly_independent_d_object)
CGAL_Kernel_pred(Linear_rank_d,
                 linear_rank_d_object)
CGAL_Kernel_pred(Orientation_d,
                 orientation_d_object)
CGAL_Kernel_pred(Coaffine_orientation_d,
                 coaffine_orientation_d_object)
CGAL_Kernel_pred(Side_of_bounded_sphere_d,
                 side_of_bounded_sphere_d_object)
CGAL_Kernel_pred(Side_of_oriented_sphere_d,
                 side_of_oriented_sphere_d_object)
CGAL_Kernel_pred(Side_of_oriented_subsphere_d,
                 side_of_oriented_subsphere_d_object)
CGAL_Kernel_pred(Oriented_side_d,
                 oriented_side_d_object)


CGAL_Kernel_cons(Linear_base_d,
                 linear_base_d_object)
CGAL_Kernel_cons(Center_of_sphere_d,
                 center_of_sphere_d_object)
CGAL_Kernel_cons(Intersection_d_,
                 intersection_d_object)
CGAL_Kernel_cons(Lift_to_paraboloid_d,
                 lift_to_paraboloid_d_object)
CGAL_Kernel_cons(Midpoint_d,
                 midpoint_d_object)
CGAL_Kernel_cons(Project_along_d_axis_d,
                 project_along_d_axis_d_object)
CGAL_Kernel_cons(Squared_distance_d,
                 squared_distance_d_object)

#undef CGAL_Kernel_pred
#undef CGAL_Kernel_cons
#undef CGAL_Kernel_obj
