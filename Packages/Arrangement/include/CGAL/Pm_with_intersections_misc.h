// Copyright (c) 2000  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Eyal flato <flato@math.tau.ac.il>
#ifndef CGAL_PM_WITH_INTERSECTIONS_MISC_H
#define CGAL_PM_WITH_INTERSECTIONS_MISC_H

#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

template <class I>
class Planar_map_with_intersections_traits_wrap : 
  public Pm_traits_wrap_2<I>
{
public:
  typedef Pm_traits_wrap_2<I>                   Base;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Point_2                Point_2;

  typedef typename I::Has_left_category         Has_left_category;
  
  Planar_map_with_intersections_traits_wrap() : Base() {}

  Planar_map_with_intersections_traits_wrap(const Base & tw) : Base(tw) {}

  void directed_curve_split(const X_monotone_curve_2 & cv,
                            const Point_2 & first_pnt, 
			    const Point_2 & split_pnt,
			    X_monotone_curve_2 & split1,
                            X_monotone_curve_2 & split2) 
  {
    CGAL_assertion(point_equal(curve_source(cv), first_pnt) || 
	   point_equal(curve_target(cv), first_pnt));
    CGAL_assertion(!point_equal(curve_source(cv), split_pnt));
    CGAL_assertion(!point_equal(curve_target(cv), split_pnt));

    curve_split(cv, split1, split2, split_pnt);
//      if (!point_equal(curve_source(split1), first_pnt)) //start flip
//        {
//  	X_monotone_curve_2 c = split2;
//  	split2 = curve_opposite(split1);
//  	split1 = curve_opposite(c);
//        }// end flip
  }

// expand (split) the segment[p1, p2] from the curve cv
//   X_monotone_curve_2 curve_portion(const X_monotone_curve_2 & cv,
//    const Point_2 & p1, const Point_2 & p2)
//   {
//     CGAL_assertion(!point_equal(p1, p2));
//     // direct cv as p1-->p2
//     X_monotone_curve_2 acv = cv;
//     bool p1_p2_right;//start flip
//     bool curve_right;
//     Point_2 ap1=p1, ap2=p2;

//      p1_p2_right = point_is_left_low(p1, p2);
//      curve_right = point_is_left_low(curve_source(cv), curve_target(cv));
//      if (curve_right != p1_p2_right) {
//             acv = curve_opposite(cv);//end flip
// //        ap1 = p2;
// //        ap2 = p1;
//      }

//     // split twice to find the portion [ap1, ap2]
//     X_monotone_curve_2 split1, split2, split3, part_cv;
//     CGAL_assertion(!point_equal(curve_target(acv), ap1));
//     if (point_equal(curve_source(acv), ap1))
//       split2 = acv;
//     else
//       directed_curve_split(acv, curve_source(acv), ap1, split1, split2);
//     CGAL_assertion(!point_equal(curve_source(split2), ap2));
//     if (point_equal(curve_target(split2), ap2))
//       part_cv = split2;
//     else
//       directed_curve_split(split2, curve_source(split2), ap2, part_cv, 
  //split3);
//     return part_cv;
//   }

  void points_swap(Point_2 & p1, Point_2 & p2)
  {
    Point_2 p = p1;
    p1 = p2;
    p2 = p;
  }

  // maps the curves to their mirror images over the y coordinate
  // and calls nearest_intersect_to_right (see there).
  bool nearest_intersection_to_left(const X_monotone_curve_2 & cv1,
                                    const X_monotone_curve_2 & cv2,
                                    const Point_2 & pt,
                                    Point_2 & p1, Point_2 & p2) const 
  {
    return nearest_intersection_to_left_imp(cv1, cv2, pt, p1, p2,
                                            Has_left_category());
  }

    
  bool nearest_intersection_to_left_imp(const X_monotone_curve_2 & cv1,
                                        const X_monotone_curve_2 & cv2,
                                        const Point_2 & pt,
                                        Point_2 & p1, Point_2 & p2,
                                        Tag_true) const
  { return Base::nearest_intersection_to_left(cv1, cv2, pt, p1, p2); }
    
  bool nearest_intersection_to_left_imp(const X_monotone_curve_2 & cv1,
                                        const X_monotone_curve_2 & cv2,
                                        const Point_2 & pt,
                                        Point_2 & p1, Point_2 & p2,
                                        Tag_false) const 
  {
    Point_2 rpt = point_reflect_in_x_and_y( pt);
    X_monotone_curve_2 rcv1 = curve_reflect_in_x_and_y( cv1);
    X_monotone_curve_2 rcv2 = curve_reflect_in_x_and_y( cv2);

    Point_2 rp1, rp2;
    bool result = nearest_intersection_to_right(rcv1, rcv2, rpt, rp1, rp2);

    p1 = point_reflect_in_x_and_y( rp1);
    p2 = point_reflect_in_x_and_y( rp2);

    return result;
  }

};

CGAL_END_NAMESPACE

#endif
