// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-17 $
// release_date  : $CGAL_Date: 2000/05/12 $
//
// file          : include/CGAL/Pm_with_intersections_misc.h
// package       : arr (1.27)
// maintainer    : Sigal Raab <raab@math.tau.ac.il>
// author(s)     : Eyal flato <flato@math.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_PM_WITH_INTERSECTIONS_MISC_H
#define CGAL_PM_WITH_INTERSECTIONS_MISC_H

#include <CGAL/Planar_map_2/Planar_map_misc.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

template <class I>
class Planar_map_with_intersections_traits_wrap : 
  public Planar_map_traits_wrap<I>
{
public:
  typedef  Planar_map_traits_wrap<I>    Base;
  typedef  typename Base::X_curve_2     X_curve_2;
  typedef  typename Base::Point_2       Point_2;

  typedef typename I::Has_left_category Has_left_category;
  
  Planar_map_with_intersections_traits_wrap() : Base() {}

  Planar_map_with_intersections_traits_wrap(const Base & tw) : Base(tw) {}

  void directed_curve_split(const X_curve_2 & cv, const Point_2 & first_pnt, 
			    const Point_2 & split_pnt,
			    X_curve_2 & split1, X_curve_2 & split2) 
  {
    CGAL_assertion(point_is_same(curve_source(cv), first_pnt) || 
	   point_is_same(curve_target(cv), first_pnt));
    CGAL_assertion(!point_is_same(curve_source(cv), split_pnt));
    CGAL_assertion(!point_is_same(curve_target(cv), split_pnt));

    curve_split(cv, split1, split2, split_pnt);
//      if (!point_is_same(curve_source(split1), first_pnt)) //start flip
//        {
//  	X_curve_2 c = split2;
//  	split2 = curve_flip(split1);
//  	split1 = curve_flip(c);
//        }// end flip
  }

// expand (split) the segment[p1, p2] from the curve cv
//   X_curve_2 curve_portion(const X_curve_2 & cv,
//    const Point_2 & p1, const Point_2 & p2)
//   {
//     CGAL_assertion(!point_is_same(p1, p2));
//     // direct cv as p1-->p2
//     X_curve_2 acv = cv;
//     bool p1_p2_right;//start flip
//     bool curve_right;
//     Point_2 ap1=p1, ap2=p2;

//      p1_p2_right = point_is_left_low(p1, p2);
//      curve_right = point_is_left_low(curve_source(cv), curve_target(cv));
//      if (curve_right != p1_p2_right) {
//             acv = curve_flip(cv);//end flip
// //        ap1 = p2;
// //        ap2 = p1;
//      }

//     // split twice to find the portion [ap1, ap2]
//     X_curve_2 split1, split2, split3, part_cv;
//     CGAL_assertion(!point_is_same(curve_target(acv), ap1));
//     if (point_is_same(curve_source(acv), ap1))
//       split2 = acv;
//     else
//       directed_curve_split(acv, curve_source(acv), ap1, split1, split2);
//     CGAL_assertion(!point_is_same(curve_source(split2), ap2));
//     if (point_is_same(curve_target(split2), ap2))
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
  bool nearest_intersection_to_left(const X_curve_2 & cv1,
                                    const X_curve_2 & cv2,
                                    const Point_2 & pt,
                                    Point_2 & p1, Point_2 & p2) const 
  {
    return nearest_intersection_to_left_imp(cv1, cv2, pt, p1, p2,
                                            Has_left_category());
  }

    
  bool nearest_intersection_to_left_imp(const X_curve_2 & cv1,
                                        const X_curve_2 & cv2,
                                        const Point_2 & pt,
                                        Point_2 & p1, Point_2 & p2,
                                        Tag_true) const
  { return Base::nearest_intersection_to_left(cv1, cv2, pt, p1, p2); }
    
  bool nearest_intersection_to_left_imp(const X_curve_2 & cv1,
                                        const X_curve_2 & cv2,
                                        const Point_2 & pt,
                                        Point_2 & p1, Point_2 & p2,
                                        Tag_false) const 
  {
    Point_2 rpt = point_reflect_in_x_and_y( pt);
    X_curve_2 rcv1 = curve_reflect_in_x_and_y( cv1);
    X_curve_2 rcv2 = curve_reflect_in_x_and_y( cv2);

    Point_2 rp1, rp2;
    bool result = nearest_intersection_to_right(rcv1, rcv2, rpt, rp1, rp2);

    p1 = point_reflect_in_x_and_y( rp1);
    p2 = point_reflect_in_x_and_y( rp2);

    return result;
  }

  // maps the curves to their mirror images over the y coordinate
  // and calls do_intersect_to_right (see there).
  bool do_intersect_to_left(const X_curve_2 & ca, const X_curve_2 & cb,
			    const Point_2 & pt) const
  { return do_intersect_to_left_imp(ca, cb, pt, Has_left_category()); }

  bool do_intersect_to_left_imp(const X_curve_2 & ca, const X_curve_2 & cb,
                                const Point_2 & pt, Tag_true) const
  { return Base::do_intersect_to_left(ca, cb, pt); }
    
  bool do_intersect_to_left_imp(const X_curve_2 & ca, const X_curve_2 & cb,
                                const Point_2 & pt, Tag_false) const
  {
      Point_2 rpt = point_reflect_in_x_and_y( pt);
      X_curve_2 rca = curve_reflect_in_x_and_y( ca);
      X_curve_2 rcb = curve_reflect_in_x_and_y( cb);

      return do_intersect_to_right(rca, rcb, rpt);
  }

};

CGAL_END_NAMESPACE

#endif
