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

#ifndef CGAL_PLANAR_MAP_MISC_H
#include <CGAL/Planar_map_2/Planar_map_misc.h>
#endif

CGAL_BEGIN_NAMESPACE

template<class Planar_map_>
class Pm_change_notification
{
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;

  virtual void add_edge(
			const typename Traits::X_curve& cv, 
			typename Planar_map::Halfedge_handle e, 
			bool original_direction, bool overlap=false)
  {
  }

  virtual void split_edge(
			  typename Planar_map::Halfedge_handle orig_edge, 
			  typename Planar_map::Halfedge_handle new_edge,
			  const typename Traits::X_curve& c1,
			  const typename Traits::X_curve& c2)
  {
  }

  virtual void split_face(
			  typename Planar_map::Face_handle orig_face, 
			  typename Planar_map::Face_handle new_face)
  {
  }
	
  virtual void add_hole(
			typename Planar_map::Face_handle in_face, 
			typename Planar_map::Halfedge_handle new_hole)
  {
  }

  virtual const typename Traits::X_curve&
  edge_support_curve(typename Planar_map::Halfedge_handle edge)
  {
    return edge->curve();
  }

  virtual bool have_support_curve()
  {
    return false;
  }

};


template <class I>
class Planar_map_with_intersections_traits_wrap : 
  public Planar_map_traits_wrap<I>
{
public:
  typedef  Planar_map_traits_wrap<I> Base;
  typedef  typename Base::X_curve    X_curve;
  typedef  typename Base::Point      Point;
  
  Planar_map_with_intersections_traits_wrap() : Base()
  {
  }

  Planar_map_with_intersections_traits_wrap(const Base& tw) : Base(tw)
  {
  }

  void directed_curve_split(const X_curve &cv, const Point &first_pnt, 
			    const Point &split_pnt,
			    X_curve &split1, X_curve &split2) 
  {
    CGAL_assertion(point_is_same(curve_source(cv), first_pnt) || 
	   point_is_same(curve_target(cv), first_pnt));
    CGAL_assertion(!point_is_same(curve_source(cv), split_pnt));
    CGAL_assertion(!point_is_same(curve_target(cv), split_pnt));

    curve_split(cv, split1, split2, split_pnt);
//      if (!point_is_same(curve_source(split1), first_pnt)) //start flip
//        {
//  	X_curve c = split2;
//  	split2 = curve_flip(split1);
//  	split1 = curve_flip(c);
//        }// end flip
  }

  // expand (split) the segment[p1, p2] from the curve cv
//   X_curve curve_portion(const X_curve &cv, const Point &p1, const Point &p2)
//   {
//     CGAL_assertion(!point_is_same(p1, p2));
//     // direct cv as p1-->p2
//     X_curve acv = cv;
//     bool p1_p2_right;//start flip
//     bool curve_right;
//     Point ap1=p1, ap2=p2;

//      p1_p2_right = point_is_left_low(p1, p2);
//      curve_right = point_is_left_low(curve_source(cv), curve_target(cv));
//      if (curve_right != p1_p2_right) {
//             acv = curve_flip(cv);//end flip
// //        ap1 = p2;
// //        ap2 = p1;
//      }

//     // split twice to find the portion [ap1, ap2]
//     X_curve split1, split2, split3, part_cv;
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

  void points_swap(Point &p1, Point &p2)
  {
    Point p = p1;
    p1 = p2;
    p2 = p;
  }

#ifndef CGAL_PMWX_TRAITS_HAVE_INTERSECT_TO_LEFT

  // maps the curves to their mirror images over the y coordinate
  // and calls nearest_intersect_to_right (see there).
  bool nearest_intersection_to_left(const X_curve& cv1,
                                    const X_curve& cv2,
                                    const Point& pt,
                                    Point& p1,
                                    Point& p2) const 
    {
      Point rpt = point_reflect_in_x_and_y( pt);
      X_curve rcv1 = curve_reflect_in_x_and_y( cv1);
      X_curve rcv2 = curve_reflect_in_x_and_y( cv2);

      Point rp1, rp2;
      bool result = nearest_intersection_to_right(rcv1, rcv2, rpt, rp1, rp2);

      p1 = point_reflect_in_x_and_y( rp1);
      p2 = point_reflect_in_x_and_y( rp2);

      return result;
    }


  // maps the curves to their mirror images over the y coordinate
  // and calls do_intersect_to_right (see there).
  bool do_intersect_to_left(const X_curve& ca, const X_curve& cb,
			    const Point& pt) const
  {
      Point rpt = point_reflect_in_x_and_y( pt);
      X_curve rca = curve_reflect_in_x_and_y( ca);
      X_curve rcb = curve_reflect_in_x_and_y( cb);

      return do_intersect_to_right(rca, rcb, rpt);
  }

#endif //CGAL_PMWX_TRAITS_HAVE_INTERSECT_TO_LEFT

};

CGAL_END_NAMESPACE

#endif
