// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Circular_kernel_3/Intersection_traits.h>

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
    equal( const typename SK::Line_arc_3 &l1,
           const typename SK::Line_arc_3 &l2)
    {
      return l1.rep() == l2.rep();
    }

    template <class SK>
    inline
    bool
    do_overlap(const typename SK::Line_arc_3 &l1,
               const typename SK::Line_arc_3 &l2,
               const bool known_equal_supporting_line = false)
    {
      if(!known_equal_supporting_line) {
        if (!non_oriented_equal<SK>(l1.supporting_line(),
                                    l2.supporting_line()))
          return false;
      }

      return SK().compare_xyz_3_object()(l1.higher_xyz_extremity(),
                             l2.lower_xyz_extremity()) >= 0
          && SK().compare_xyz_3_object()(l1.lower_xyz_extremity(),
                             l2.higher_xyz_extremity()) <= 0;
    }

    template < class SK >
    void
    split(const typename SK::Line_arc_3 &l,
          const typename SK::Circular_arc_point_3 &p,
          typename SK::Line_arc_3 &l1,
          typename SK::Line_arc_3 &l2)
    {
      typedef typename SK::Line_arc_3  Line_arc_3;
      // The point must be on the line arc
      CGAL_kernel_precondition(SK().has_on_3_object()(l, p));
      // It doesn't make sense to split an arc on an extremity
      CGAL_kernel_precondition(l.source() != p);
      CGAL_kernel_precondition(l.target() != p);
      if(SK().compare_xyz_3_object()(l.source(),p) == SMALLER) {
        l1 = Line_arc_3(l.supporting_line(),l.source(),p);
        l2 = Line_arc_3(l.supporting_line(),p,l.target());
      } else {
        l1 = Line_arc_3(l.supporting_line(),p,l.target());
        l2 = Line_arc_3(l.supporting_line(),l.source(),p);
      }
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Line_arc_3 & l1,
                const typename SK::Line_arc_3 & l2,
                OutputIterator res)
    {
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK::Line_3 Line_3;
      typedef typename SK::Line_arc_3 Line_arc_3;
      typedef typename SK3_Intersection_traits<SK, Line_arc_3, Line_arc_3>::type result_type;

      typename Intersection_traits<SK, Line_3, Line_3>::result_type o =
        SK().intersect_3_object()(l1.supporting_line(),
        l2.supporting_line());

      if(!o)
        return res;

      if(const Point_3* inters_p = CGAL::Intersections::internal::intersect_get<Point_3>(o)) {
          Circular_arc_point_3 p = *inters_p;
        if(!SK().has_on_3_object()(l1,p,true)) return res;
        if(!SK().has_on_3_object()(l2,p,true)) return res;
          *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(p,1u));
      } else if( CGAL::Intersections::internal::intersect_get<Line_3>(o) ) {
        if(SK().compare_xyz_3_object()(l1.lower_xyz_extremity(),
                                       l2.lower_xyz_extremity()) < 0) {
          int comparison =
            SK().compare_xyz_3_object()(l2.lower_xyz_extremity(),
                                        l1.higher_xyz_extremity());
          if(comparison < 0) {
            if(SK().compare_xyz_3_object()(l1.higher_xyz_extremity(),
                                           l2.higher_xyz_extremity()) <= 0) {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>
                (Line_arc_3(l1.supporting_line(),
                            l2.lower_xyz_extremity(),
                            l1.higher_xyz_extremity()));
            } else {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>(l2);
            }
          } else if (comparison == 0) {
            *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(l2.lower_xyz_extremity(),1u));
          }
        }
        else {
          int comparison =
            SK().compare_xyz_3_object()(l1.lower_xyz_extremity(),
                                        l2.higher_xyz_extremity());
          if(comparison < 0){
            if(SK().compare_xyz_3_object()(l1.higher_xyz_extremity(),
                                           l2.higher_xyz_extremity()) <= 0) {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>(l1);
            } else {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>
                (Line_arc_3(l1.supporting_line(),
                            l1.lower_xyz_extremity(),
                            l2.higher_xyz_extremity() ));
            }
          }
          else if (comparison == 0){
            *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(l1.lower_xyz_extremity(),1u));
          }
        }
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Line_3 & l,
                const typename SK::Line_arc_3 & la,
                OutputIterator res)
    {
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK::Line_3 Line_3;
      typedef typename SK::Line_arc_3 Line_arc_3;
      typedef typename SK3_Intersection_traits<SK, Line_3, Line_arc_3>::type result_type;

      typename Intersection_traits<SK, Line_3, Line_3>::result_type o =
        SK().intersect_3_object()(l, la.supporting_line());

      if(!o)
        return res;

      if(const Line_3* inters_l = CGAL::Intersections::internal::intersect_get<Line_3>(o)) {
        *res++ = CGAL::internal::sk3_intersection_return<result_type>(la);
      } else if(const Point_3* inters_p = CGAL::Intersections::internal::intersect_get<Point_3>(o)) {
        Circular_arc_point_3 p = *inters_p;
        if(!SK().has_on_3_object()(la,p,true)) return res;
        *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(p,1u));
      }

      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c,
                const typename SK::Line_arc_3 & l,
                OutputIterator res)
    {
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<
        typename SK3_Intersection_traits<SK, typename SK::Line_3, typename SK::Circle_3
                                  >::type> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      solutions_container solutions;
      SK().intersect_3_object()(l.supporting_line(), c,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
        const Solution* sol = CGAL::Intersections::internal::intersect_get<Solution>(solutions[0]);
        if(SK().has_on_3_object()(l,(*sol).first,true))
           *res++ = solutions[0];
      } else {
         const Solution* sol1 = CGAL::Intersections::internal::intersect_get<Solution>(solutions[0]);
         const Solution* sol2 = CGAL::Intersections::internal::intersect_get<Solution>(solutions[1]);

         if(SK().has_on_3_object()(l,(*sol1).first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(l,(*sol2).first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Sphere_3 & s,
                const typename SK::Line_arc_3 & l,
               OutputIterator res)
    {
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<
        typename SK3_Intersection_traits<SK, typename SK::Line_3, typename SK::Sphere_3>::type
      > solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;
      solutions_container solutions;
      SK().intersect_3_object()(l.supporting_line(), s,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
        const Solution* sol = CGAL::Intersections::internal::intersect_get<Solution>(solutions[0]);
        if(SK().has_on_3_object()(l,(*sol).first,true))
           *res++ = solutions[0];
      } else {
         const Solution* sol1 = CGAL::Intersections::internal::intersect_get<Solution>(solutions[0]);
         const Solution* sol2 = CGAL::Intersections::internal::intersect_get<Solution>(solutions[1]);
         if(SK().has_on_3_object()(l,(*sol1).first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(l,(*sol2).first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Plane_3 & p,
                const typename SK::Line_arc_3 & l,
                OutputIterator res)
    {
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      if(SK().has_on_3_object()(p,l.supporting_line())) {
        *res++ = result_type(l);
      }
      const Point_3* sol;
      typename Intersection_traits<SK, typename SK::Plane_3, typename SK::Line_3>
        ::result_type o = SK().intersect_3_object()(p,l.supporting_line());

      if(!o)
        return res;
      if((sol = CGAL::Intersections::internal::intersect_get<Point_3>(o))) {
        if(!SK().has_on_3_object()(l,*sol)) return res;
      Circular_arc_point_3 point = sol;
        *res++ = result_type(std::make_pair(point,1u));
      }
      return res;
    }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H
