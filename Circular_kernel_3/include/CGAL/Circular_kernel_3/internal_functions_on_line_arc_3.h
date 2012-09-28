// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H

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

      Point_3 inters_p;
      Line_3 inters_l;

      Object o = SK().intersect_3_object()(l1.supporting_line(),
        l2.supporting_line());

      if(assign(inters_p, o)) {
        Circular_arc_point_3 p = inters_p;
        if(!SK().has_on_3_object()(l1,p,true)) return res;
        if(!SK().has_on_3_object()(l2,p,true)) return res;
        *res++ = make_object(std::make_pair(p,1u));
      } else if(assign(inters_l, o)) {

        if(SK().compare_xyz_3_object()(l1.lower_xyz_extremity(), 
                                       l2.lower_xyz_extremity()) < 0) {
	  int comparison = 
            SK().compare_xyz_3_object()(l2.lower_xyz_extremity(),
                                        l1.higher_xyz_extremity());
	  if(comparison < 0) {
	    if(SK().compare_xyz_3_object()(l1.higher_xyz_extremity(), 
                                           l2.higher_xyz_extremity()) <= 0) {
	      *res++ = make_object
	        (Line_arc_3(l1.supporting_line(),
                            l2.lower_xyz_extremity(),
                            l1.higher_xyz_extremity()));
	    } else {
	      *res++ = make_object(l2);
	    }
	  } else if (comparison == 0) {
	    *res++ = make_object(std::make_pair(l2.lower_xyz_extremity(),1u));
	  } 
        }
        else {
	  int comparison = 
            SK().compare_xyz_3_object()(l1.lower_xyz_extremity(),
                                        l2.higher_xyz_extremity());
	  if(comparison < 0){
	    if(SK().compare_xyz_3_object()(l1.higher_xyz_extremity(),
                                           l2.higher_xyz_extremity()) <= 0) {
	      *res++ = make_object(l1);
	    } else {
	      *res++ = make_object
	        (Line_arc_3(l1.supporting_line(), 
                            l1.lower_xyz_extremity(), 
                            l2.higher_xyz_extremity() ));
	    }
	  }
	  else if (comparison == 0){
	    *res++ = make_object(std::make_pair(l1.lower_xyz_extremity(),1u));
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

      Point_3 inters_p;
      Line_3 inters_l;

      Object o = SK().intersect_3_object()(l, la.supporting_line());

      if(assign(inters_l, o)) {
        *res++ = make_object(la);
      } else if(assign(inters_p, o)) {
        Circular_arc_point_3 p = inters_p;
        if(!SK().has_on_3_object()(la,p,true)) return res;
        *res++ = make_object(std::make_pair(p,1u));
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
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;
      solutions_container solutions;
      SK().intersect_3_object()(l.supporting_line(), c, 
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(l,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(l,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(l,sol2.first,true))
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
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;
      solutions_container solutions;
      SK().intersect_3_object()(l.supporting_line(), s,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(l,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(l,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(l,sol2.first,true))
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
        *res++ = make_object(l);
      }
      Point_3 sol;
      if(!assign(sol,SK().intersect_3_object()(p,l.supporting_line())))
        return res;
      if(!SK().has_on_3_object()(l,sol)) return res;
      Circular_arc_point_3 point = sol;
      *res++ = make_object(std::make_pair(point,1u));
      return res;
    }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H
