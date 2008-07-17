// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H

#include <CGAL/Circular_kernel_3/internal_function_has_on_spherical_kernel.h>

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
    equal( const typename SK::Circular_arc_3 &c1,
           const typename SK::Circular_arc_3 &c2)
    {
      return c1.rep() == c2.rep();
    }

    template <class SK>
    inline
    bool
    do_overlap(const typename SK::Circular_arc_3 &c1,
               const typename SK::Circular_arc_3 &c2,
               const bool known_equal_supporting_circle = false)
    { 
      if(!known_equal_supporting_circle) {
        if(!non_oriented_equal<SK>(c1.supporting_circle(), 
                                   c2.supporting_circle()))
          return false;
      }
      if(c1.rep().is_full()) return true;
      if(c2.rep().is_full()) return true;
      if((SK().has_on_3_object()(c1,c2.target(),true)) || 
         (SK().has_on_3_object()(c1,c2.source(),true))) return true;
      return SK().has_on_3_object()(c2,c1.source(),true);
    }

    template < class SK >
    void
    split(const typename SK::Circular_arc_3 &c,
	  const typename SK::Circular_arc_point_3 &p,
	  typename SK::Circular_arc_3 &c1,
	  typename SK::Circular_arc_3 &c2)
    {
      // The point must be on the circular arc 
      CGAL_kernel_precondition(SK().has_on_3_object()(c, p));
      typedef typename SK::Circular_arc_3  Circular_arc_3;
      // It doesn't make sense to split an arc on an extremity
      CGAL_kernel_precondition(c.source() != p);
      CGAL_kernel_precondition(c.target() != p);
      const Circular_arc_3 &rc1 = 
        Circular_arc_3(c.supporting_circle(), c.source(), p);
      const Circular_arc_3 &rc2 = 
        Circular_arc_3(c.supporting_circle(), p, c.target());
      if ( SK().compare_xyz_3_object()(rc1.source(), rc2.source()) != 
           SMALLER) {
        c1 = rc2; c2 = rc1;
      } else { c1 = rc1; c2 = rc2; }
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Line_3 & l, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      solutions_container solutions;

      SK().intersect_3_object()(l, ca.supporting_circle(),
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(ca,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      if(non_oriented_equal<SK>(c, ca.supporting_circle())) {
        *res++ = make_object(ca);
      }

      solutions_container solutions;

      SK().intersect_3_object()(ca.supporting_circle(), c, 
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(ca,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Sphere_3 & s,
                const typename SK::Circular_arc_3 & c,
	       OutputIterator res)
    {
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      if(SK().has_on_3_object()(s, c.supporting_circle())) {
        *res++ = make_object(c);
      }

      solutions_container solutions;

      SK().intersect_3_object()(c.supporting_circle(), s,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(c,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(c,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(c,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Plane_3 & p, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    {
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;
      if(SK().has_on_3_object()(p,ca.supporting_circle())) {
        *res++ = make_object(ca);
      }
      solutions_container solutions;

      SK().intersect_3_object()(ca.supporting_circle(), p,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(ca,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Line_arc_3 & la, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      solutions_container solutions;

      SK().intersect_3_object()(la.supporting_line(), ca.supporting_circle(),
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         Solution sol;
         assign(sol, solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true) &&
            SK().has_on_3_object()(la,sol.first,true))
           *res++ = solutions[0];
      } else {
         Solution sol1, sol2;
         assign(sol1, solutions[0]);
         assign(sol2, solutions[1]);
         if(SK().has_on_3_object()(ca,sol1.first,true) &&
            SK().has_on_3_object()(la,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true) &&
            SK().has_on_3_object()(la,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circular_arc_3 & a1, 
                const typename SK::Circular_arc_3 & a2, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK::Circular_arc_3 Circular_arc_3;
      typedef std::vector<CGAL::Object> solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      if(non_oriented_equal<SK>(a1.supporting_circle(), a2.supporting_circle())) {
        if(a1.rep().is_full()) {
          *res++ = make_object(a2); 
          //return res;
        }
        else if(a2.rep().is_full()) {
          *res++ = make_object(a1); 
          //return res;
        } else {
          bool t2_in_a1 = SK().has_on_3_object()(a1,a2.target(),true);
          bool s2_in_a1 = SK().has_on_3_object()(a1,a2.source(),true);
          if(t2_in_a1 && s2_in_a1) {
            bool t1_in_a2 = SK().has_on_3_object()(a2,a1.target(),true);
            bool s1_in_a2 = SK().has_on_3_object()(a2,a1.source(),true);
            if(t1_in_a2 && s1_in_a2) {
              const Comparison_result comp = 
                SK().compare_xyz_3_object()(a1.source(), a2.source());
              if(comp < 0) {
                if(a1.source() == a2.target()) {
                  *res++ = make_object(std::make_pair(a1.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a1.source(),a2.target());
	          *res++ = make_object(arc);
                }
                if(a2.source() == a1.target()) {
                  *res++ = make_object(std::make_pair(a2.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a2.source(),a1.target());
	          *res++ = make_object(arc);
                }
              } else if(comp > 0) {
                if(a2.source() == a1.target()) {
                  *res++ = make_object(std::make_pair(a2.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a2.source(),a1.target());
	          *res++ = make_object(arc);
                }
                if(a1.source() == a2.target()) {
                  *res++ = make_object(std::make_pair(a1.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a1.source(),a2.target());
	          *res++ = make_object(arc);
                } 
              } else { 
                *res++ = make_object(a1);
              }
            } else {
              *res++ = make_object(a2);
            }
          } else if(t2_in_a1) {
            if(a1.source() == a2.target()) 
              *res++ = make_object(std::make_pair(a1.source(),1u));
            else {
              const Circular_arc_3 & arc =
	        Circular_arc_3(a1.supporting_circle(),a1.source(),a2.target());
	      *res++ = make_object(arc);
            } //return res;
          } else if(s2_in_a1) {
            if(a2.source() == a1.target()) {
              *res++ = make_object(std::make_pair(a2.source(),1u));
            } else {
              const Circular_arc_3 & arc =
	        Circular_arc_3(a1.supporting_circle(),a2.source(),a1.target());
	      *res++ = make_object(arc);
            }
          } else if(SK().has_on_3_object()(a2,a1.source(),true)) {
              *res++ = make_object(a1);
          } 
        }
      } else {
        solutions_container solutions;

        SK().intersect_3_object()(a1.supporting_circle(), a2.supporting_circle(), 
                                  std::back_inserter(solutions) );
        if(solutions.size() == 0) return res;
        if(solutions.size() == 1) {
          Solution sol;
          assign(sol, solutions[0]);
          if(SK().has_on_3_object()(a1,sol.first,true) &&
             SK().has_on_3_object()(a2,sol.first,true))
            *res++ = solutions[0];
        } else {
          Solution sol1, sol2;
          assign(sol1, solutions[0]);
          assign(sol2, solutions[1]);
          if(SK().has_on_3_object()(a1,sol1.first,true) &&
             SK().has_on_3_object()(a2,sol1.first,true))
            *res++ = solutions[0];
          if(SK().has_on_3_object()(a1,sol2.first,true) &&
             SK().has_on_3_object()(a2,sol2.first,true))
            *res++ = solutions[1];
        }
      }
      return res;
    }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H
