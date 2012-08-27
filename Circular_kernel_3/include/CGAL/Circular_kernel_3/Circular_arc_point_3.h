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
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_3_H

#include <iostream>

//#include <CGAL/global_functions_on_roots_and_polynomials_2_2.h> 
// fixme, devrait
// appeler fonction de global_functions_on_circular_arcs

namespace CGAL {
  namespace internal {

template <class SK >
class Circular_arc_point_3
{
  typedef typename SK::FT                         FT;
  typedef typename SK::Root_of_2                  Root_of_2;
  typedef typename SK::Point_3                    Point_3;
  typedef typename SK::Algebraic_kernel           Algebraic_kernel;
  typedef typename Algebraic_kernel::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
  typedef typename Algebraic_kernel::Polynomial_1_3             Polynomial_1_3;
  typedef typename Algebraic_kernel::Polynomials_for_line_3     Polynomials_for_line_3;
  typedef typename SK::Line_3                     Line_3;
  typedef typename SK::Plane_3                    Plane_3;
  typedef typename SK::Sphere_3                   Sphere_3;
  typedef typename SK::Circle_3                   Circle_3;
  typedef typename SK::Root_for_spheres_2_3       Root_for_spheres_2_3;

  typedef Root_for_spheres_2_3  Rep__;
  typedef typename SK::template Handle<Rep__>::type  Base;

  Base base;

public: 

  Circular_arc_point_3() {}
	
  Circular_arc_point_3(const Root_of_2 & x,
                       const Root_of_2 & y,
                       const Root_of_2 & z)
  : base(x,y,z){}

  Circular_arc_point_3(const Root_for_spheres_2_3 & np)
  : base(np){}

  Circular_arc_point_3(const Point_3 & p)
  : base(p.x(),p.y(),p.z()){}

  Circular_arc_point_3(const Sphere_3 &s1, 
                       const Sphere_3 &s2,
                       const Sphere_3 &s3,
                       const bool less_xyz = true) {
    std::vector<Object> sols;
    SK().intersect_3_object()(s1, s2, s3, std::back_inserter(sols));
    // s1,s2,s3 must intersect
    CGAL_kernel_precondition(sols.size() != 0);
    if(sols.size() == 1) {
      // the intersection must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
      CGAL_kernel_precondition(pair!=NULL);
      *this = pair->first.rep();
    } else {
      // the intersections must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[less_xyz?0:1]);
      CGAL_kernel_precondition(pair!=NULL);
      *this = pair->first.rep();
    } 
  }

  Circular_arc_point_3(const Plane_3 &p, 
                       const Sphere_3 &s1,
                       const Sphere_3 &s2,
                       const bool less_xyz = true) {
    std::vector<Object> sols;
    SK().intersect_3_object()(p, s1, s2, std::back_inserter(sols));
    // s1,s2,s3 must intersect
    CGAL_kernel_precondition(sols.size() != 0);
    if(sols.size() == 1) {
      // the intersection must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
      CGAL_kernel_precondition(pair!=NULL);
      *this = pair->first.rep();
    } else {
      // the intersections must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[less_xyz?0:1]);
      CGAL_kernel_precondition(pair!=NULL);
      *this = pair->first.rep();
    } 
  }

  Circular_arc_point_3(const Plane_3 &p1, 
                       const Plane_3 &p2,
                       const Sphere_3 &s,
                       const bool less_xyz = true) {
    std::vector<Object> sols;
    SK().intersect_3_object()(p1, p2, s, std::back_inserter(sols));
    // s1,s2,s3 must intersect
    CGAL_kernel_precondition(sols.size() != 0);
    if(sols.size() == 1) {
      // the intersection must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
      CGAL_kernel_precondition(pair!=NULL);      
      *this = pair->first.rep();
    } else {
      // the intersections must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[less_xyz?0:1]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    }
  }

  Circular_arc_point_3(const Line_3 &l,
                       const Sphere_3 &s,
                       const bool less_xyz = true) {
    std::vector<Object> sols;
    SK().intersect_3_object()(l, s, std::back_inserter(sols));
    // s1,s2,s3 must intersect
    CGAL_kernel_precondition(sols.size() != 0);
    if(sols.size() == 1) {
      // the intersection must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    } else {
      // the intersections must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[less_xyz?0:1]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    }
  }

  Circular_arc_point_3(const Circle_3 &c,
                       const Plane_3 &p,
                       const bool less_xyz = true) {
    std::vector<Object> sols;
    SK().intersect_3_object()(c, p, std::back_inserter(sols));
    // s1,s2,s3 must intersect
    CGAL_kernel_precondition(sols.size() != 0);
    if(sols.size() == 1) {
      // the intersection must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    } else {
      // the intersections must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[less_xyz?0:1]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    }
  }

  Circular_arc_point_3(const Circle_3 &c,
                       const Sphere_3 &s,
                       const bool less_xyz = true) {
    std::vector<Object> sols;
    SK().intersect_3_object()(c, s, std::back_inserter(sols));
    // s1,s2,s3 must intersect
    CGAL_kernel_precondition(sols.size() != 0);
    if(sols.size() == 1) {
      // the intersection must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    } else {
      // the intersections must be a point
      const std::pair<typename SK::Circular_arc_point_3, unsigned>* pair=
        object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[less_xyz?0:1]);
      CGAL_kernel_precondition(pair!=NULL);            
      *this = pair->first.rep();
    }
  }

  const Root_of_2 & x() const { return get(base).x(); }
  const Root_of_2 & y() const { return get(base).y(); }
  const Root_of_2 & z() const { return get(base).z(); }
	  
  const Root_for_spheres_2_3 & coordinates() const { return get(base); }

  const CGAL::Bbox_3 bbox() const {
    return get(base).bbox();
  }

  bool operator==(const Circular_arc_point_3 &) const;
  bool operator!=(const Circular_arc_point_3 &) const;

};

template < class SK >
CGAL_KERNEL_INLINE
bool
Circular_arc_point_3<SK>::operator==(const Circular_arc_point_3<SK> &t) const
{
  if (CGAL::identical(base, t.base))
      return true;
  return x() == t.x() &&
         y() == t.y() &&
         z() == t.z();
}

template < class SK >
CGAL_KERNEL_INLINE
bool
Circular_arc_point_3<SK>::operator!=(const Circular_arc_point_3<SK> &t) const
{
  return !(*this == t);
}
    
template < typename SK >
std::ostream &
print(std::ostream & os, const Circular_arc_point_3<SK> &p)
{
  return os << "CirclArcEndPoint_3(" << p.x() << ", " << p.y() << ')' << std::endl;
}

  } // namespace internal
} // namespace CGAL

#endif // CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_3_H
