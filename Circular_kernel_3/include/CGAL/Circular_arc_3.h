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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_ARC_3_H
#define CGAL_CIRCULAR_ARC_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/result_of.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
  template <class SK> 
    class Circular_arc_3
    : public SK::Kernel_base::Circular_arc_3
  {
    
    typedef typename SK::RT                          RT;
    typedef typename SK::FT                          FT;
    typedef typename SK::Line_3                      Line_3;
    typedef typename SK::Point_3                     Point_3;
    typedef typename SK::Plane_3                     Plane_3;
    typedef typename SK::Circle_3                    Circle_3;
    typedef typename SK::Sphere_3                    Sphere_3;
    typedef typename SK::Segment_3                   Segment_3;
    typedef typename SK::Circular_arc_point_3        Circular_arc_point_3;
    typedef typename SK::Kernel_base::Circular_arc_3 RCircular_arc_3; 
   
  
  public:
    typedef  RCircular_arc_3 Rep;
    typedef  SK   R; 
    
    const Rep& rep() const
      {
	return *this;
      }
    
    Rep& rep()
      {
	return *this;
      }
    
    Circular_arc_3()
      : RCircular_arc_3(typename R::Construct_circular_arc_3()())
      {}

    Circular_arc_3(const Circle_3& c, 
               const Circular_arc_point_3& s, 
               const Circular_arc_point_3& t)
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(c,s,t))
      {}

    explicit Circular_arc_3(const Circle_3& c)
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(c))
      {}
        
    Circular_arc_3(const Circle_3& c,const Circular_arc_point_3& pt)
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(c,pt))
      {}        

    // Not Documented
    Circular_arc_3(const Circle_3 &c, 
                   const Sphere_3 &s1, bool less_xyz_s1,
                   const Sphere_3 &s2, bool less_xyz_s2) 
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(c,s1,less_xyz_s1,
                                                                 s2,less_xyz_s2))
      {}

    // Not Documented
    Circular_arc_3(const Sphere_3 &s1, bool less_xyz_s1,
                   const Sphere_3 &s2, bool less_xyz_s2,
                   const Circle_3 &c) 
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(s1,less_xyz_s1,
                                                               s2,less_xyz_s2,
                                                               c))
      {}

    // Not Documented
    Circular_arc_3(const Circle_3 &c, 
                   const Plane_3 &p1, bool less_xyz_p1,
                   const Plane_3 &p2, bool less_xyz_p2) 
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(c,p1,less_xyz_p1,
                                                                 p2,less_xyz_p2))
      {}

    // Not Documented
    Circular_arc_3(const Plane_3 &p1, bool less_xyz_p1,
                   const Plane_3 &p2, bool less_xyz_p2,
                   const Circle_3 &c) 
      : RCircular_arc_3(typename R::Construct_circular_arc_3()(p1,less_xyz_p1,
                                                               p2,less_xyz_p2,
                                                               c))
      {}

	  Circular_arc_3(const Point_3 &start,
	                 const Point_3 &middle,
	                 const Point_3 &end)
	    : RCircular_arc_3(typename 
			      R::Construct_circular_arc_3()(start, middle, end)) 
	  {}

    Circular_arc_3(const RCircular_arc_3 &a)
     : RCircular_arc_3(a)
      {}

    typename cpp11::result_of<typename R::Construct_circular_source_vertex_3(Circular_arc_3)>::type
    source() const
    {
      return typename R::Construct_circular_source_vertex_3()(*this);
    }

    typename cpp11::result_of<typename R::Construct_circular_target_vertex_3(Circular_arc_3)>::type
    target() const
    {
      return typename R::Construct_circular_target_vertex_3()(*this);
    }

    typename cpp11::result_of<typename R::Construct_circle_3(Circular_arc_3)>::type
    supporting_circle() const
    {
      return typename R::Construct_circle_3()(*this);
    }

    Sphere_3 diametral_sphere() const
    {
      return typename R::Construct_sphere_3()(*this);
    }

    Point_3 center() const
    {
      return typename R::Construct_sphere_3()(*this).center();
    }

    FT squared_radius() const
    {
      return typename R::Construct_sphere_3()(*this).squared_radius();
    }

    Plane_3 supporting_plane() const
    {
      return typename R::Construct_plane_3()(*this);
    }

    Bbox_3 bbox() const
    { return typename R::Construct_bbox_3()(*this); }

  };

  template < typename SK >
  inline
  bool
  operator==(const Circular_arc_3<SK> &p,
	     const Circular_arc_3<SK> &q)
  {
    return SK().equal_3_object()(p, q);
  }
  
  template < typename SK >
  inline
  bool
  operator!=(const Circular_arc_3<SK> &p,
	     const Circular_arc_3<SK> &q)
  {
    return ! (p == q);
  }
  
  template < typename SK >
  std::ostream &
  operator<<(std::ostream & os, const Circular_arc_3<SK> &a)
  {
    return os << a.supporting_circle() << " "
	      << a.source() << " "
	      << a.target() << " ";
  }
  
  template < typename SK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_3<SK> &a)
  {
    typename SK::Circle_3 s;
    typename SK::Circular_arc_point_3 p1;
    typename SK::Circular_arc_point_3 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Circular_arc_3<SK>(s, p1, p2);
    return is;
  }

}

#endif
