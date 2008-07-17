// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_3/include/CGAL/Circular_kernel_3/Circular_arc_3.h $
// $Id: Circular_arc_3.h 40627 2007-10-16 15:00:59Z sloriot $
//
// Author(s) : Loriot Sebastien <Sebastien.Loriot@sophia.inria.fr>
//                   Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                   Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//                   Pedro Machado    <tashimir@gmail.com>
//                   Julien Hazebrouck
//                   Damien Leroy

#ifndef CGAL_CIRCULAR_ARC_ON_REFERENCE_SPHERE_3_H
#define CGAL_CIRCULAR_ARC_ON_REFERENCE_SPHERE_3_H

//~ #include <CGAL/Circular_arc_3.h>
//~ #include <CGAL/utility.h>
//~ #include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_3.h>

namespace CGAL {
  template <class SK> 
    class Circular_arc_on_reference_sphere_3
    //~ : public SK::Kernel_base::Circular_arc_on_reference_sphere_3, public CGAL::Circular_arc_3<SK>
    : public SK::Kernel_base::Circular_arc_on_reference_sphere_3
  {
    
    typedef typename SK::RT                          RT;
    typedef typename SK::FT                          FT;
    typedef typename SK::Line_3                      Line_3;
    typedef typename SK::Point_3                     Point_3;
    typedef typename SK::Plane_3                     Plane_3;
    typedef typename SK::Circle_on_reference_sphere_3                    Circle_3;
    typedef typename SK::Sphere_with_radius_3        Sphere_3;
    typedef typename SK::Segment_3                   Segment_3;
    typedef typename SK::Circular_arc_point_on_reference_sphere_3        Circular_arc_point_3;
    typedef typename SK::Kernel_base::Circular_arc_on_reference_sphere_3 RCircular_arc_on_reference_sphere_3; 
    typedef typename CGAL::Circular_arc_3<SK>                  Circular_arc_3;
  
  public:
    typedef  RCircular_arc_on_reference_sphere_3 Rep;
    typedef  SK   R; 
    
    const Rep& rep() const
      {
	return *this;
      }
    
    Rep& rep()
      {
	return *this;
      }
    
    Circular_arc_on_reference_sphere_3(const RCircular_arc_on_reference_sphere_3 &a)
     : RCircular_arc_on_reference_sphere_3(a)
      {}      
      
    Circular_arc_on_reference_sphere_3()
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()())
      {}

    Circular_arc_on_reference_sphere_3(const Circle_3& c, 
               const Circular_arc_point_3& s, 
               const Circular_arc_point_3& t)
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c,s,t)){}

    Circular_arc_on_reference_sphere_3(const Circle_3& c, 
               const Point_3& s, 
               const Circular_arc_point_3& t)
      :RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c,s,t)){}

    Circular_arc_on_reference_sphere_3(const Circle_3& c, 
               const Circular_arc_point_3& s, 
               const Point_3& t)
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c,s,t)){}

    Circular_arc_on_reference_sphere_3(const Circle_3& c, 
               const Point_3& s, 
               const Point_3& t)
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c,s,t)){}

    Circular_arc_on_reference_sphere_3(const Circle_3& c)
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c)){}

    Circular_arc_on_reference_sphere_3(const Circle_3 &c, 
                   const Sphere_3 &s1, bool less_xyz_s1,
                   const Sphere_3 &s2, bool less_xyz_s2) 
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c,s1,less_xyz_s1,s2,less_xyz_s2)){}

    Circular_arc_on_reference_sphere_3(const Sphere_3 &s1, bool less_xyz_s1,
                   const Sphere_3 &s2, bool less_xyz_s2,
                   const Circle_3 &c) 
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(s1,less_xyz_s1,s2,less_xyz_s2,c)){}

    Circular_arc_on_reference_sphere_3(const Circle_3 &c, 
                   const Plane_3 &p1, bool less_xyz_p1,
                   const Plane_3 &p2, bool less_xyz_p2) 
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(c,p1,less_xyz_p1,p2,less_xyz_p2)){}

    Circular_arc_on_reference_sphere_3(const Plane_3 &p1, bool less_xyz_p1,
                   const Plane_3 &p2, bool less_xyz_p2,
                   const Circle_3 &c) 
      : RCircular_arc_on_reference_sphere_3(typename R::Construct_circular_arc_on_reference_sphere_3()(p1,less_xyz_p1,p2,less_xyz_p2,c)){}

    typename Qualified_result_of
    <typename R::Compute_reference_sphere_3,Circular_arc_on_reference_sphere_3>::type
    reference_sphere() const
    {
      return typename R::Compute_reference_sphere_3()(*this);
    }
  };

}

#endif
