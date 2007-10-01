// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Pedro Machado    <tashimir@gmail.com>
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_CIRCLE_3_H
#define CGAL_CIRCLE_3_H

namespace CGAL {
  template <class SK> 
    class Circle_3
    : public SK::Kernel_base::Circle_3
  {
    
    typedef typename SK::RT                    RT;
    typedef typename SK::FT                    FT;
    typedef typename SK::Point_3               Point_3;
    typedef typename SK::Plane_3               Plane_3;
    typedef typename SK::Sphere_3              Sphere_3;
    typedef typename SK::Vector_3              Vector_3;
    typedef typename SK::Direction_3           Direction_3;
    typedef typename SK::Kernel_base::Circle_3 RCircle_3; 
   
  
  public:
    typedef  RCircle_3 Rep;
    typedef  SK   R; 
    
    const Rep& rep() const
      {
	return *this;
      }
    
    Rep& rep()
      {
	return *this;
      }
    
    Circle_3()
      : RCircle_3(typename R::Construct_circle_3()())
      {}

    Circle_3(const Point_3& c, const FT& sr, const Plane_3& p)
      : RCircle_3(typename R::Construct_circle_3()(c,sr,p))
      {}

    Circle_3(const Point_3& c, const FT& sr, const Direction_3& d) 
      : RCircle_3(typename R::Construct_circle_3()(c,sr,d))
      {}

    Circle_3(const Point_3& c, const FT& sr, const Vector_3& v) 
      : RCircle_3(typename R::Construct_circle_3()(c,sr,v))
      {}

    Circle_3(const Sphere_3& s1, const Sphere_3& s2)
      : RCircle_3(typename R::Construct_circle_3()(s1,s2))
      {}

    Circle_3(const Sphere_3& s, const Plane_3& p)
      : RCircle_3(typename R::Construct_circle_3()(s,p))
      {}

    Circle_3(const Plane_3& p, const Sphere_3& s)
      : RCircle_3(typename R::Construct_circle_3()(p,s))
      {}

    Circle_3(const RCircle_3& r)
      : RCircle_3(r)
      {}

    typename Qualified_result_of
    <typename R::Construct_diametral_sphere_3, Circle_3>::type
    //const Sphere_3 &
    diametral_sphere() const
    {
      return typename R::Construct_diametral_sphere_3()(*this);
    }

    Point_3 center() const
    {
      return typename R::Construct_diametral_sphere_3()(*this).center();
    }

    FT squared_radius() const
    {
      return typename R::Construct_diametral_sphere_3()(*this).squared_radius();
    }

    typename Qualified_result_of
    <typename R::Construct_supporting_plane_3, Circle_3>::type
    //const Plane_3 &
    supporting_plane() const
    {
      return typename R::Construct_supporting_plane_3()(*this);
    }

    typename Qualified_result_of
    <typename R::Construct_supporting_sphere_3, Circle_3>::type
    //const Plane_3 &
    supporting_sphere() const
    {
      return typename R::Construct_supporting_sphere_3()(*this);
    }
    
    Bbox_3 bbox() const
    { return typename R::Construct_bbox_3()(*this); }

  };

  template < typename SK >
  inline
  bool
  operator==(const Circle_3<SK> &p,
	     const Circle_3<SK> &q)
  {
    return SK().equal_3_object()(p, q);
  }
  
  template < typename SK >
  inline
  bool
  operator!=(const Circle_3<SK> &p,
	     const Circle_3<SK> &q)
  {
    return ! (p == q);
  }
  
  template < typename SK >
  std::ostream &
  operator<<(std::ostream & os, const Circle_3<SK> &c)
  {
    
    return os << c.supporting_plane() << " "
	      << c.diametral_sphere() << " ";
  }
  
  template < typename SK >
  std::istream &
  operator>>(std::istream & is, Circle_3<SK> &c)
  {
    typename SK::Plane_3 p;
    typename SK::Sphere_3 s;

    is >> p >> s ;
    if (is)
      c = Circle_3<SK>(p, s);
    return is;
  }

}

#endif
