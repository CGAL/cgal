// Copyright (c) 2001-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_MESH_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_TRAITS_2_H

namespace CGAL {

template <class K>
class Delaunay_mesh_traits_2 : public K
{
  double B;
public:

  Delaunay_mesh_traits_2(const double bound = 0.125) : B(bound) {};

  typedef double Quality;

  inline
  double bound() const { return B; };

  inline 
  void set_bound(const double bound) { B = bound; };

  class Is_bad
  {
  protected:
    const double B;
  public:
    typedef typename K::Point_2 Point_2;
    typedef Delaunay_mesh_traits_2<K> Traits;
      
    Is_bad(const double bound) : B(bound) {};
      
    bool operator()(const Point_2& pa,
		    const Point_2& pb,
		    const Point_2& pc,
		    Quality& q) const
    {
      typedef typename K::Point_2 Point_2;
      typedef typename K::Triangle_2 Triangle_2;
      typedef typename K::Compute_area_2 Compute_area_2;
      typedef typename K::Compute_squared_distance_2
	Compute_squared_distance_2;
      typedef typename K::FT FT;

      K k;
      Compute_area_2 area_2 = k.compute_area_2_object();
      Compute_squared_distance_2 squared_distance = 
	k.compute_squared_distance_2_object();

      Triangle_2 t = k.construct_triangle_2_object()(pa,pb,pc);
      double area = 2*CGAL::to_double(area_2(t));
      area=area*area;

      double
	a = CGAL::to_double(squared_distance(pb, pc)),
	b = CGAL::to_double(squared_distance(pc, pa)),
	c = CGAL::to_double(squared_distance(pa, pb));

      if(a<b)
	if(a<c)
	  q = area/(b*c);
	else
	  q = area/(a*b);
      else
	if(b<c)
	  q = area/(a*c);
	else
	  q = area/(a*b);

      return ( q < B );
    };
  };

  Is_bad is_bad_object() const
    { return Is_bad(B); }
};

}; // end namespace CGAL

#endif
