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

#ifndef CGAL_DELAUNAY_MESH_AREA_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_AREA_TRAITS_2_H

#include <CGAL/Delaunay_mesh_size_traits_2.h>

namespace CGAL {

template <class K>
class Delaunay_mesh_area_traits_2 : public virtual Delaunay_mesh_traits_2<K>,
				    private Delaunay_mesh_size_traits_2<K>
/* This class "is a" Delaunay_mesh_traits_2<K> and is implemented by
   Delaunay_mesh_size_traits_2<K>. Delaunay_mesh_traits_2<K> is a
   virtual base class of Delaunay_mesh_size_traits_2<K>. */
{
public:
  typedef Delaunay_mesh_traits_2<K> Base;
  typedef Delaunay_mesh_size_traits_2<K> Private_base;

  typedef typename Delaunay_mesh_size_traits_2<K>::Quality Quality;

  Delaunay_mesh_area_traits_2(const double aspect_bound = 0.125, 
			      const double area_bound = 0)
    : Private_base(aspect_bound, area_bound) {};

  inline
  double area_bound() const { return sizebound; };

  inline
  void set_area_bound(const double ab) { sizebound = ab; };

  class Is_bad: public Private_base::Is_bad
  {
  public:
    typedef typename Private_base::Is_bad Is_bad_base;

    typedef typename K::Point_2 Point_2;
    typedef typename K::Triangle_2 Triangle_2;

    Is_bad(const double aspect_bound,
	   const double area_bound)
      : Is_bad_base(aspect_bound, area_bound) {};

    bool operator()(const typename Base::Point_2& pa,
		    const typename Base::Point_2& pb,
		    const typename Base::Point_2& pc,
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
      area=area*area; // squared area

      double
	a = CGAL::to_double(squared_distance(pb, pc)),
	b = CGAL::to_double(squared_distance(pc, pa)),
	c = CGAL::to_double(squared_distance(pa, pb));

      double min_sine; // squared minimum sine

      if(a<b)
	if(a<c)
	  min_sine = area/(b*c);
	else
	  min_sine = area/(a*b);
      else
	if(b<c)
	  min_sine = area/(a*c);
	else
	  min_sine = area/(a*b);
      
      q.first = min_sine;
      q.second = area;

      if( min_sine < B ) return true;
      if( SB == 0 ) return false;
      return ( area > SB );
    };
  };

  Is_bad is_bad_object() const
  { return Is_bad(bound(), area_bound()); }
};

} //end namespace

#endif
