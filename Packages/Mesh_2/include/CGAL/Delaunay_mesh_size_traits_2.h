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

#ifndef CGAL_DELAUNAY_MESH_SIZE_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_SIZE_TRAITS_2_H

#include <CGAL/Delaunay_mesh_traits_2.h>
#include <utility>
#include <ostream>

namespace CGAL {

template <class K>
class Delaunay_mesh_size_traits_2 : public virtual Delaunay_mesh_traits_2<K>
{
protected:
  double sizebound;

public:
  typedef Delaunay_mesh_traits_2<K> Base;

  Delaunay_mesh_size_traits_2(const double aspect_bound = 0.125, 
		     const double size_bound = 0)
    : Base(aspect_bound), sizebound(size_bound) {};

  inline
  double size_bound() const { return sizebound; };

  inline
  void set_size_bound(const double sb) { sizebound = sb; };

  // first: squared_minimum_sine
  // second: size
  struct Quality : public std::pair<double, double>
  {
    typedef std::pair<double, double> Base;

    Quality() : Base() {};
    Quality(double _sine, double _size) : Base(_sine, _size) {};

    double size() const { return second; }
    double sine() const { return first; }

    bool operator<(const Quality& q) const
      {
	if( size() > q.size() )
	  return true;
	else if( size() == q.size() )
	  return( sine() < q.sine() );
	else
	  return false;
      }

    friend std::ostream& operator<<(std::ostream& out, const Quality& q)
      {
	return out << "(size=" << q.size() << ", sine=" << q.sine() <<
	  ")";
      }
  };

  class Is_bad: public Base::Is_bad
  {
  protected:
    const double SB; // squared size bound on edge length
  public:
    typedef typename Base::Is_bad::Point_2 Point_2;
    typedef typename Base::Is_bad::Traits Traits;

    Is_bad(const double aspect_bound,
	   const double size_bound)
      : Base::Is_bad(aspect_bound), SB(size_bound * size_bound) {};

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
      double max_lenght; // squared max edge length
      
      if(a<b)
	{
	  if(a<c) min_sine = area/(b*c);
	  else min_sine = area/(a*b);
	  
	  if(b<c) max_lenght = c;
	  else max_lenght = b;
	}
      else
	{
	  if(b<c) min_sine = area/(a*c);
	  else min_sine = area/(a*b);
	  
	  if(a<c) max_lenght = c;
	  else max_lenght = a;
	}
      
      q.first = min_sine;
      q.second = max_lenght;
      
      if( min_sine < B ) return true;
      if( SB == 0 ) return false;
      return ( max_lenght > SB );
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(bound(), size_bound()); }
};

} //end namespace

#endif
