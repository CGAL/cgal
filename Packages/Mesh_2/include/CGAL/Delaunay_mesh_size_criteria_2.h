// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_DELAUNAY_MESH_SIZE_CRITERIA_2_H
#define CGAL_DELAUNAY_MESH_SIZE_CRITERIA_2_H

#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>

namespace CGAL {

template <class CDT>
class Delaunay_mesh_size_criteria_2 : 
    public virtual Delaunay_mesh_criteria_2<CDT>
{
protected:
  double sizebound;

public:
  typedef Delaunay_mesh_criteria_2<CDT> Base;

  Delaunay_mesh_size_criteria_2(const double aspect_bound = 0.125, 
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

    Is_bad(const double aspect_bound,
	   const double size_bound)
      : Base::Is_bad(aspect_bound), SB(size_bound * size_bound) {};

    bool operator()(const typename CDT::Face_handle& fh,
		    Quality& q) const
    {
      typedef typename CDT::Geom_traits Geom_traits;
      typedef typename Geom_traits::Triangle_2 Triangle_2;
      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2
	Compute_squared_distance_2;
      typedef typename Geom_traits::Construct_triangle_2
	Construct_triangle_2;
      typedef typename Geom_traits::FT FT;

      Geom_traits traits; /** @warning traits with data!! */

      Compute_area_2 area_2 = 
        traits.compute_area_2_object();
      Compute_squared_distance_2 squared_distance = 
	traits.compute_squared_distance_2_object();
      Construct_triangle_2 triangle =
        traits.construct_triangle_2_object();

      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      Triangle_2 t = triangle(pa,pb,pc);
      double area = 2*CGAL::to_double(area_2(t));
      area=area*area;

      double
	a = CGAL::to_double(squared_distance(pb, pc)),
	b = CGAL::to_double(squared_distance(pc, pa)),
	c = CGAL::to_double(squared_distance(pa, pb));
      
      double min_sine; // squared minimum sine
      double max_length; // squared max edge length
      
      if(a<b)
	{
	  if(a<c) min_sine = area/(b*c);
	  else min_sine = area/(a*b);
	  
	  if(b<c) max_length = c;
	  else max_length = b;
	}
      else
	{
	  if(b<c) min_sine = area/(a*c);
	  else min_sine = area/(a*b);
	  
	  if(a<c) max_length = c;
	  else max_length = a;
	}
      
      q.first = min_sine;
      q.second = 1; // normalized by size bound to deal with
                    // size field
      
      if( min_sine < this->B ) return true;
      if( SB == 0 ) return false;
      q.second = max_length / SB;
      return ( max_length > SB ); /** @todo If max_length > SB
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(this->bound(), size_bound()); }
};

} // end namespace CGAL

#endif
