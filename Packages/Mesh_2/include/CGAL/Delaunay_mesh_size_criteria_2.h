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

    const double& size() const { return second; }
    const double& sine() const { return first; }

    // q1<q2 means q1 is prioritised over q2
    bool operator<(const Quality& q) const
    {
      if( size() > 1 )
	if( q.size() > 1 )
	  return ( size() > q.size() );
	else
	  return true; // *this is big but not q
      else
	if( q.size() >  1 )
	  return false; // q is big but not *this
      return( sine() < q.sine() );
    }

    friend std::ostream& operator<<(std::ostream& out, const Quality& q)
    {
      return out << "(size=" << q.size()
		 << ", sine=" << q.sine() << ")";
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
      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2
	Compute_squared_distance_2;
      typedef typename Geom_traits::Construct_triangle_2
	Construct_triangle_2;
      typedef typename Geom_traits::FT FT;

      Geom_traits traits; /** @warning traits with data!! */

      Compute_squared_distance_2 squared_distance = 
	traits.compute_squared_distance_2_object();

      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      double
	a = CGAL::to_double(squared_distance(pb, pc)),
	b = CGAL::to_double(squared_distance(pc, pa)),
	c = CGAL::to_double(squared_distance(pa, pb));
      
      double max_length; // squared max edge length
      double other_length_1, other_length_2;
      
      if(a<b)
	{
	  other_length_1 = a;
	  if(b<c) {
	    max_length = c;
	    other_length_2 = b;
	  }
	  else { // c<=b
	    max_length = b;
	    other_length_2 = c;
	  }
	}
      else // b<=a
	{
	  other_length_1 = b;
	  if(a<c) {
	    max_length = c;
	    other_length_2 = a;
	  }
	  else { // c<=a
	    max_length = a;
	    other_length_2 = c;
	  }
	}

      q.second = 0;
      if( SB != 0 )
	{
	  q.second = max_length / SB; // normalized by size bound to deal
				      // with size field
	  if( q.size() > 1 )
	    {
	      q.first = 1; // (do not compute sine)
	      return true;
	    }
	}

      Compute_area_2 area_2 = traits.compute_area_2_object();

      double area = 2*CGAL::to_double(area_2(pa, pb, pc));

      q.first = (area * area) / (other_length_1 * other_length_2); // (sine)
      
      return ( q.sine() < this->B );
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(this->bound(), size_bound()); }
};

} // end namespace CGAL

#endif
