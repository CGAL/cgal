// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_MESH_SIZE_CRITERIA_2_H
#define CGAL_DELAUNAY_MESH_SIZE_CRITERIA_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>

namespace CGAL {

template <class CDT>
class Delaunay_mesh_size_criteria_2 : 
    public virtual Delaunay_mesh_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits Geom_traits;
  double sizebound;

public:
  typedef Delaunay_mesh_criteria_2<CDT> Base;

  Delaunay_mesh_size_criteria_2(const double aspect_bound = 0.125, 
                                const double size_bound = 0,
                                const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound, traits), sizebound(size_bound) {}

  inline
  double size_bound() const { return sizebound; }

  inline
  void set_size_bound(const double sb) { sizebound = sb; }

  // first: squared_minimum_sine
  // second: size
  struct Quality : public std::pair<double, double>
  {
    typedef std::pair<double, double> Base;

    Quality() : Base() {};
    Quality(double _sine, double _size) : Base(_sine, _size) {}

    const double& size() const { return second; }
    const double& sine() const { return first; }

    // q1<q2 means q1 is prioritised over q2
    // ( q1 == *this, q2 == q )
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

    std::ostream& operator<<(std::ostream& out) const
    {
      return out << "(size=" << size()
		 << ", sine=" << sine() << ")";
    }
  };

  class Is_bad: public Base::Is_bad
  {
  protected:
    const double squared_size_bound; // squared size bound on edge length
  public:
    typedef typename Base::Is_bad::Point_2 Point_2;

    Is_bad(const double aspect_bound,
	   const double size_bound,
           const Geom_traits& traits)
      : Base::Is_bad(aspect_bound, traits),
        squared_size_bound(size_bound * size_bound) {}

    Mesh_2::Face_badness operator()(const Quality q) const
    {
      if( q.size() > 1 )
	return Mesh_2::IMPERATIVELY_BAD;
      if( q.sine() < this->B )
	return Mesh_2::BAD;
      else
	return Mesh_2::NOT_BAD;
    }

    Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
				    Quality& q) const
    {
      typedef typename CDT::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2
	Compute_squared_distance_2;

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
      
      double max_sq_length; // squared max edge length
      double second_max_sq_length;
      
      if(a<b)
	{
	  if(b<c) {
	    max_sq_length = c;
	    second_max_sq_length = b;
	  }
	  else { // c<=b
	    max_sq_length = b;
	    second_max_sq_length = ( a < c ? c : a );
	  }
	}
      else // b<=a
	{
	  if(a<c) {
	    max_sq_length = c;
	    second_max_sq_length = a;
	  }
	  else { // c<=a
	    max_sq_length = a;
	    second_max_sq_length = ( b < c ? c : b );
	  }
	}

      q.second = 0;
      if( squared_size_bound != 0 )
	{
          //	  std::cerr << squared_size_bound << std::endl;
	  q.second = max_sq_length / squared_size_bound;
	    // normalized by size bound to deal
	    // with size field
	  if( q.size() > 1 )
	    {
	      q.first = 1; // (do not compute sine)
	      return Mesh_2::IMPERATIVELY_BAD;
	    }
	}

      Compute_area_2 area_2 = traits.compute_area_2_object();

      double area = 2*CGAL::to_double(area_2(pa, pb, pc));

      q.first = (area * area) / (max_sq_length * second_max_sq_length); // (sine)
      
      if( q.sine() < this->B )
	return Mesh_2::BAD;
      else
	return Mesh_2::NOT_BAD;
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(this->bound(), size_bound(), 
                  this->traits /* from the bad class */); }
};

} // end namespace CGAL

#endif
