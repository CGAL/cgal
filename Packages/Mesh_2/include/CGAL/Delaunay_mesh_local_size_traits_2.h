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

#ifndef CGAL_DELAUNAY_MESH_LOCAL_SIZE_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_LOCAL_SIZE_TRAITS_2_H

#include <CGAL/Delaunay_mesh_size_traits_2.h>

namespace CGAL {

template <class K>
class Delaunay_mesh_local_size_traits_2 : public Delaunay_mesh_size_traits_2<K>
{
public:
  typedef Delaunay_mesh_size_traits_2<K> Base;
  typedef typename Base::Base PreviousBase;
  typedef typename K::FT FT;
  typedef typename K::Point_2 Point;

private:
  bool local;
  Point _p;

public:
  Delaunay_mesh_local_size_traits_2(const double aspect_bound = 0.125, 
			   const double size_bound = 0,
			   const bool is_local_size = false,
			   const Point p = Point())
    : Base(aspect_bound, size_bound), local(is_local_size), _p(p) {};

  inline
  Point point() const { return _p; };

  inline 
  void set_point(const Point p) { _p = p; };

  inline
  bool is_local_size() const { return local; };

  inline 
  void set_local_size(bool local_size) { local = local_size; };

  class Is_bad: public Base::Is_bad
  {
  public:
    typedef typename Base::Is_bad Baseclass;
    typedef typename Baseclass::Point_2 Point_2;
    typedef typename Baseclass::Traits Traits;
    typedef typename Base::Base PreviousBase;
    typedef typename PreviousBase::Is_bad PreviousBaseclass;

    typedef typename Base::Quality Quality;

  private:
    const bool local;
    const Point_2 p;

  public:
    Is_bad(const double aspect_bound,
	   const double size_bound,
	   const bool l,
	   const Point_2 _p)
      : Base::Is_bad(aspect_bound, size_bound), local(l), p(_p) {};

    bool operator()(const Point_2& a,
		    const Point_2& b,
		    const Point_2& c,
		    Quality& q) const
    {
      if(!local)
	return Baseclass::operator()(a,b,c,q);
      else
	{
	  typename Traits::Orientation_2 orient = 
	    Traits().orientation_2_object();
	  
	  bool is_non_locally_bad = Baseclass::operator()(a,b,c,q);

	  if( q.first < B )
	    return true;

	  Orientation 
	    o1 = orient(a,b,p),
	    o2 = orient(b,c,p),
	    o3 = orient(c,a,p);
	  
	  if((o1==o2) && (o2==o3))
	    return is_non_locally_bad;
	  else
	    return false;
	};
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(bound(), size_bound(), local, point()); };
};

}; //end namespace

#endif
