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

#ifndef CGAL_DELAUNAY_MESH_LOCAL_SIZE_CRITERIA_2_H
#define CGAL_DELAUNAY_MESH_LOCAL_SIZE_CRITERIA_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Delaunay_mesh_size_criteria_2.h>

namespace CGAL {

template <class CDT>
class Delaunay_mesh_local_size_criteria_2
  : public Delaunay_mesh_size_criteria_2<CDT>
{
public:
  typedef Delaunay_mesh_size_criteria_2<CDT> Base;
  typedef typename Base::Base PreviousBase;
  typedef typename CDT::Geom_traits Geom_traits;
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_2 Point;
  typedef typename Geom_traits::Segment_2 Segment;

  typedef typename CDT::Face_handle Face_handle;

private:
  bool local;
  Segment _s;
  Geom_traits traits;

public:
  Delaunay_mesh_local_size_criteria_2(const double aspect_bound = 0.125, 
                                      const double size_bound = 0,
                                      const bool is_local_size = false,
                                      const Segment s = Segment(),
                                      const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound, size_bound), local(is_local_size), _s(s)
    , traits(traits) {}

  inline
  Segment segment() const { return _s; }

  inline 
  void set_segment(const Segment s) { _s = s; }

  inline
  bool is_local_size() const { return local; }

  inline 
  void set_local_size(bool local_size) { local = local_size; }

  class Is_bad: public Base::Is_bad
  {
  public:
    typedef typename Base::Is_bad Baseclass;
    typedef typename Baseclass::Point_2 Point_2;
    typedef typename CDT::Geom_traits::Segment_2 Segment_2;
    typedef typename CDT::Geom_traits::Triangle_2 Triangle_2;
    typedef typename Base::Base PreviousBase;
    typedef typename PreviousBase::Is_bad PreviousBaseclass;
    typedef Geom_traits Traits;

    typedef typename Base::Quality Quality;

  private:
    const bool local;
    const Segment_2 s;

  public:
    Is_bad(const double aspect_bound,
	   const double size_bound,
	   const bool l,
	   const Segment_2 _s,
           const Geom_traits& traits)
      : Base::Is_bad(aspect_bound, size_bound, traits), local(l), s(_s) {}

    Mesh_2::Face_badness operator()(Quality q) const
    {
      return Base::Is_bad::operator()(q);
    }
    
    Mesh_2::Face_badness operator()(const Face_handle& fh,
				    Quality& q) const
    {
      if(!local)
	return Baseclass::operator()(fh,q);
      else
	{
	  typename Geom_traits::Do_intersect_2 do_intersect = 
	    Geom_traits().do_intersect_2_object();
	  
	  Mesh_2::Face_badness is_non_locally_bad = 
	    Baseclass::operator()(fh,q);

          const Point_2& a = fh->vertex(0)->point();
          const Point_2& b = fh->vertex(1)->point();
          const Point_2& c = fh->vertex(2)->point();

	  if(do_intersect(Triangle_2(a,b,c), s))
	    return is_non_locally_bad;
	  else
	    if( q.sine() < this->B )
	      return Mesh_2::BAD;
	    else
	      return Mesh_2::NOT_BAD;
	}
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(this->bound(), this->size_bound(), local, segment(), traits); }
};

} //end namespace

#endif
