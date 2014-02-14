// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/Mesh_3/mesh_standard_cell_criteria.h $
// $Id: mesh_standard_cell_criteria.h 53504 2009-12-18 17:15:58Z stayeb $
//
//
// Author(s)     : Bogdanov Mikhail
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_PERIODIC_MESH_STANDARD_CELL_CRITERIA_H
#define CGAL_PERIODIC_MESH_STANDARD_CELL_CRITERIA_H


#include <CGAL/Mesh_3/mesh_standard_cell_criteria.h>


namespace CGAL {

namespace Mesh_3 {

namespace Periodic_mesh_3 {

template <typename Tr, typename Visitor_>
class Cell_radius_edge_criterion
  : public Abstract_criterion<Tr,Visitor_>
{
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Cell_radius_edge_criterion<Tr, Visitor_> Self;

public:
  // Constructor
  Cell_radius_edge_criterion(const Tr& tr_, const FT& radius_edge_bound)
    : sq_radius_edge_bound_(radius_edge_bound*radius_edge_bound), tr(tr_) 
  { };

  // Destructor
  ~Cell_radius_edge_criterion() { };


protected:
  virtual void do_accept(Visitor_& v) const
  {
    v.visit(*this);
  }

  virtual Self* do_clone() const
  {
    // Call copy ctor on this
    return new Self(*this);
  }

  virtual Badness do_is_bad(const Cell_handle& ch) const
  {
    typedef typename Tr::Point Point_3;
    typedef typename Tr::Geom_traits Geom_traits;
    typedef typename Geom_traits::Compute_squared_radius_3 Radius;
    typedef typename Geom_traits::Compute_squared_distance_3 Distance;

    const Point_3 p = tr.point(ch, 0);
    const Point_3 q = tr.point(ch, 1);
    const Point_3 r = tr.point(ch, 2);
    const Point_3 s = tr.point(ch, 3);

    Radius radius = Geom_traits().compute_squared_radius_3_object();
    Distance distance = Geom_traits().compute_squared_distance_3_object();

    const FT size = radius(p, q, r, s);

    FT min_sq_length = distance(p, q);
    min_sq_length = (CGAL::min)(min_sq_length, distance(p, r));
    min_sq_length = (CGAL::min)(min_sq_length, distance(p, s));
    min_sq_length = (CGAL::min)(min_sq_length, distance(q, r));
    min_sq_length = (CGAL::min)(min_sq_length, distance(q, s));
    min_sq_length = (CGAL::min)(min_sq_length, distance(r, s));

    if ( size > min_sq_length*sq_radius_edge_bound_  )
    {
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
      std::cerr << "bad cell (radius-edge bound): radius-edge["
                << size/min_sq_length << "] bound[" << sq_radius_edge_bound_
                << "]\n" ;
#endif
      return Badness(Quality( (sq_radius_edge_bound_*min_sq_length)/size ));
    }
    else
      return Badness();
  }

private:
  FT sq_radius_edge_bound_;

  const Tr& tr;
};  // end class Cell_radius_edge_criterion



template <typename Tr, typename Visitor_>
class Cell_radius_criterion
  : public Abstract_criterion<Tr, Visitor_>
{
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Abstract_criterion<Tr, Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Cell_radius_criterion<Tr, Visitor_> Self;

public:
  // Constructor
  Cell_radius_criterion(const Tr& tr_, const FT& radius_bound)
    : sq_radius_bound_(radius_bound*radius_bound), tr(tr_)   { };

  // Destructor
  ~Cell_radius_criterion() { };

protected:
  virtual void do_accept(Visitor_& v) const
  {
    v.visit(*this);
  }

  virtual Self* do_clone() const
  {
    // Call copy ctor on this
    return new Self(*this);
  }

  virtual Badness do_is_bad(const Cell_handle& ch) const
  {
    typedef typename Tr::Point Point_3;
    typedef typename Tr::Geom_traits Geom_traits;
    typedef typename Geom_traits::Compute_squared_radius_3 Radius;

	const Point_3 p = tr.point(ch, 0);
	const Point_3 q = tr.point(ch, 1);
	const Point_3 r = tr.point(ch, 2);
	const Point_3 s = tr.point(ch, 3);

    Radius radius = Geom_traits().compute_squared_radius_3_object();

    const FT size = radius(p, q, r, s);

    if ( size > sq_radius_bound_ )
    {
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
      std::cerr << "bad cell (radius bound): size[" << size
                << "] bound[" << sq_radius_bound_ << "]\n" ;
#endif
      return Badness(Quality(sq_radius_bound_/size));
    }
    else
      return Badness();
  }

private:
  FT sq_radius_bound_;

  const Tr& tr;
};  // end class Cell_radius_criterion

}  // end namespace Periodic_mesh_3

}  // end namespace Mesh_3

}  // end namespace CGAL


#endif // CGAL_PERIODIC_MESH_STANDARD_CELL_CRITERIA_H
