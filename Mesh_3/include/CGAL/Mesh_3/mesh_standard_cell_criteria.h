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
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_STANDARD_CELL_CRITERIA_H
#define CGAL_MESH_3_MESH_STANDARD_CELL_CRITERIA_H


#include <CGAL/Mesh_3/mesh_standard_criteria.h>


namespace CGAL {

namespace Mesh_3 {


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
  Cell_radius_edge_criterion(const FT& radius_edge_bound)
    : sq_radius_edge_bound_(radius_edge_bound*radius_edge_bound) 
  {}

  // Destructor
  ~Cell_radius_edge_criterion() {}


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

    const Point_3& p = ch->vertex(0)->point();
    const Point_3& q = ch->vertex(1)->point();
    const Point_3& r = ch->vertex(2)->point();
    const Point_3& s = ch->vertex(3)->point();

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
  Cell_radius_criterion(const FT& radius_bound)
    : sq_radius_bound_(radius_bound*radius_bound)   {}

  // Destructor
  ~Cell_radius_criterion() {}

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

    const Point_3& p = ch->vertex(0)->point();
    const Point_3& q = ch->vertex(1)->point();
    const Point_3& r = ch->vertex(2)->point();
    const Point_3& s = ch->vertex(3)->point();

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

};  // end class Cell_radius_criterion



template <typename Tr>
class Cell_criterion_visitor
  : public Criterion_visitor<Tr, typename Tr::Cell_handle>
{
  typedef Criterion_visitor<Tr, typename Tr::Cell_handle> Base;
  typedef Cell_criterion_visitor<Tr> Self;

public:
  typedef Abstract_criterion<Tr, Self> Criterion;
  typedef typename Base::Quality Cell_quality;
  typedef typename Base::Badness Cell_badness;
  typedef typename Base::Handle Handle;
  typedef Handle Cell_handle;

  // Constructor
  Cell_criterion_visitor(const Cell_handle& ch)
    : Base(ch) {}

  // Destructor
  ~Cell_criterion_visitor() {}

  void visit(const Criterion& criterion)
  {
    Base::do_visit(criterion);
  }

};  // end class Cell_criterion_visitor


}  // end namespace Mesh_3

}  // end namespace CGAL


#endif // CGAL_MESH_3_MESH_STANDARD_CELL_CRITERIA_H
