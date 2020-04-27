// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Mesh_3/mesh_standard_criteria.h>
#include <CGAL/utils.h> // for CGAL::min
#include <CGAL/number_utils.h> // for CGAL::square
#include <CGAL/enum.h>


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
  typedef typename Base::Is_bad  Is_bad;

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

  virtual Is_bad do_is_bad(const Tr& tr, const Cell_handle& ch) const
  {
    typedef typename Tr::Geom_traits    Geom_traits;
    typedef typename Tr::Bare_point     Bare_point;
    typedef typename Tr::Weighted_point Weighted_point;

    typedef typename Geom_traits::Compute_squared_distance_3 Distance;
    typedef typename Geom_traits::Compute_squared_radius_3   Radius;
    typedef typename Geom_traits::Construct_point_3          Construct_point_3;

    Distance distance = tr.geom_traits().compute_squared_distance_3_object();
    Radius sq_radius = tr.geom_traits().compute_squared_radius_3_object();
    Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    const Weighted_point& wp = tr.point(ch, 0);
    const Weighted_point& wq = tr.point(ch, 1);
    const Weighted_point& wr = tr.point(ch, 2);
    const Weighted_point& ws = tr.point(ch, 3);
    const Bare_point& p = cp(wp);
    const Bare_point& q = cp(wq);
    const Bare_point& r = cp(wr);
    const Bare_point& s = cp(ws);

    const FT size = sq_radius(p, q, r, s);

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
      return Is_bad(Quality( (sq_radius_edge_bound_*min_sq_length)/size ));
    }
    else
      return Is_bad();
  }

private:
  FT sq_radius_edge_bound_;

};  // end class Cell_radius_edge_criterion



template <typename Tr, typename Visitor_>
class Cell_size_criterion
  : public Abstract_criterion<Tr, Visitor_>
{
};

template <typename Tr, typename Visitor_>
class Cell_uniform_size_criterion
  : public Cell_size_criterion<Tr, Visitor_>
{
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Abstract_criterion<Tr, Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Is_bad  Is_bad;

  typedef Cell_uniform_size_criterion<Tr, Visitor_> Self;

public:
  // Constructor
  Cell_uniform_size_criterion(const FT& radius_bound)
    : sq_radius_bound_(radius_bound*radius_bound)   {}

  // Destructor
  ~Cell_uniform_size_criterion() {}

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

  virtual Is_bad do_is_bad(const Tr& tr, const Cell_handle& ch) const
  {
    typedef typename Tr::Geom_traits     Geom_traits;
    typedef typename Tr::Bare_point      Bare_point;
    typedef typename Tr::Weighted_point  Weighted_point;

    typedef typename Geom_traits::Compute_squared_radius_3 Radius;
    typedef typename Geom_traits::Construct_point_3        Construct_point_3;
    Radius sq_radius = tr.geom_traits().compute_squared_radius_3_object();
    Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    const Weighted_point& wp = tr.point(ch, 0);
    const Weighted_point& wq = tr.point(ch, 1);
    const Weighted_point& wr = tr.point(ch, 2);
    const Weighted_point& ws = tr.point(ch, 3);
    const Bare_point& p = cp(wp);
    const Bare_point& q = cp(wq);
    const Bare_point& r = cp(wr);
    const Bare_point& s = cp(ws);

    const FT size = sq_radius(p, q, r, s);

    if ( size > sq_radius_bound_ )
    {
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
      std::cerr << "bad cell (radius bound): size[" << size
                << "] bound[" << sq_radius_bound_ << "]\n" ;
#endif
      return Is_bad(Quality(sq_radius_bound_/size));
    }
    else
      return Is_bad();
  }

private:
  FT sq_radius_bound_;

};  // end class Cell_uniform_size_criterion


template <typename Tr, typename Visitor_, typename SizingField>
class Cell_variable_size_criterion
: public Cell_size_criterion<Tr, Visitor_>
{
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Geom_traits::FT  FT;
  typedef typename Tr::Vertex::Index    Index;

  typedef Abstract_criterion<Tr, Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Is_bad  Is_bad;

  typedef Cell_variable_size_criterion<Tr, Visitor_, SizingField> Self;
  typedef SizingField Sizing_field;

public:
  // Constructor
  Cell_variable_size_criterion(const Sizing_field& s)
    : size_(s)   {}

  // Destructor
  ~Cell_variable_size_criterion() {}

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

  virtual Is_bad do_is_bad(const Tr& tr, const Cell_handle& ch) const
  {
    typedef typename Tr::Geom_traits      Geom_traits;
    typedef typename Tr::Bare_point       Bare_point;
    typedef typename Tr::Weighted_point   Weighted_point;

    typedef typename Geom_traits::Compute_squared_radius_3 Radius;
    typedef typename Geom_traits::Construct_point_3        Construct_point_3;
    Radius sq_radius = tr.geom_traits().compute_squared_radius_3_object();
    Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    const Weighted_point& wp = tr.point(ch, 0);
    const Weighted_point& wq = tr.point(ch, 1);
    const Weighted_point& wr = tr.point(ch, 2);
    const Weighted_point& ws = tr.point(ch, 3);
    const Bare_point& p = cp(wp);
    const Bare_point& q = cp(wq);
    const Bare_point& r = cp(wr);
    const Bare_point& s = cp(ws);

    const FT size = sq_radius(p, q, r, s);
    const FT sq_bound = CGAL::square( size_(tr.dual(ch), 3,
                                            Index(ch->subdomain_index())) );

    if ( size > sq_bound )
    {
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
      std::cerr << "bad cell (radius bound): size[" << size
      << "] bound[" << sq_bound << "]\n" ;
#endif
      return Is_bad(Quality(sq_bound/size));
    }
    else
      return Is_bad();
  }

private:
  Sizing_field size_;

};  // end class Cell_variable_size_criterion

/// New cell criterion that disallows a cell to have points on different
/// surfaces, if they are all of dimension 2.
template <typename C3t3_, typename Visitor_>
class No_bridge_cell_criterion
  : public Abstract_criterion<typename C3t3_::Triangulation, Visitor_>
{
  typedef C3t3_ C3t3;
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Is_bad  Is_bad;

  typedef No_bridge_cell_criterion<C3t3, Visitor_> Self;

public:
  // Constructor
  No_bridge_cell_criterion()
  {}

  // Destructor
  ~No_bridge_cell_criterion() {}


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

  virtual Is_bad do_is_bad(const Tr& /*tr*/, const Cell_handle& ch) const
  {
    typedef typename Tr::Vertex_handle Vertex_handle;

    const Vertex_handle vp = ch->vertex(0);
    const Vertex_handle vq = ch->vertex(1);
    const Vertex_handle vr = ch->vertex(2);
    const Vertex_handle vs = ch->vertex(3);

    if(vp->in_dimension() != 2 ||
       vq->in_dimension() != 2 ||
       vr->in_dimension() != 2 ||
       vs->in_dimension() != 2)
      return Is_bad();

    typedef typename C3t3::Index Index;
    const Index& vp_index = vp->index();

    if(vq->index() != vp_index ||
       vr->index() != vp_index ||
       vs->index() != vp_index)
      return Is_bad(Quality(1));

    return Is_bad();
  }
};  // end class No_bridge_cell_criterion


template <typename Tr>
class Cell_criterion_visitor
  : public Criterion_visitor<Tr, typename Tr::Cell_handle>
{
  typedef Criterion_visitor<Tr, typename Tr::Cell_handle> Base;
  typedef Cell_criterion_visitor<Tr> Self;

public:
  typedef Abstract_criterion<Tr, Self> Criterion;
  typedef typename Base::Quality Cell_quality;
  typedef typename Base::Is_bad  Is_cell_bad;
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


template <typename Tr>
class Cell_criteria_visitor_with_features
  : public Criterion_visitor<Tr, typename Tr::Cell_handle>
{
  typedef Criterion_visitor<Tr, typename Tr::Cell_handle> Base;
  typedef Cell_criteria_visitor_with_features<Tr> Self;


  typedef Abstract_criterion<Tr, Self>                  Criterion;
  typedef Mesh_3::Cell_size_criterion<Tr, Self>         Cell_size_criterion;
  typedef Mesh_3::Cell_radius_edge_criterion<Tr, Self>  Cell_radius_edge_criterion;

  typedef typename Tr::Geom_traits    Gt;
  typedef typename Gt::FT             FT;
  typedef typename Tr::Weighted_point Weighted_point;


public:
  typedef typename Base::Quality  Cell_quality;
  typedef typename Base::Is_bad   Is_cell_bad;
  typedef typename Base::Handle   Handle;
  typedef Handle                  Cell_handle;

  // Constructor
  Cell_criteria_visitor_with_features(const Tr& tr, const Cell_handle& ch)
    : Base(tr, ch)
    , wp_nb_(0)
    , do_spheres_intersect_(false)
    , ratio_(0)
    , size_ratio_(0.5*0.5*4.)
  {
    typename Gt::Compare_weighted_squared_radius_3 compare =
      tr.geom_traits().compare_weighted_squared_radius_3_object();
    typename Gt::Compute_weight_3 cw =
      tr.geom_traits().compute_weight_3_object();
    typename Gt::Compute_squared_radius_smallest_orthogonal_sphere_3 sq_radius =
      tr.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    int k1 = 0;
    int k2 = 1;
    int k3 = 2;
    int k4 = 3;

    // Get number of weighted points, and ensure that they will be accessible
    // using k1...ki, if i is the number of weighted points.
    const Weighted_point& wpk1 = tr.point(ch, k1);
    if(compare(wpk1, FT(0)) == CGAL::SMALLER) // 0 < wpk1's weight
    {
      ++wp_nb_;
    }

    const Weighted_point& wpk2 = tr.point(ch, k2);
    if(compare(wpk2, FT(0)) == CGAL::SMALLER)
    {
      if ( 0 == wp_nb_ ) { std::swap(k1,k2); }
      ++wp_nb_;
    }

    const Weighted_point& wpk3 = tr.point(ch, k3);
    if(compare(wpk3, FT(0)) == CGAL::SMALLER)
    {
      if ( 0 == wp_nb_ ) { std::swap(k1,k3); }
      if ( 1 == wp_nb_ ) { std::swap(k2,k3); }
      ++wp_nb_;
    }

    const Weighted_point& wpk4 = tr.point(ch, k4);
    if(compare(wpk4, FT(0)) == CGAL::SMALLER)
    {
      if ( 0 == wp_nb_ ) { std::swap(k1,k4); }
      if ( 1 == wp_nb_ ) { std::swap(k2,k4); }
      if ( 2 == wp_nb_ ) { std::swap(k3,k4); }
      ++wp_nb_;
    }

    const Weighted_point& p1 = tr.point(ch, k1);
    const Weighted_point& p2 = tr.point(ch, k2);
    const Weighted_point& p3 = tr.point(ch, k3);
    const Weighted_point& p4 = tr.point(ch, k4);

    switch ( wp_nb_ )
    {
      case 1:
      {
        FT r12 = sq_radius(p1,p2);
        FT r13 = sq_radius(p1,p3);
        FT r14 = sq_radius(p1,p4);
        FT r = (std::max)((std::max)(r12,r13),r14);
        ratio_ = r / cw(p1);
        break;
      }

      case 2:
      {
        FT r13 = sq_radius(p1,p3);
        FT r14 = sq_radius(p1,p4);
        FT r1 = (std::max)(r13,r14) / cw(p1);
        FT r23 = sq_radius(p2,p3);
        FT r24 = sq_radius(p2,p4);
        FT r2 = (std::max)(r23,r24) / cw(p2);
        ratio_ = (std::max)(r1,r2);

        do_spheres_intersect_ = (compare(p1,p2,FT(0)) != CGAL::LARGER);
        break;
      }

      case 3:
      {
        FT r14 = sq_radius(p1,p4) / cw(p1);
        FT r24 = sq_radius(p2,p4) / cw(p2);
        FT r34 = sq_radius(p3,p4) / cw(p3);
        ratio_ = (std::max)((std::max)(r14,r24),r34);

        do_spheres_intersect_ = (compare(p1,p2,p3,FT(0)) != CGAL::LARGER);
        break;
      }

      case 4:
        do_spheres_intersect_ = (compare(p1,p2,p3,p4,FT(0)) != CGAL::LARGER);
        break;

      default:
        break;
    }
  }

  // Destructor
  ~Cell_criteria_visitor_with_features() {}

  // visit functions
  void visit(const Cell_size_criterion& criterion)
  {
    if (   ratio_ < size_ratio_
        && (do_spheres_intersect_ || 1 == wp_nb_) )
    {
      Base::increment_counter();
      return;
    }

    Base::do_visit(criterion);
  }

  void visit(const Cell_radius_edge_criterion& criterion)
  {
    if (   (wp_nb_ >= 2 && do_spheres_intersect_)
        || 1 == wp_nb_ )
    {
      Base::increment_counter();
      return;
    }

    Base::do_visit(criterion);
  }

  void visit(const Criterion& criterion)
  {
    Base::do_visit(criterion);
  }

private:
  int wp_nb_;
  bool do_spheres_intersect_;
  FT ratio_;
  FT size_ratio_;
};  // end class Cell_criterion_visitor



}  // end namespace Mesh_3

}  // end namespace CGAL


#endif // CGAL_MESH_3_MESH_STANDARD_CELL_CRITERIA_H
