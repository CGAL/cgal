// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_STANDARD_CELL_CRITERIA_H
#define CGAL_MESH_3_MESH_STANDARD_CELL_CRITERIA_H

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

    Radius sq_radius = Geom_traits().compute_squared_radius_3_object();
    Distance distance = Geom_traits().compute_squared_distance_3_object();

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
      return Badness(Quality( (sq_radius_edge_bound_*min_sq_length)/size ));
    }
    else
      return Badness();
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
  typedef typename Base::Badness Badness;

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

  virtual Badness do_is_bad(const Cell_handle& ch) const
  {
    typedef typename Tr::Point Point_3;
    typedef typename Tr::Geom_traits Geom_traits;
    typedef typename Geom_traits::Compute_squared_radius_3 Radius;

    const Point_3& p = ch->vertex(0)->point();
    const Point_3& q = ch->vertex(1)->point();
    const Point_3& r = ch->vertex(2)->point();
    const Point_3& s = ch->vertex(3)->point();

    Radius sq_radius = Geom_traits().compute_squared_radius_3_object();

    const FT size = sq_radius(p, q, r, s);

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
  typedef typename Base::Badness Badness;
  
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
  
  virtual Badness do_is_bad(const Cell_handle& ch) const
  {
    typedef typename Tr::Point Point_3;
    typedef typename Tr::Geom_traits Geom_traits;
    typedef typename Geom_traits::Compute_squared_radius_3 Radius;
    
    const Point_3& p = ch->vertex(0)->point();
    const Point_3& q = ch->vertex(1)->point();
    const Point_3& r = ch->vertex(2)->point();
    const Point_3& s = ch->vertex(3)->point();
    
    Radius sq_radius = Geom_traits().compute_squared_radius_3_object();
    
    const FT size = sq_radius(p, q, r, s);
    const FT sq_bound = CGAL::square( size_(ch->weighted_circumcenter(),
                                            3,
                                            Index(ch->subdomain_index())) );
    
    if ( size > sq_bound )
    {
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
      std::cerr << "bad cell (radius bound): size[" << size
      << "] bound[" << sq_bound << "]\n" ;
#endif
      return Badness(Quality(sq_bound/size));
    }
    else
      return Badness();
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
  typedef typename Base::Badness Badness;

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

  virtual Badness do_is_bad(const Cell_handle& ch) const
  {
    typedef typename Tr::Vertex_handle Vertex_handle;

    const Vertex_handle vp = ch->vertex(0);
    const Vertex_handle vq = ch->vertex(1);
    const Vertex_handle vr = ch->vertex(2);
    const Vertex_handle vs = ch->vertex(3);

    if(vp->in_dimension() != 2 || 
       vq->in_dimension() != 2 || 
       vr->in_dimension() != 2 || 
       vs->in_dimension() != 2) return Badness();

    typedef typename C3t3::Index Index;
    const Index& vp_index = vp->index();

    if(vq->index() != vp_index ||
       vr->index() != vp_index ||
       vs->index() != vp_index) return Badness(Quality(1));

    return Badness();
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


template <typename Tr>
class Cell_criteria_visitor_with_features
  : public Criterion_visitor<Tr, typename Tr::Cell_handle>
{
  typedef Criterion_visitor<Tr, typename Tr::Cell_handle> Base;
  typedef Cell_criteria_visitor_with_features<Tr> Self;
  

  typedef Abstract_criterion<Tr, Self>                  Criterion;
  typedef Mesh_3::Cell_size_criterion<Tr, Self>         Cell_size_criterion;
  typedef Mesh_3::Cell_radius_edge_criterion<Tr, Self>  Cell_radius_edge_criterion;

  typedef typename Tr::Geom_traits  Gt;
  typedef typename Gt::FT           FT;
  typedef typename Tr::Point        Point_3;
  
  
public:  
  typedef typename Base::Quality  Cell_quality;
  typedef typename Base::Badness  Cell_badness;
  typedef typename Base::Handle   Handle;
  typedef Handle                  Cell_handle;
  
  // Constructor
  Cell_criteria_visitor_with_features(const Cell_handle& ch)
    : Base(ch)
    , wp_nb_(0)
    , do_spheres_intersect_(false)
    , ratio_(0)
    , size_ratio_(0.5*0.5*4.)
  {
    typename Gt::Compute_squared_radius_smallest_orthogonal_sphere_3 sq_radius =
      Gt().compute_squared_radius_smallest_orthogonal_sphere_3_object();
    
    typename Gt::Compare_weighted_squared_radius_3 compare =
      Gt().compare_weighted_squared_radius_3_object();
    
    int k1 = 0;
    int k2 = 1;
    int k3 = 2;
    int k4 = 3;
    
    // Get number of weighted points, and ensure that they will be accessible
    // using k1...ki, if i is the number of weighted points.
    if(ch->vertex(k1)->point().weight() > FT(0))
    { 
      ++wp_nb_;
    }
    
    if(ch->vertex(k2)->point().weight() > FT(0))
    { 
      if ( 0 == wp_nb_ ) { std::swap(k1,k2); }
      ++wp_nb_;
    }
    
    if(ch->vertex(k3)->point().weight() > FT(0))
    { 
      if ( 0 == wp_nb_ ) { std::swap(k1,k3); }
      if ( 1 == wp_nb_ ) { std::swap(k2,k3); }
      ++wp_nb_;
    }
    
    if(ch->vertex(k4)->point().weight() > FT(0))
    { 
      if ( 0 == wp_nb_ ) { std::swap(k1,k4); }
      if ( 1 == wp_nb_ ) { std::swap(k2,k4); }
      if ( 2 == wp_nb_ ) { std::swap(k3,k4); }
      ++wp_nb_;
    }
    
    const Point_3& p1 = ch->vertex(k1)->point();
    const Point_3& p2 = ch->vertex(k2)->point();
    const Point_3& p3 = ch->vertex(k3)->point();
    const Point_3& p4 = ch->vertex(k4)->point();
    
    switch ( wp_nb_ )
    {
      case 1:
      {
        FT r12 = sq_radius(p1,p2);
        FT r13 = sq_radius(p1,p3);
        FT r14 = sq_radius(p1,p4);
        FT r = (std::max)((std::max)(r12,r13),r14);
        ratio_ = r / p1.weight();
        break;
      }
        
      case 2:
      {
        FT r13 = sq_radius(p1,p3);
        FT r14 = sq_radius(p1,p4);
        FT r1 = (std::max)(r13,r14) / p1.weight();
        FT r23 = sq_radius(p2,p3);
        FT r24 = sq_radius(p2,p4);
        FT r2 = (std::max)(r23,r24) / p2.weight();
        ratio_ = (std::max)(r1,r2);
        
        do_spheres_intersect_ = (compare(p1,p2,FT(0)) != CGAL::LARGER);
        break;
      }
        
      case 3:
      {
        FT r14 = sq_radius(p1,p4) / p1.weight();
        FT r24 = sq_radius(p2,p4) / p2.weight();
        FT r34 = sq_radius(p3,p4) / p3.weight();
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

  template <typename C3t3>
  void visit(const No_bridge_cell_criterion<C3t3, Self>& criterion)
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
