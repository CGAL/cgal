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


#ifndef CGAL_MESH_3_MESH_STANDARD_FACET_CRITERIA_H
#define CGAL_MESH_3_MESH_STANDARD_FACET_CRITERIA_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/number_utils.h> // for to_double
#include <CGAL/Mesh_3/mesh_standard_criteria.h>
#include <cmath>


namespace CGAL {

namespace Mesh_3 {

namespace details {

  template<typename K>
  inline
  typename K::FT
  min_3(const typename K::FT& a,
        const typename K::FT& b,
        const typename K::FT& c)
  {
    return (std::min)(a, (std::min)(b,c));
  }

} // end namespace details



// Aspect_ratio Criterion class
template <typename Tr, typename Visitor_>
class Aspect_ratio_criterion :
  public Mesh_3::Abstract_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Aspect_ratio_criterion<Tr,Visitor_> Self;

public:
  // Nb: the default bound of the criterion is such that the criterion
  // is always fulfilled
  Aspect_ratio_criterion(const FT angle_min = 0.)
  { // TODO: document that FT must constructible from a double!
    set_angle_min(angle_min);
  }

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

  void set_angle_min(const FT angle_min)
  {
    if(angle_min == FT(0))
    {
      B_ = 0;
    }
    else
    {
      B_ = std::sin (CGAL_PI * CGAL::to_double(angle_min) / 180);
      B_ = B_ * B_;
    }
  }

  virtual Badness do_is_bad (const Facet& f) const
  {
    CGAL_assertion (f.first->is_facet_on_surface(f.second));
    CGAL_assertion (B_ != 0);

    typedef typename Tr::Geom_traits    Gt;
    typedef typename Tr::Bare_point     Bare_point;

    const typename Gt::Construct_triangle_3 triangle =
        Gt().construct_triangle_3_object();
    const typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    const typename Gt::Compute_squared_area_3 area =
        Gt().compute_squared_area_3_object();
    const typename Gt::Construct_point_3 wp2p =
        Gt().construct_point_3_object();

    const Bare_point& p1 = wp2p(f.first->vertex((f.second+1)&3)->point());
    const Bare_point& p2 = wp2p(f.first->vertex((f.second+2)&3)->point());
    const Bare_point& p3 = wp2p(f.first->vertex((f.second+3)&3)->point());

    const FT triangle_area = area(triangle(p1,p2,p3));
    const FT d12 = distance(p1,p2);
    const FT d13 = distance(p1,p3);
    const FT d23 = distance(p2,p3);
    const FT min_d123 = details::min_3<Gt>(d12,d13,d23);

    const FT aspect_ratio = 4 * triangle_area * min_d123 / (d12*d13*d23);

    CGAL_assertion (aspect_ratio >= 0 && aspect_ratio <= 1);

    if ( aspect_ratio < B_ )
    {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (angle criterion): aspect_ratio[" << aspect_ratio
                << "] bound[" << B_ << "]" << std::endl;
#endif
      return Badness(Quality(aspect_ratio));
    }
    else
      return Badness();
  }

private:
  FT B_;

};  // end Aspect_ratio_criterion


// Curvature_adapted size Criterion class
template <typename Tr, typename Visitor_>
class Curvature_size_criterion :
  public Mesh_3::Abstract_criterion<Tr, Visitor_>
{};

template <typename Tr, typename Visitor_>
class Uniform_curvature_size_criterion :
  public Curvature_size_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Uniform_curvature_size_criterion<Tr,Visitor_> Self;

public:
  // Nb: the default bound of the criterion is such that the criterion
  // is always fulfilled
  Uniform_curvature_size_criterion(const FT b = 0) : B_(b * b) {}

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

  virtual Badness do_is_bad (const Facet& f) const
  {
    CGAL_assertion(f.first->is_facet_on_surface(f.second));
    CGAL_assertion (B_ != 0);

    typedef typename Tr::Geom_traits    Gt;
    typedef typename Tr::Weighted_point Weighted_point;
    typedef typename Tr::Bare_point Bare_point;

    typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    typename Gt::Construct_weighted_circumcenter_3 circumcenter =
        Gt().construct_weighted_circumcenter_3_object();

    const Weighted_point& p1 = f.first->vertex((f.second+1)&3)->point();
    const Weighted_point& p2 = f.first->vertex((f.second+2)&3)->point();
    const Weighted_point& p3 = f.first->vertex((f.second+3)&3)->point();

    const Bare_point c = circumcenter(p1,p2,p3);

    const FT sq_dist = distance(c, f.first->get_facet_surface_center(f.second));

    if ( sq_dist > B_ )
    {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (curvature size): sq_dist[" << sq_dist
                << "] bound[" << B_ << "]\n";
#endif
      return Badness(Quality(B_/sq_dist));
    }
    else
      return Badness();
  }

private:
  FT B_;

};  // end Uniform_curvature_size_criterion

// Variable size Criterion class
template <typename Tr, typename Visitor_, typename SizingField>
class Variable_curvature_size_criterion :
  public Curvature_size_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet            Facet;
  typedef typename Tr::Geom_traits::FT  FT;
  typedef typename Tr::Vertex::Index    Index;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Variable_curvature_size_criterion<Tr,Visitor_,SizingField> Self;
  typedef SizingField Sizing_field;

public:
  // Nb: the default bound of the criterion is such that the criterion
  // is always fulfilled
  Variable_curvature_size_criterion(const Sizing_field& s) : size_(s) {}

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

  virtual Badness do_is_bad (const Facet& f) const
  {
    CGAL_assertion (f.first->is_facet_on_surface(f.second));

    typedef typename Tr::Geom_traits    Gt;
    typedef typename Tr::Weighted_point Weighted_point;
    typedef typename Tr::Bare_point Bare_point;

    typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    typename Gt::Construct_weighted_circumcenter_3 circumcenter =
        Gt().construct_weighted_circumcenter_3_object();

    const Weighted_point& p1 = f.first->vertex((f.second+1)&3)->point();
    const Weighted_point& p2 = f.first->vertex((f.second+2)&3)->point();
    const Weighted_point& p3 = f.first->vertex((f.second+3)&3)->point();

    const Bare_point c = circumcenter(p1,p2,p3);
    const Bare_point& ball_center = f.first->get_facet_surface_center(f.second);

    const FT sq_dist = distance(c, ball_center);

    const Index& index = f.first->get_facet_surface_center_index(f.second);

    const FT sq_bound = CGAL::square(size_(ball_center, 2, index));
    CGAL_assertion(sq_bound > FT(0));

    if ( sq_dist > sq_bound )
    {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (curvature size): sq_dist[" << sq_dist
      << "] bound[" << sq_bound << "]\n";
#endif
      return Badness(Quality(sq_bound/sq_dist));
    }
    else
      return Badness();
  }

private:
  Sizing_field size_;

};  // end Variable_curvature_size_criterion

// Size Criterion base class
template < typename Tr, typename Visitor_ >
class Facet_size_criterion :
  public Mesh_3::Abstract_criterion<Tr, Visitor_>
{
};
  
// Variable size Criterion class
template <typename Tr, typename Visitor_, typename SizingField>
class Variable_size_criterion :
  public Facet_size_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet            Facet;
  typedef typename Tr::Geom_traits::FT  FT;
  typedef typename Tr::Vertex::Index    Index;
  
  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;
  
  typedef Variable_size_criterion<Tr,Visitor_,SizingField> Self;
  typedef SizingField Sizing_field;
  
public:
  // Nb: the default bound of the criterion is such that the criterion
  // is always fulfilled
  Variable_size_criterion(const Sizing_field& s) : size_(s) {}
  
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
  
  virtual Badness do_is_bad (const Facet& f) const
  {
    CGAL_assertion (f.first->is_facet_on_surface(f.second));
    
    typedef typename Tr::Geom_traits    Gt;
    typedef typename Tr::Bare_point Bare_point;
    
    typename Gt::Compute_squared_distance_3 distance =
      Gt().compute_squared_distance_3_object();

    typename Gt::Construct_point_3 wp2p = Gt().construct_point_3_object();
    
    const Bare_point& p1 = wp2p(f.first->vertex((f.second+1)&3)->point());
    const Bare_point& ball_center = f.first->get_facet_surface_center(f.second);
    const Index& index = f.first->get_facet_surface_center_index(f.second);
    
    const FT sq_radius = distance(p1,ball_center);
    const FT sq_bound = CGAL::square(size_(ball_center, 2, index));
    CGAL_assertion(sq_bound > FT(0));
    
    if ( sq_radius > sq_bound )
    {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (uniform size): sq_radius[" << sq_radius
      << "] bound[" << sq_bound << "]\n";
#endif
      return Badness(Quality(sq_bound/sq_radius));
    }
    else
      return Badness();
  }
  
private:
  Sizing_field size_;
  
};  // end Variable_size_criterion
  
  
  
// Uniform size Criterion class
template <typename Tr, typename Visitor_>
class Uniform_size_criterion :
  public Facet_size_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Uniform_size_criterion<Tr,Visitor_> Self;
  
public:
  // Nb: the default bound of the criterion is such that the criterion
  // is always fulfilled
  Uniform_size_criterion(const FT b = 1e20) : B_(b * b) {}

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

  virtual Badness do_is_bad (const Facet& f) const
  {
    CGAL_assertion (f.first->is_facet_on_surface(f.second));
    CGAL_assertion (B_ != 0);

    typedef typename Tr::Geom_traits    Gt;
    typedef typename Tr::Bare_point     Bare_point;

    typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    typename Gt::Construct_point_3 wp2p =
        Gt().construct_point_3_object();

    const Bare_point p1 = wp2p(f.first->vertex((f.second+1)&3)->point());

    const FT sq_radius = distance(
        p1, f.first->get_facet_surface_center(f.second));

    if ( sq_radius > B_ )
    {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (uniform size): sq_radius[" << sq_radius
                << "] bound[" << B_ << "]\n";
#endif
      return Badness(Quality(B_/sq_radius));
    }
    else
      return Badness();
  }

private:
  FT B_;

};  // end Uniform_size_criterion



template <typename Tr, typename Visitor_>
class Facet_on_surface_criterion :
  public Mesh_3::Abstract_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Facet_on_surface_criterion<Tr,Visitor_> Self;

public:
  /// Constructor
  Facet_on_surface_criterion() {}
  /// Destructor
  ~Facet_on_surface_criterion() {}

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

  virtual Badness do_is_bad (const Facet& f) const
  {
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Cell_handle Cell_handle;

    const Cell_handle& ch = f.first;
    const int i = f.second;
    const Vertex_handle& v1 = ch->vertex((i+1)&3);
    const Vertex_handle& v2 = ch->vertex((i+2)&3);
    const Vertex_handle& v3 = ch->vertex((i+3)&3);

    // Look if vertex are on surface
    if ( (v1->in_dimension() > 2) ||
         (v2->in_dimension() > 2) ||
         (v3->in_dimension() > 2) )
    {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (on surface criterion)" << std::endl;
#endif
      return Badness(Quality(1));
    }
    else
      return Badness();
  }
}; // end class Facet_on_surface_criterion

  
template <typename Tr, typename Visitor_>
class Facet_on_same_surface_criterion :
public Mesh_3::Abstract_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;
  
  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;
  
  typedef Facet_on_same_surface_criterion<Tr,Visitor_> Self;
  
public:
  /// Constructor
  Facet_on_same_surface_criterion() {}
  /// Destructor
  virtual ~Facet_on_same_surface_criterion() {}

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
  
  virtual Badness do_is_bad (const Facet& f) const
  {
    typedef typename Tr::Vertex_handle  Vertex_handle;
    typedef typename Tr::Cell_handle    Cell_handle;
    typedef typename Tr::Vertex::Index  Index;
    
    const Cell_handle& ch = f.first;
    const int& i = f.second;
    
    const Vertex_handle& v1 = ch->vertex((i+1)&3);
    const Vertex_handle& v2 = ch->vertex((i+2)&3);
    const Vertex_handle& v3 = ch->vertex((i+3)&3);
    
    Index index = Index();
    bool is_index_initialized = false;
    
    if ( v1->in_dimension() == 2 )
    { 
      index = v1->index();
      is_index_initialized = true;
    }
    
    if ( v2->in_dimension() == 2 )
    {
      if ( is_index_initialized )
      {
        if ( !(v2->index() == index) )
        {
          return Badness(Quality(1));
        }
      }
      else
      {
        index = v2->index();
        is_index_initialized = true;        
      }
    }
    
    if ( v3->in_dimension() == 2 )
    {
      if ( is_index_initialized && !(v3->index() == index) )
      {
        return Badness(Quality(1));
      } 
    }
    
    return  Badness();			
  }
  
}; // end class Facet_on_same_surface_criterion



template <typename Tr>
class Facet_criterion_visitor
  : public Mesh_3::Criterion_visitor<Tr, typename Tr::Facet>
{
  typedef Mesh_3::Criterion_visitor<Tr, typename Tr::Facet> Base;
  typedef Facet_criterion_visitor<Tr> Self;

public:
  typedef Mesh_3::Abstract_criterion<Tr, Self> Criterion;
  typedef typename Base::Quality Facet_quality;
  typedef typename Base::Badness Facet_badness;
  typedef typename Base::Handle Handle;
  typedef Handle Facet;

  // Constructor
  Facet_criterion_visitor(const Facet& f)
    : Base(f) {}

  // Destructor
  ~Facet_criterion_visitor() {}

  void visit(const Criterion& criterion)
  {
    Base::do_visit(criterion);
  }

};  // end class Facet_criterion_visitor
  
  
  
template <typename Tr>
class Facet_criterion_visitor_with_features
  : public Mesh_3::Criterion_visitor<Tr, typename Tr::Facet>
{
  typedef Mesh_3::Criterion_visitor<Tr, typename Tr::Facet> Base;
  typedef Facet_criterion_visitor_with_features<Tr> Self;
  
  typedef Mesh_3::Abstract_criterion<Tr, Self>                Criterion;
  typedef Mesh_3::Curvature_size_criterion<Tr, Self>          Curvature_size_criterion;
  typedef Mesh_3::Aspect_ratio_criterion<Tr, Self>            Aspect_ratio_criterion;
  typedef Mesh_3::Facet_on_surface_criterion<Tr, Self>        Facet_on_surface_criterion;
  typedef Mesh_3::Facet_size_criterion<Tr, Self>              Facet_size_criterion;
  typedef Mesh_3::Facet_on_same_surface_criterion<Tr, Self>   Facet_on_same_surface_criterion;

  typedef typename Tr::Geom_traits  Gt;
  typedef typename Gt::FT           FT;

public:  
  typedef typename Base::Quality  Facet_quality;
  typedef typename Base::Badness  Facet_badness;
  typedef typename Base::Handle   Handle;
  typedef Handle                  Facet;
  
  // Constructor
  Facet_criterion_visitor_with_features(const Facet& fh)
    : Base(fh)
    , wp_nb_(0)
    , do_spheres_intersect_(false)
    , ratio_(0.)
    , approx_ratio_(0.1*0.1*4.)
    , angle_ratio_(0.5*0.5*4.)
    , size_ratio_(0.4*0.4*4.)
  {
    typedef typename Tr::Geom_traits    Gt;
    typedef typename Tr::Weighted_point Weighted_point;
    typedef typename Tr::Cell_handle    Cell_handle;
    
    typename Gt::Compare_weighted_squared_radius_3 compare =
      Gt().compare_weighted_squared_radius_3_object();
    
    typename Gt::Compute_squared_radius_smallest_orthogonal_sphere_3 sq_radius =
      Gt().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    const Cell_handle& c = fh.first;
    const int& k = fh.second;
    
    int k1 = (k+1)&3;
    int k2 = (k+2)&3;
    int k3 = (k+3)&3;
    
    // Get number of weighted points, and ensure that they will be accessible
    // using k1...ki, if i is the number of weighted points.
    if(c->vertex(k1)->point().weight() > FT(0))
    { 
      ++wp_nb_;
    }
    
    if(c->vertex(k2)->point().weight() > FT(0))
    { 
      if ( 0 == wp_nb_ ) { std::swap(k1,k2); }
      ++wp_nb_;
    }
    
    if(c->vertex(k3)->point().weight() > FT(0))
    { 
      if ( 0 == wp_nb_ ) { std::swap(k1,k3); }
      if ( 1 == wp_nb_ ) { std::swap(k2,k3); }
      ++wp_nb_;
    }
    
    const Weighted_point& p1 = c->vertex(k1)->point();
    const Weighted_point& p2 = c->vertex(k2)->point();
    const Weighted_point& p3 = c->vertex(k3)->point();
    
    // Compute ratio
    switch ( wp_nb_ )
    {
      case 1:
      {
        FT r = (std::max)(sq_radius(p1,p2),sq_radius(p1,p3));
        ratio_ = r / p1.weight();
        break;
      }
        
      case 2:
      {
        FT r13 = sq_radius(p1,p3) / p1.weight();
        FT r23 = sq_radius(p2,p3) / p2.weight();
        ratio_ = (std::max)(r13, r23);
        
        do_spheres_intersect_ = (compare(p1,p2,FT(0)) != CGAL::LARGER);
        break;
      }

      case 3:
      {
        do_spheres_intersect_ = (compare(p1,p2,p3,FT(0)) != CGAL::LARGER);
        break;
      }  
      
      default: break;
    }
  }
  
  // Destructor
  ~Facet_criterion_visitor_with_features() {}

  // visit functions
  void visit(const Criterion& criterion)
  {
    if ( 3 == wp_nb_ && do_spheres_intersect_ )
    { 
      Base::increment_counter();
      return;
    }
    
    Base::do_visit(criterion);
  }
  
  void visit(const Curvature_size_criterion& criterion)
  {
    if (   ratio_ < approx_ratio_
        && (do_spheres_intersect_ || 1 == wp_nb_ ) )
    {
      Base::increment_counter();
      return;
    }
    
    Base::do_visit(criterion);
  }
  
  void visit(const Aspect_ratio_criterion& criterion)
  {
    if (   ratio_ < angle_ratio_
        && (do_spheres_intersect_ || 1 == wp_nb_) )
    {
      Base::increment_counter();
      return;
    }
    
    Base::do_visit(criterion);
  }
  
  void visit(const Facet_size_criterion& criterion)
  {
    if (   ratio_ < size_ratio_
        && (do_spheres_intersect_ || 1 == wp_nb_) )
    {
      Base::increment_counter();
      return;
    }
    
    Base::do_visit(criterion);
  }
  
private:
  int wp_nb_;
  bool do_spheres_intersect_;
  FT ratio_;
  FT approx_ratio_;
  FT angle_ratio_;
  FT size_ratio_;
  
};  // end class Facet_criterion_visitor
  

}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_STANDARD_FACET_CRITERIA_H
