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


#ifndef CGAL_MESH_3_MESH_STANDARD_FACET_CRITERIA_H
#define CGAL_MESH_3_MESH_STANDARD_FACET_CRITERIA_H


#include <CGAL/Mesh_3/mesh_standard_criteria.h>


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

    typedef typename Tr::Geom_traits Gt;
    typedef typename Tr::Point Point_3;

    const typename Gt::Construct_triangle_3 triangle =
        Gt().construct_triangle_3_object();
    const typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    const typename Gt::Compute_squared_area_3 area =
        Gt().compute_squared_area_3_object();

    const Point_3& p1 = f.first->vertex((f.second+1)&3)->point();
    const Point_3& p2 = f.first->vertex((f.second+2)&3)->point();
    const Point_3& p3 = f.first->vertex((f.second+3)&3)->point();

    const FT triangle_area = area(triangle(p1,p2,p3));

    const FT d12 = distance(p1,p2);
    const FT d13 = distance(p1,p3);
    const FT d23 = distance(p2,p3);
    const FT min_d123 = details::min_3<Gt>(d12,d13,d23);

    const FT aspect_ratio = 4 * triangle_area * min_d123 / (d12*d13*d23);

    CGAL_assertion (aspect_ratio >= 0 && aspect_ratio <= 1);

    if ( aspect_ratio < B_ )
      return Badness(Quality(aspect_ratio));
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
{
private:
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Curvature_size_criterion<Tr,Visitor_> Self;

public:
  // Nb: the default bound of the criterion is such that the criterion
  // is always fulfilled
  Curvature_size_criterion(const FT b = 0) : B_(b * b) {}

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

    typedef typename Tr::Geom_traits Gt;
    typedef typename Tr::Point Point_3;

    typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    typename Gt::Construct_circumcenter_3 circumcenter =
        Gt().construct_circumcenter_3_object();

    const Point_3& p1 = f.first->vertex((f.second+1)&3)->point();
    const Point_3& p2 = f.first->vertex((f.second+2)&3)->point();
    const Point_3& p3 = f.first->vertex((f.second+3)&3)->point();

    const Point_3 c = circumcenter(p1,p2,p3);

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

};  // end Curvature_size_criterion



// Uniform size Criterion class
template <typename Tr, typename Visitor_>
class Uniform_size_criterion :
  public Mesh_3::Abstract_criterion<Tr, Visitor_>
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
  Uniform_size_criterion(const FT b = 1000) : B_(b * b) {}

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

    typedef typename Tr::Geom_traits Gt;
    typedef typename Tr::Point Point_3;

    typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();

    const Point_3& p1 = f.first->vertex((f.second+1)&3)->point();

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
    typedef typename Tr::Vertex::Index Index;

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
      return Badness(Quality(1));
    }
    else
      return Badness();
  }
}; // end class Facet_on_surface_criterion


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

}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_STANDARD_FACET_CRITERIA_H
