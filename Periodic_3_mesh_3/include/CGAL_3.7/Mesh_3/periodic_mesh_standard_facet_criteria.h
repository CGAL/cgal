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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/Mesh_3/mesh_standard_facet_criteria.h $
// $Id: mesh_standard_facet_criteria.h 53419 2009-12-15 14:26:19Z lrineau $
//
//
// Author(s)     : Mikhail Bogdanov
//
//******************************************************************************
// File Description :
//
//******************************************************************************


#ifndef CGAL_PERIODIC_MESH_STANDARD_FACET_CRITERIA_H
#define CGAL_PERIODIC_MESH_STANDARD_FACET_CRITERIA_H

#include <CGAL/Mesh_3/mesh_standard_criteria.h>


namespace CGAL {

namespace Mesh_3 {

namespace Periodic_mesh_3 {

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
  Aspect_ratio_criterion(const Tr& tr_, const FT angle_min = 0.) : tr(tr_)
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
    typedef typename Tr::Point Point;

    const typename Gt::Construct_triangle_3 triangle =
        Gt().construct_triangle_3_object();
    const typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    const typename Gt::Compute_squared_area_3 area =
        Gt().compute_squared_area_3_object();

    CGAL_assertion(tr.tds().is_cell(
      f.first->vertex(0),f.first->vertex(1),
      f.first->vertex(2),f.first->vertex(3)));
    const Point p1 = tr.point(f.first, (f.second+1)&3);
    const Point p2 = tr.point(f.first, (f.second+2)&3);
    const Point p3 = tr.point(f.first, (f.second+3)&3);

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

  const Tr& tr;
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
  Curvature_size_criterion(const Tr& tr_, const FT b = 0) : B_(b * b), tr(tr_) {}

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
    typedef typename Tr::Point Point;

    typename Gt::Compute_squared_distance_3 distance =
        Gt().compute_squared_distance_3_object();
    typename Gt::Construct_weighted_circumcenter_3 circumcenter =
        Gt().construct_weighted_circumcenter_3_object();

    const Point& p1 = tr.point(f.first, (f.second+1)&3);
    const Point& p2 = tr.point(f.first, (f.second+2)&3);
    const Point& p3 = tr.point(f.first, (f.second+3)&3);

    const Point c = tr.canonicalize_point(circumcenter(p1,p2,p3));

    // TODO: normally facet_surface_center comes from the oracle and
    // should not need to be canonicalized!
    const FT sq_dist = distance(c, tr.canonicalize_point(f.first->get_facet_surface_center(f.second)));

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

  const Tr& tr;
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
  Uniform_size_criterion(const Tr& tr_, const FT b = 1000) : B_(b * b), tr(tr_) {}

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
    typedef typename Tr::Point Point;

    typename Gt::Compute_squared_distance_3 distance =
      Gt().compute_squared_distance_3_object();

    Point v1 = f.first->vertex ((f.second+1)&3)->point();
    Point surface_center =  f.first->get_facet_surface_center (f.second);

    Point surface_centers[27];
    typedef typename Tr::Offset Offset;
    for( int i = 0; i < 3; i++ ) {
      for( int j = 0; j < 3; j++) {
        for( int k = 0; k < 3; k++ ) {
          surface_centers[9*i+3*j+k] = tr.point(
            std::make_pair(surface_center, Offset(i-1,j-1,k-1)));
        }
      }
    }
    FT min_distance = distance(v1, surface_centers[0]);
    for( int i = 1; i < 27; i++ ) {
      FT current_distance = distance(v1, surface_centers[i]);
      if ( current_distance < min_distance )
      {
        min_distance = current_distance;
      }
    }
    const FT sq_radius = min_distance;

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

  const Tr& tr;
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
  Facet_on_surface_criterion() {};
  /// Destructor
  ~Facet_on_surface_criterion() {};

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

} // end namespace Periodic_mesh_3

}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_PERIODIC_MESH_STANDARD_FACET_CRITERIA_H
