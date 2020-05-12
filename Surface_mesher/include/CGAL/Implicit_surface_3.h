// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_IMPLICIT_SURFACE_3_H
#define CGAL_IMPLICIT_SURFACE_3_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesher/Implicit_surface_oracle_3.h>

#include <functional>
#include <functional>

namespace CGAL {

  template<
    typename GT,
    typename Function_ = std::function<typename GT::FT(typename GT::Point_3)>
    // The type of the argument `Function` will be ignored anyway.
    // The parameter is here only for backward-compatibility.
    >
  class Implicit_surface_3
  {
  public:
    typedef GT Geom_traits;
    typedef typename Geom_traits::Sphere_3 Sphere_3;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef std::function<FT(Point)> Function;
    typedef Implicit_surface_3<Geom_traits, Function_> Self;

    Function& function() { return func; }

    typedef Surface_mesher::Implicit_surface_oracle_3<
      Geom_traits,
      Self> Surface_mesher_traits_3;

    Implicit_surface_3(Function f,
                       const Sphere_3 bounding_sphere,
                       const FT error_bound = FT(1e-3),
                       Geom_traits gt = Geom_traits())
      : func(f),
        sphere(bounding_sphere),
        gt(gt)
    {
      squared_error = error_bound * error_bound;
      squared_error = squared_error *
        gt.compute_squared_radius_3_object()(bounding_sphere);
    }

    FT operator()(Point p) const
    {
      return func(p);
    }

    const FT& squared_error_bound() const
    {
      return squared_error;
    }

    const Sphere_3& bounding_sphere() const
    {
      return sphere;
    }

    const Sphere_3& bounding_sphere_squared_radius() const
    {
      return gt.compute_squared_radius_3_object()(sphere);
    }

    template <typename Vertex_handle>
    bool vertices_not_on_same_surface_patch(const Vertex_handle& v1,
                                            const Vertex_handle& v2,
                                            const Vertex_handle& v3) const
    {
      return
        v1->point().element_index() != v2->point().element_index() ||
        v1->point().element_index() != v3->point().element_index();
    }

    const Function& function() const
    {
      return func;
    }

  private:
    Function func;
    Sphere_3 sphere;
    FT squared_error;
    Geom_traits gt;
  }; // end Implicit_surface_3


  template <typename GT, typename Function>
  Implicit_surface_3<GT, Function>
  make_implicit_surface_3(GT, Function f,
                          typename GT::Sphere_3 sphere,
                          typename GT::FT error_bound)
  {
    typedef Implicit_surface_3<GT> surface;
    return surface(f, sphere, error_bound);
  }

//   template <typename GT, typename Function>
//   struct Surface_mesh_traits_generator_3<Implicit_surface_3<GT, Function> >
//   {
//     typedef Implicit_surface_3<GT, Function> Surface_type;
//     typedef typename Surface_mesher::Implicit_surface_oracle_3<GT,
//                                                              Surface_type> Type;
//     typedef Type type; // Boost meta-programming compatibility
//   };

  // non documented class
  template <typename FT, typename Point>
  class Implicit_function_wrapper : public CGAL::cpp98::unary_function<Point, FT>
  {
    typedef FT (*Implicit_function)(FT, FT, FT);

    Implicit_function function;

  public:
    Implicit_function_wrapper(Implicit_function f) : function(f) {}

    FT operator()(Point p) const
    {
      return function(p.x(), p.y(), p.z());
    }
  };

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_IMPLICIT_SURFACE_3_H
