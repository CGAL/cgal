// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2011       GeometryFactory Sarl (France)
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POISSON_IMPLICIT_SURFACE_3_H
#define CGAL_POISSON_IMPLICIT_SURFACE_3_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>


#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesher/Poisson_implicit_surface_oracle_3.h>

#include <functional>

namespace CGAL {

  template<
    typename GT,
    typename Function_
    >
  class Poisson_implicit_surface_3 
  {
  public:
    typedef GT Geom_traits;
    typedef typename Geom_traits::Sphere_3 Sphere_3;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef Function_ Function;
    typedef Poisson_implicit_surface_3<Geom_traits, Function> Self;

    Function& function() { return func; }

    typedef Surface_mesher::Poisson_implicit_surface_oracle_3<
      Geom_traits,
      Self> Surface_mesher_traits_3;

    Poisson_implicit_surface_3(Function f,
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
  }; // end Poisson_implicit_surface_3


  template <typename GT, typename Function>
  Poisson_implicit_surface_3<GT, Function>
  make_implicit_surface_3(GT, Function f,
			  typename GT::Sphere_3 sphere,
			  typename GT::FT error_bound)
  {
    typedef Poisson_implicit_surface_3<GT, Function> surface;
    return surface(f, sphere, error_bound);
  }

//   template <typename GT, typename Function>
//   struct Surface_mesh_traits_generator_3<Poisson_implicit_surface_3<GT, Function> >
//   {
//     typedef Poisson_implicit_surface_3<GT, Function> Surface_type;
//     typedef typename Surface_mesher::Poisson_implicit_surface_oracle_3<GT,
// 							     Surface_type> Type;
//     typedef Type type; // Boost meta-programming compatibility
//   };

  // non documented class
  template <typename FT, typename Point>
  class Poisson_implicit_function_wrapper : public std::unary_function<Point, FT> 
  {
    typedef FT (*Poisson_implicit_function)(FT, FT, FT);

    Poisson_implicit_function function;

  public:
    Poisson_implicit_function_wrapper(Poisson_implicit_function f) : function(f) {}

    FT operator()(Point p) const
    {
      return function(p.x(), p.y(), p.z());
    }
  };

} // end namespace CGAL

#endif // CGAL_POISSON_IMPLICIT_SURFACE_3_H
