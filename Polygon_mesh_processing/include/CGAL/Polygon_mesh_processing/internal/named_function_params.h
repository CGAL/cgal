// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_PMP_BGL_NAMED_FUNCTION_PARAMS_H
#define CGAL_PMP_BGL_NAMED_FUNCTION_PARAMS_H

#include <CGAL/boost/graph/named_function_params.h>

namespace CGAL{

//namespace Polygon_mesh_processing{

  enum density_control_factor_t     { density_control_factor      };
  enum use_delaunay_triangulation_t { use_delaunay_triangulation  };
  enum fairing_continuity_t         { fairing_continuity };
  enum sparse_linear_solver_t       { sparse_linear_solver };
  enum vertex_point_map_t           { vertex_point_map };
  enum less_halfedge_t              { less_halfedge };

  //internal
  enum weight_calculator_t          { weight_calculator };
  enum all_default_t                { all_default };

  template <typename T, typename Tag, typename Base = boost::no_property>
  struct pmp_bgl_named_params
    : CGAL::cgal_bgl_named_params<T, Tag, Base>
  {
    typedef CGAL::cgal_bgl_named_params<T, Tag, Base> base;
    typedef pmp_bgl_named_params self;

    pmp_bgl_named_params(T v = T()) : base(v) {}
    pmp_bgl_named_params(T v, const Base& b) : base(v, b) {}

    pmp_bgl_named_params<bool, all_default_t, self>
    all_default() const
    {
      typedef pmp_bgl_named_params<bool, all_default_t, self> Params;
      return Params(*this);
    }

    template <typename Double>
    pmp_bgl_named_params<Double, density_control_factor_t, self>
    density_control_factor(const Double& d) const
    {
      typedef pmp_bgl_named_params<Double, density_control_factor_t, self> Params;
      return Params(d, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, use_delaunay_triangulation_t, self>
    use_delaunay_triangulation(const Boolean b) const
    {
      typedef pmp_bgl_named_params<Boolean, use_delaunay_triangulation_t, self> Params;
      return Params(b, *this);
    }

    template <typename UnsignedInt>
    pmp_bgl_named_params<UnsignedInt, fairing_continuity_t, self>
    fairing_continuity(const UnsignedInt& ui) const
    {
      typedef pmp_bgl_named_params<UnsignedInt, fairing_continuity_t, self> Params;
      return Params(ui, *this);
    }

    template<typename Solver>
    pmp_bgl_named_params<Solver, sparse_linear_solver_t, self>
    sparse_linear_solver(const Solver& s) const
    {
      typedef pmp_bgl_named_params<Solver, sparse_linear_solver_t, self> Params;
      return Params(s, *this);
    }

    template<typename WeightCalc>
    pmp_bgl_named_params<WeightCalc, weight_calculator_t, self>
    weight_calculator(const WeightCalc& w) const
    {
      typedef pmp_bgl_named_params<WeightCalc, weight_calculator_t, self> Params;
      return Params(w, *this);
    }

    template<typename VPMap>
    pmp_bgl_named_params<VPMap, vertex_point_map_t, self>
    vertex_point_map(const VPMap& vpmap) const
    {
      typedef pmp_bgl_named_params<VPMap, vertex_point_map_t, self> Params;
      return Params(vpmap, *this);
    }

    template<typename Less>
    pmp_bgl_named_params<Less, less_halfedge_t, self>
    less_halfedge(const Less& less) const
    {
      typedef pmp_bgl_named_params<Less, less_halfedge_t, self> Params;
      return Params(less, *this);
    }

  };


namespace parameters{

  pmp_bgl_named_params<bool, all_default_t>
  all_default()
  {
    typedef pmp_bgl_named_params<bool, all_default_t> Params;
    return Params();
  }

  template <typename Double>
  pmp_bgl_named_params<Double, density_control_factor_t>
  density_control_factor(const Double& d)
  {
    typedef pmp_bgl_named_params<Double, density_control_factor_t> Params;
    return Params(d);
  }

  template <typename Boolean>
  pmp_bgl_named_params<Boolean, use_delaunay_triangulation_t>
  use_delaunay_triangulation(const Boolean b)
  {
    typedef pmp_bgl_named_params<Boolean, use_delaunay_triangulation_t> Params;
    return Params(b);
  }

  template <typename UnsignedInt>
  pmp_bgl_named_params<UnsignedInt, fairing_continuity_t>
  fairing_continuity(const UnsignedInt& ui)
  {
    typedef pmp_bgl_named_params<UnsignedInt, fairing_continuity_t> Params;
    return Params(ui);
  }

  template<typename Solver>
  pmp_bgl_named_params<Solver, sparse_linear_solver_t>
  sparse_linear_solver(const Solver& s)
  {
    typedef pmp_bgl_named_params<Solver, sparse_linear_solver_t> Params;
    return Params(s);
  }

  template<typename WeightCalc>
  pmp_bgl_named_params<WeightCalc, weight_calculator_t>
  weight_calculator(const WeightCalc& w)
  {
    typedef pmp_bgl_named_params<WeightCalc, weight_calculator_t> Params;
    return Params(w);
  }

  template<typename VPMap>
  pmp_bgl_named_params<VPMap, vertex_point_map_t>
  vertex_point_map(const VPMap& vpmap)
  {
    typedef pmp_bgl_named_params<VPMap, vertex_point_map_t> Params;
    return Params(vpmap);
  }

  template<typename Less>
  pmp_bgl_named_params<Less, less_halfedge_t>
  less_halfedge(const Less& less)
  {
    typedef pmp_bgl_named_params<Less, less_halfedge_t> Params;
    return Params(less);
  }


} //namespace parameters
//} //namespace Polygon_mesh_processing
} //namespace CGAL

#endif //CGAL_PMP_BGL_NAMED_FUNCTION_PARAMS_H
