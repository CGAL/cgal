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

  enum density_control_factor_t     { density_control_factor      };
  enum use_delaunay_triangulation_t { use_delaunay_triangulation  };
  enum fairing_continuity_t         { fairing_continuity };
  enum sparse_linear_solver_t       { sparse_linear_solver };
  enum vertex_point_map_t           { vertex_point_map };
  enum less_halfedge_t              { less_halfedge };
  enum geom_traits_t                { geom_traits };

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

    //overload
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

    template<typename K>
    pmp_bgl_named_params<K, geom_traits_t, self>
    kernel(const K& k) const
    {
      typedef pmp_bgl_named_params<K, geom_traits_t, self> Params;
      return Params(k, *this);
    }

    //overload
    template <typename EdgeIsConstrained>
    pmp_bgl_named_params<EdgeIsConstrained, edge_is_constrained_t, self>
    edge_is_constrained_map(const EdgeIsConstrained& em) const
    {
      typedef pmp_bgl_named_params<EdgeIsConstrained, edge_is_constrained_t, self> Params;
      return Params(em, *this);
    }

    template <typename EdgeIsConstrainedParams>
    pmp_bgl_named_params<EdgeIsConstrainedParams, edge_is_constrained_params_t, self>
    edge_is_constrained_map_params(const EdgeIsConstrainedParams& em) const
    {
      typedef pmp_bgl_named_params<EdgeIsConstrainedParams,
                                   edge_is_constrained_params_t, self> Params;
      return Params(em, *this);
    }

    //overload
    template <typename IndexMap>
    pmp_bgl_named_params<IndexMap, boost::face_index_t, self>
      face_index_map(const IndexMap& p) const
    {
      typedef pmp_bgl_named_params<IndexMap, boost::face_index_t, self> Params;
      return Params(p, *this);
    }

  };

namespace Polygon_mesh_processing{

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

  template<typename K>
  pmp_bgl_named_params<K, geom_traits_t>
  kernel(const K& k)
  {
    typedef pmp_bgl_named_params<K, geom_traits_t> Params;
    return Params(k);
  }

  //overload
  template <typename EdgeIsConstrained>
  pmp_bgl_named_params<EdgeIsConstrained, edge_is_constrained_t>
  edge_is_constrained_map(const EdgeIsConstrained& em)
  {
    typedef pmp_bgl_named_params<EdgeIsConstrained, edge_is_constrained_t> Params;
    return Params(em);
  }

  template <typename EdgeIsConstrainedParams>
  pmp_bgl_named_params<EdgeIsConstrainedParams, edge_is_constrained_params_t>
  edge_is_constrained_map_params(const EdgeIsConstrainedParams& em)
  {
    typedef pmp_bgl_named_params<EdgeIsConstrainedParams,
                                 edge_is_constrained_params_t> Params;
    return Params(em);
  }

  template <typename IndexMap>
  pmp_bgl_named_params<IndexMap, boost::face_index_t>
  face_index_map(IndexMap const& p)
  {
    typedef pmp_bgl_named_params<IndexMap, boost::face_index_t> Params;
    return Params(p);
  }



} //namespace parameters
} //namespace Polygon_mesh_processing
} //namespace CGAL

#if BOOST_VERSION >= 105100
// partial specializations hate inheritance and we need to repeat
// those here. this is rather fragile.
namespace boost {
  template <typename T, typename Tag, typename Base, typename Def>
  struct lookup_named_param_def<Tag, CGAL::pmp_bgl_named_params<T, Tag, Base>, Def> {
    typedef T type;
    static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def&) {
      return p.m_value;
    }
  };

  template <typename Tag1, typename T, typename Tag, typename Base, typename Def>
  struct lookup_named_param_def<Tag1, CGAL::pmp_bgl_named_params<T, Tag, Base>, Def> {
    typedef typename lookup_named_param_def<Tag1, Base, Def>::type type;
    static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def& def) {
      return lookup_named_param_def<Tag1, Base, Def>::get(p.m_base, def);
    }
  };
} // boost
#endif


#endif //CGAL_PMP_BGL_NAMED_FUNCTION_PARAMS_H
