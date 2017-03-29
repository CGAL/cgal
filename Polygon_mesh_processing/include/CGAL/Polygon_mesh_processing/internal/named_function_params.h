// Copyright (c) 2015 GeometryFactory (France).
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

#include <CGAL/license/Polygon_mesh_processing/core.h>


#include <CGAL/boost/graph/named_function_params.h>

#define CGAL_PMP_NP_TEMPLATE_PARAMETERS T, typename Tag, typename Base
#define CGAL_PMP_NP_CLASS CGAL::pmp_bgl_named_params<T,Tag,Base>

namespace CGAL{

  enum density_control_factor_t     { density_control_factor      };
  enum use_delaunay_triangulation_t { use_delaunay_triangulation  };
  enum fairing_continuity_t         { fairing_continuity };
  enum sparse_linear_solver_t       { sparse_linear_solver };
  enum geom_traits_t                { geom_traits };
  enum number_of_iterations_t       { number_of_iterations };
  enum number_of_relaxation_steps_t { number_of_relaxation_steps };
  enum protect_constraints_t        { protect_constraints };
  enum relax_constraints_t          { relax_constraints };
  enum vertex_is_constrained_t      { vertex_is_constrained };
  enum face_patch_t                 { face_patch };
  enum random_uniform_sampling_t    { random_uniform_sampling };
  enum grid_sampling_t              { grid_sampling };
  enum monte_carlo_sampling_t       { monte_carlo_sampling };
  enum do_sample_edges_t            { do_sample_edges };
  enum do_sample_vertices_t         { do_sample_vertices };
  enum do_sample_faces_t            { do_sample_faces };
  enum number_of_points_on_faces_t  { number_of_points_on_faces };
  enum number_of_points_per_face_t  { number_of_points_per_face };
  enum grid_spacing_t               { grid_spacing };
  enum nb_points_per_area_unit_t    { nb_points_per_area_unit };
  enum number_of_points_per_edge_t  { number_of_points_per_edge };
  enum number_of_points_on_edges_t  { number_of_points_on_edges };
  enum nb_points_per_distance_unit_t{ nb_points_per_distance_unit };

  //to be documented
  enum face_normal_t                { face_normal };

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

    template <typename FaceNormalMap>
    pmp_bgl_named_params<FaceNormalMap, face_normal_t, self>
    face_normal_map(const FaceNormalMap& m) const
    {
      typedef pmp_bgl_named_params<FaceNormalMap, face_normal_t, self> Params;
      return Params(m, *this);
    }

    //overload
    template <typename PointMap>
    pmp_bgl_named_params<PointMap, boost::vertex_point_t, self>
    vertex_point_map(const PointMap& p) const
    {
      typedef pmp_bgl_named_params<PointMap, boost::vertex_point_t, self> Params;
      return Params(p, *this);
    }

    template<typename K>
    pmp_bgl_named_params<K, geom_traits_t, self>
    geom_traits(const K& k) const
    {
      typedef pmp_bgl_named_params<K, geom_traits_t, self> Params;
      return Params(k, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, number_of_iterations_t, self>
    number_of_iterations(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, number_of_iterations_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, number_of_relaxation_steps_t, self>
    number_of_relaxation_steps(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, number_of_relaxation_steps_t, self> Params;
      return Params(n, *this);
    }

    template<typename Boolean>
    pmp_bgl_named_params<Boolean, protect_constraints_t, self>
    protect_constraints(const Boolean b) const
    {
      typedef pmp_bgl_named_params<Boolean, protect_constraints_t, self> Params;
      return Params(b, *this);
    }

    template<typename Boolean>
    pmp_bgl_named_params<Boolean, relax_constraints_t, self>
    relax_constraints(const Boolean b) const
    {
      typedef pmp_bgl_named_params<Boolean, relax_constraints_t, self> Params;
      return Params(b, *this);
    }

    template <typename VertexIsConstrained>
    pmp_bgl_named_params<VertexIsConstrained, vertex_is_constrained_t, self>
    vertex_is_constrained_map(const VertexIsConstrained& vm) const
    {
      typedef pmp_bgl_named_params<VertexIsConstrained, vertex_is_constrained_t, self> Params;
      return Params(vm, *this);
    }

    template <typename FacetPatch>
    pmp_bgl_named_params<FacetPatch, face_patch_t, self>
    face_patch_map(const FacetPatch& fp) const
    {
      typedef pmp_bgl_named_params<FacetPatch, face_patch_t, self> Params;
      return Params(fp, *this);
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

    //overload
    template <typename IndexMap>
    pmp_bgl_named_params<IndexMap, boost::vertex_index_t, self>
    vertex_index_map(const IndexMap& p) const
    {
      typedef pmp_bgl_named_params<IndexMap, boost::vertex_index_t, self> Params;
      return Params(p, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, random_uniform_sampling_t, self>
    use_random_uniform_sampling(const Boolean& p) const
    {
      typedef pmp_bgl_named_params<Boolean, random_uniform_sampling_t, self> Params;
      return Params(p, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, grid_sampling_t, self>
    use_grid_sampling(const Boolean& p) const
    {
      typedef pmp_bgl_named_params<Boolean, grid_sampling_t, self> Params;
      return Params(p, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, monte_carlo_sampling_t, self>
    use_monte_carlo_sampling(const Boolean& p) const
    {
      typedef pmp_bgl_named_params<Boolean, monte_carlo_sampling_t, self> Params;
      return Params(p, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, do_sample_edges_t, self>
    sample_edges(const Boolean& p) const
    {
      typedef pmp_bgl_named_params<Boolean, do_sample_edges_t, self> Params;
      return Params(p, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, do_sample_vertices_t, self>
    sample_vertices(const Boolean& p) const
    {
      typedef pmp_bgl_named_params<Boolean, do_sample_vertices_t, self> Params;
      return Params(p, *this);
    }

    template <typename Boolean>
    pmp_bgl_named_params<Boolean, do_sample_faces_t, self>
    sample_faces(const Boolean& p) const
    {
      typedef pmp_bgl_named_params<Boolean, do_sample_faces_t, self> Params;
      return Params(p, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, number_of_points_on_faces_t, self>
    number_of_points_on_faces(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, number_of_points_on_faces_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, number_of_points_per_face_t, self>
    number_of_points_per_face(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, number_of_points_per_face_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, grid_spacing_t, self>
    grid_spacing(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, grid_spacing_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, nb_points_per_area_unit_t, self>
    number_of_points_per_area_unit(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, nb_points_per_area_unit_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, number_of_points_per_edge_t, self>
    number_of_points_per_edge(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, number_of_points_per_edge_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, number_of_points_on_edges_t, self>
    number_of_points_on_edges(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, number_of_points_on_edges_t, self> Params;
      return Params(n, *this);
    }

    template<typename NT>
    pmp_bgl_named_params<NT, nb_points_per_distance_unit_t, self>
    number_of_points_per_distance_unit(const NT& n) const
    {
      typedef pmp_bgl_named_params<NT, nb_points_per_distance_unit_t, self> Params;
      return Params(n, *this);
    }
  };

namespace Polygon_mesh_processing{

namespace parameters{

  pmp_bgl_named_params<bool, all_default_t>
  inline all_default()
  {
    typedef pmp_bgl_named_params<bool, all_default_t> Params;
    return Params();
  }

  template <typename T, typename Tag, typename Base>
  pmp_bgl_named_params<T,Tag,Base>
  inline no_parameters(pmp_bgl_named_params<T,Tag,Base>)
  {
    typedef pmp_bgl_named_params<T,Tag,Base> Params;
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

  template <typename FaceNormalMap>
  pmp_bgl_named_params<FaceNormalMap, face_normal_t>
  face_normal_map(const FaceNormalMap& m)
  {
    typedef pmp_bgl_named_params<FaceNormalMap, face_normal_t> Params;
    return Params(m);
  }

  //overload
  template <typename PointMap>
  pmp_bgl_named_params<PointMap, boost::vertex_point_t>
  vertex_point_map(const PointMap& p)
  {
    typedef pmp_bgl_named_params<PointMap, boost::vertex_point_t> Params;
    return Params(p);
  }

  template<typename K>
  pmp_bgl_named_params<K, geom_traits_t>
  geom_traits(const K& k)
  {
    typedef pmp_bgl_named_params<K, geom_traits_t> Params;
    return Params(k);
  }

  template<typename NT>
  pmp_bgl_named_params<NT, number_of_iterations_t>
    number_of_iterations(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, number_of_iterations_t> Params;
    return Params(n);
  }

  template<typename NT>
  pmp_bgl_named_params<NT, number_of_relaxation_steps_t>
  number_of_relaxation_steps(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, number_of_relaxation_steps_t> Params;
    return Params(n);
  }

  template <typename Boolean>
  pmp_bgl_named_params<Boolean, protect_constraints_t>
  protect_constraints(const Boolean b)
  {
    typedef pmp_bgl_named_params<Boolean, protect_constraints_t> Params;
    return Params(b);
  }

  template<typename Boolean>
  pmp_bgl_named_params<Boolean, relax_constraints_t>
  relax_constraints(const Boolean b)
  {
    typedef pmp_bgl_named_params<Boolean, relax_constraints_t> Params;
    return Params(b);
  }

  template <typename VertexIsConstrained>
  pmp_bgl_named_params<VertexIsConstrained, vertex_is_constrained_t>
  vertex_is_constrained_map(const VertexIsConstrained& vm)
  {
    typedef pmp_bgl_named_params<VertexIsConstrained, vertex_is_constrained_t> Params;
    return Params(vm);
  }

  template <typename FacetPatch>
  pmp_bgl_named_params<FacetPatch, face_patch_t>
  face_patch_map(const FacetPatch& fp)
  {
    typedef pmp_bgl_named_params<FacetPatch, face_patch_t> Params;
    return Params(fp);
  }

  //overload
  template <typename EdgeIsConstrained>
  pmp_bgl_named_params<EdgeIsConstrained, edge_is_constrained_t>
  edge_is_constrained_map(const EdgeIsConstrained& em)
  {
    typedef pmp_bgl_named_params<EdgeIsConstrained, edge_is_constrained_t> Params;
    return Params(em);
  }

  //overload
  template <typename EdgeIsConstrainedParams>
  pmp_bgl_named_params<EdgeIsConstrainedParams, edge_is_constrained_params_t>
  edge_is_constrained_map_params(const EdgeIsConstrainedParams& em)
  {
    typedef pmp_bgl_named_params<EdgeIsConstrainedParams,
                                 edge_is_constrained_params_t> Params;
    return Params(em);
  }

  //overload
  template <typename IndexMap>
  pmp_bgl_named_params<IndexMap, boost::face_index_t>
  face_index_map(IndexMap const& p)
  {
    typedef pmp_bgl_named_params<IndexMap, boost::face_index_t> Params;
    return Params(p);
  }

  //overload
  template <typename IndexMap>
  pmp_bgl_named_params<IndexMap, boost::vertex_index_t>
  vertex_index_map(const IndexMap& p)
  {
    typedef pmp_bgl_named_params<IndexMap, boost::vertex_index_t> Params;
    return Params(p);
  }

  //overload
  template <typename Boolean>
  pmp_bgl_named_params<Boolean, random_uniform_sampling_t>
  use_random_uniform_sampling(const Boolean& p)
  {
    typedef pmp_bgl_named_params<Boolean, random_uniform_sampling_t> Params;
    return Params(p);
  }

  //overload
  template <typename Boolean>
  pmp_bgl_named_params<Boolean, grid_sampling_t>
  use_grid_sampling(const Boolean& p)
  {
    typedef pmp_bgl_named_params<Boolean, grid_sampling_t> Params;
    return Params(p);
  }

  //overload
  template <typename Boolean>
  pmp_bgl_named_params<Boolean, monte_carlo_sampling_t>
  use_monte_carlo_sampling(const Boolean& p)
  {
    typedef pmp_bgl_named_params<Boolean, monte_carlo_sampling_t> Params;
    return Params(p);
  }

  //overload
  template <typename Boolean>
  pmp_bgl_named_params<Boolean, do_sample_edges_t>
  sample_edges(const Boolean& p)
  {
    typedef pmp_bgl_named_params<Boolean, do_sample_edges_t> Params;
    return Params(p);
  }

  //overload
  template <typename Boolean>
  pmp_bgl_named_params<Boolean, do_sample_vertices_t>
  sample_vertices(const Boolean& p)
  {
    typedef pmp_bgl_named_params<Boolean, do_sample_vertices_t> Params;
    return Params(p);
  }

  //overload
  template <typename Boolean>
  pmp_bgl_named_params<Boolean, do_sample_faces_t>
  sample_faces(const Boolean& p)
  {
    typedef pmp_bgl_named_params<Boolean, do_sample_faces_t> Params;
    return Params(p);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, number_of_points_on_faces_t>
  number_of_points_on_faces(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, number_of_points_on_faces_t> Params;
    return Params(n);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, number_of_points_per_face_t>
  number_of_points_per_face(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, number_of_points_per_face_t> Params;
    return Params(n);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, grid_spacing_t>
  grid_spacing(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, grid_spacing_t> Params;
    return Params(n);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, nb_points_per_area_unit_t>
  number_of_points_per_area_unit(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, nb_points_per_area_unit_t> Params;
    return Params(n);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, number_of_points_per_edge_t>
  number_of_points_per_edge(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, number_of_points_per_edge_t> Params;
    return Params(n);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, number_of_points_on_edges_t>
  number_of_points_on_edges(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, number_of_points_on_edges_t> Params;
    return Params(n);
  }

  //overload
  template<typename NT>
  pmp_bgl_named_params<NT, nb_points_per_distance_unit_t>
  number_of_points_per_distance_unit(const NT& n)
  {
    typedef pmp_bgl_named_params<NT, nb_points_per_distance_unit_t> Params;
    return Params(n);
  }

} //namespace parameters
} //namespace Polygon_mesh_processing

} //namespace CGAL

// partial specializations hate inheritance and we need to repeat
// those here. this is rather fragile.
namespace boost {
#if BOOST_VERSION < 105100
  template <class Tag1, class Tag2, class T1, class Base>
  inline
  typename property_value< CGAL::pmp_bgl_named_params<T1,Tag1,Base>, Tag2>::type
  get_param(const CGAL::pmp_bgl_named_params<T1,Tag1,Base>& p, Tag2 tag2)
  {
    enum { match = detail::same_property<Tag1,Tag2>::value };
    typedef typename
      boost::property_value< CGAL::pmp_bgl_named_params<T1,Tag1,Base>, Tag2>::type T2;
    T2* t2 = 0;
    typedef detail::property_value_dispatch<match> Dispatcher;
    return Dispatcher::const_get_value(p, t2, tag2);
  }
#endif

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


#endif //CGAL_PMP_BGL_NAMED_FUNCTION_PARAMS_H
