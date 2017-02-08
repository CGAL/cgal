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

#include <CGAL/license/Polygon_mesh_processing.h>


#include <CGAL/boost/graph/named_function_params.h>

#define CGAL_PMP_NP_TEMPLATE_PARAMETERS T, typename Tag, typename Base
#define CGAL_PMP_NP_CLASS CGAL::pmp_bgl_named_params<T,Tag,Base>
namespace CGAL{
template <typename T, typename Tag, typename Base = boost::no_property>
struct pmp_bgl_named_params;
}

#define Cgal_add_pmp_parameter(X_t, X)        \
namespace CGAL{                               \
namespace Polygon_mesh_processing{            \
namespace parameters{                         \
template<typename K>                          \
pmp_bgl_named_params<K, X_t>                  \
X(const K& k)                                 \
{                                             \
  typedef pmp_bgl_named_params<K, X_t> Params;\
  return Params(k);                           \
}                                             \
}                                             \
}                                             \
}

#define Cgal_add_pmp_parameter_and_enum(X_t, X) \
  namespace CGAL{                               \
  enum X_t { X };                               \
  namespace Polygon_mesh_processing{            \
  namespace parameters{                         \
  template<typename K>                          \
  pmp_bgl_named_params<K, X_t>                  \
  X(const K& k)                                 \
  {                                             \
    typedef pmp_bgl_named_params<K, X_t> Params;\
    return Params(k);                           \
  }                                             \
  }                                             \
  }                                             \
}

#define Cgal_add_pmp_parameter_in_struct(X_t, X)       \
  template<typename K>                                 \
  pmp_bgl_named_params<K, X_t, self>                   \
  X(const K& k) const                                  \
  {                                                    \
    typedef pmp_bgl_named_params<K, X_t, self> Params; \
    return Params(k, *this);                           \
  }

Cgal_add_pmp_parameter_and_enum(geom_traits_t, geom_traits)
Cgal_add_pmp_parameter_and_enum(density_control_factor_t, density_control_factor)
Cgal_add_pmp_parameter_and_enum(use_delaunay_triangulation_t, use_delaunay_triangulation)
Cgal_add_pmp_parameter_and_enum(fairing_continuity_t, fairing_continuity)
Cgal_add_pmp_parameter_and_enum(sparse_linear_solver_t, sparse_linear_solver)
Cgal_add_pmp_parameter_and_enum(number_of_iterations_t, number_of_iterations)
Cgal_add_pmp_parameter_and_enum(number_of_relaxation_steps_t, number_of_relaxation_steps)
Cgal_add_pmp_parameter_and_enum(protect_constraints_t, protect_constraints)
Cgal_add_pmp_parameter_and_enum(relax_constraints_t, relax_constraints)
Cgal_add_pmp_parameter_and_enum(vertex_is_constrained_t, vertex_is_constrained)
Cgal_add_pmp_parameter_and_enum(face_patch_t, face_patch)
Cgal_add_pmp_parameter_and_enum(random_uniform_sampling_t, random_uniform_sampling)
Cgal_add_pmp_parameter_and_enum(grid_sampling_t, grid_sampling)
Cgal_add_pmp_parameter_and_enum(monte_carlo_sampling_t, monte_carlo_sampling)
Cgal_add_pmp_parameter_and_enum(do_sample_edges_t, do_sample_edges)
Cgal_add_pmp_parameter_and_enum(do_sample_vertices_t, do_sample_vertices)
Cgal_add_pmp_parameter_and_enum(do_sample_faces_t, do_sample_faces)
Cgal_add_pmp_parameter_and_enum(number_of_points_on_faces_t, number_of_points_on_faces)
Cgal_add_pmp_parameter_and_enum(number_of_points_per_face_t, number_of_points_per_face)
Cgal_add_pmp_parameter_and_enum(grid_spacing_t, grid_spacing)
Cgal_add_pmp_parameter_and_enum(number_of_points_per_edge_t, number_of_points_per_edge)
Cgal_add_pmp_parameter_and_enum(number_of_points_on_edges_t, number_of_points_on_edges)

//internal
Cgal_add_pmp_parameter_and_enum(weight_calculator_t, weight_calculator)
Cgal_add_pmp_parameter(boost::vertex_point_t, vertex_point_map)
Cgal_add_pmp_parameter(edge_is_constrained_t, edge_is_constrained_map)
Cgal_add_pmp_parameter(edge_is_constrained_params_t, edge_is_constrained_map_params)
Cgal_add_pmp_parameter(boost::face_index_t, face_index_map)
Cgal_add_pmp_parameter(boost::vertex_index_t, vertex_index_map)

#undef Cgal_add_pmp_parameter
#undef Cgal_add_pmp_parameter_and_enum
namespace CGAL{

  enum all_default_t { all_default} ; //cannot use macro because of inline
  enum nb_points_per_area_unit_t    { nb_points_per_area_unit }; //cannot use macro because names diverge
  enum nb_points_per_distance_unit_t{ nb_points_per_distance_unit }; //cannot use macro because names diverge
  //to be documented
  enum face_normal_t { face_normal }; //cannot use macro because names diverge

  template <typename T, typename Tag, typename Base>
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

    Cgal_add_pmp_parameter_in_struct(density_control_factor_t, density_control_factor)
    Cgal_add_pmp_parameter_in_struct(use_delaunay_triangulation_t, use_delaunay_triangulation)
    Cgal_add_pmp_parameter_in_struct(fairing_continuity_t, fairing_continuity)
    Cgal_add_pmp_parameter_in_struct(sparse_linear_solver_t, sparse_linear_solver)
    Cgal_add_pmp_parameter_in_struct(number_of_iterations_t, number_of_iterations)
    Cgal_add_pmp_parameter_in_struct(number_of_relaxation_steps_t, number_of_relaxation_steps)
    Cgal_add_pmp_parameter_in_struct(protect_constraints_t, protect_constraints)
    Cgal_add_pmp_parameter_in_struct(relax_constraints_t, relax_constraints)
    Cgal_add_pmp_parameter_in_struct(vertex_is_constrained_t, vertex_is_constrained)
    Cgal_add_pmp_parameter_in_struct(face_patch_t, face_patch)
    Cgal_add_pmp_parameter_in_struct(random_uniform_sampling_t, random_uniform_sampling)
    Cgal_add_pmp_parameter_in_struct(grid_sampling_t, grid_sampling)
    Cgal_add_pmp_parameter_in_struct(monte_carlo_sampling_t, monte_carlo_sampling)
    Cgal_add_pmp_parameter_in_struct(do_sample_edges_t, do_sample_edges)
    Cgal_add_pmp_parameter_in_struct(do_sample_vertices_t, do_sample_vertices)
    Cgal_add_pmp_parameter_in_struct(do_sample_faces_t, do_sample_faces)
    Cgal_add_pmp_parameter_in_struct(number_of_points_on_faces_t, number_of_points_on_faces)
    Cgal_add_pmp_parameter_in_struct(number_of_points_per_face_t, number_of_points_per_face)
    Cgal_add_pmp_parameter_in_struct(grid_spacing_t, grid_spacing)
    Cgal_add_pmp_parameter_in_struct(nb_points_per_area_unit_t, nb_points_per_area_unit)
    Cgal_add_pmp_parameter_in_struct(number_of_points_per_edge_t, number_of_points_per_edge)
    Cgal_add_pmp_parameter_in_struct(number_of_points_on_edges_t, number_of_points_on_edges)
    Cgal_add_pmp_parameter_in_struct(nb_points_per_distance_unit_t, nb_points_per_distance_unit)
    Cgal_add_pmp_parameter_in_struct(face_normal_t, face_normal_map)
    Cgal_add_pmp_parameter_in_struct(weight_calculator_t, weight_calculator)
    Cgal_add_pmp_parameter_in_struct(geom_traits_t, geom_traits)
    Cgal_add_pmp_parameter_in_struct(boost::vertex_point_t, vertex_point_map)
    Cgal_add_pmp_parameter_in_struct(edge_is_constrained_t, edge_is_constrained_map)
    Cgal_add_pmp_parameter_in_struct(edge_is_constrained_params_t, edge_is_constrained_map_params)
    Cgal_add_pmp_parameter_in_struct(boost::face_index_t, face_index_map)
    Cgal_add_pmp_parameter_in_struct(boost::vertex_index_t, vertex_index_map)
    #undef Cgal_add_pmp_parameter_in_struct
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

  template <typename FaceNormalMap>
    pmp_bgl_named_params<FaceNormalMap, face_normal_t>
    face_normal_map(const FaceNormalMap& m)
    {
      typedef pmp_bgl_named_params<FaceNormalMap, face_normal_t> Params;
      return Params(m);
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
