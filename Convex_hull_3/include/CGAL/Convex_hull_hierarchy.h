// Copyright (c) 2025 Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque
//

#ifndef CGAL_CONVEX_HULL_HIERARCHY_3_H
#define CGAL_CONVEX_HULL_HIERARCHY_3_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/property_map.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

#ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
      size_t nb_visited=0;
#endif

namespace CGAL{

  /// \ingroup PkgConvexHull3Ref
  /// This class stores a convex hull with a data structure optimized for finding the extreme point of the convex
  /// hull in a given direction. In particular, this operation is called by `CGAL::Convex_hull_3::do_intersect()` and therefore, this class
  ///  is optimized for very fast intersection tests.
  ///
  /// @tparam PolygonMesh The polygon mesh structure used to construct each level of the hierarchy. Must be a model of `MutableFaceGraph`.
  ///         An internal property map for  ` CGAL::vertex_point_t` must be available.
template < class PolygonMesh >
struct Convex_hull_hierarchy{
  // parameterization of the hierarchy
  /// @private
  constexpr static size_t RATIO = 24;
  /// @private
  constexpr static size_t MINSIZE_FOR_NEXT_LEVEL = RATIO*6;
  /// @private
  constexpr static size_t MAXSIZE_FOR_NAIVE_SEARCH = RATIO;

  /// The mesh type.
  using Mesh = PolygonMesh;

  /// @private
  using VPM = typename boost::property_map<Mesh, CGAL::vertex_point_t>::const_type;

  /// @private
  using Point_3 = typename boost::property_traits<VPM>::value_type;

  /// @private
  using GT = typename Kernel_traits<Point_3>::Kernel;

  /// @private
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  /// @private
  using V2VMap = typename boost::property_map<PolygonMesh, dynamic_vertex_property_t<vertex_descriptor> >::type;

  // /// returns the number of the higher level of the hierarchy.
  /// @private
  std::size_t maxlevel() const{
    return hierarchy_sm.size()-1;
  }

  // /// returns the mesh of the corresponding level of the hierarchy.
  /// @private
  const PolygonMesh& mesh(std::size_t level) const{
    return hierarchy_sm[level];
  }

  /**
  * constructor taking the points associated to the vertices of  `g`.
  *
  * @tparam VertexListGraph: a model of `VertexListGraph`
  * @tparam NamedParameters: a sequence of named parameters
  *
  * @param g the graph
  * @param np an optional sequence of `Named Parameters` among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of  `g`}
  *     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
  *     \cgalParamDefault{If this parameter is omitted, an internal property map for ` CGAL::vertex_point_t`  must be available in `VertexListGraph`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{random_seed}
  *     \cgalParamDescription{Define the seed used by}
  *     \cgalParamType{unsigned int}
  *     \cgalParamDefault{The seed used by `CGAL::get_default_random()`.}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <typename Graph,
            typename NamedParameters=parameters::Default_named_parameters>
  Convex_hull_hierarchy(const Graph &g, const NamedParameters& np = parameters::default_values()){
    PolygonMesh ch;
    convex_hull_3(g, ch, np);
    hierarchy_sm.reserve(4);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy();
  };

   /**
  * constructor taking the points in the range `[first, last)`.
  *
  * @tparam RangeIterator must be an input iterator with a value type equivalent to `Traits::Point_3`
  * @tparam Traits must be a model of the concept ConvexHullTraits_3. For the purposes of checking the postcondition that the convex hull is valid,
  *         `Traits` must also be a model of the concept `IsStronglyConvexTraits_3`. Furthermore, `Traits` must define a type `PolygonMesh` that is a model of `MutableFaceGraph`.
  */
  template <typename RangeIterator, typename Traits = Convex_hull_traits_3<GT , Mesh> >
  Convex_hull_hierarchy(RangeIterator begin, RangeIterator end, const Traits &traits = Traits()){
    PolygonMesh ch;
    convex_hull_3(begin, end, ch, traits);
    hierarchy_sm.reserve(4);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(traits);
  };

  ~Convex_hull_hierarchy(){
    to_base_maps.clear();
    next_in_hierarchy_maps.clear();
    for(int i=0; i<maxlevel(); ++i)
      clear(hierarchy_sm[i]);
  }

  /**
  * constructs the furthest point of the convex hull along the direction.
  *
  * @tparam `Direction_3` model of `Kernel::Direction_3`. The kernel does not require to be the same than the one using by the Mesh
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param dir the direction
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \return Return a `Point_3` using the same kernel as the direction.
  */
  template <typename Direction_3,
            typename NamedParameters=parameters::Default_named_parameters>
  typename Kernel_traits<Direction_3>::Kernel::Point_3 extreme_point_3(const Direction_3 &dir, const NamedParameters &np=parameters::default_values) const {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    using Default_GT = typename Kernel_traits<Direction_3>::Kernel;
    using IGT=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_t,
        NamedParameters,
        Default_GT
      > ::type;
    IGT gt = choose_parameter<IGT>(get_parameter(np, internal_np::geom_traits));
    using FT= typename IGT::FT;

    using Default_geom_traits_converter = Cartesian_converter<GT, IGT>;
    using GTC=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_converter_t,
        NamedParameters,
        Default_geom_traits_converter
      > ::type;
    GTC converter = choose_parameter<GTC>(get_parameter(np, internal_np::geom_traits_converter));

    VPM vpm = get_const_property_map(vertex_point, hierarchy_sm[0]);

    size_t level=maxlevel();

    const PolygonMesh &sm = hierarchy_sm[level];

    vertex_descriptor argmax=*vertices(sm).begin();
    Vector_3 vec=gt.construct_vector_3_object()(ORIGIN, converter(get(vpm, get(to_base_maps[level], argmax))));
    FT sp_max=gt.compute_scalar_product_3_object()(vec,dir.vector());
    if(vertices(sm).size() <= MAXSIZE_FOR_NAIVE_SEARCH){
      //If maxlevel is small, we simply go through all its vertices
      for(auto vh=++(vertices(sm).begin()); vh!=vertices(sm).end(); ++vh){
        vertex_descriptor v=*vh;
  #ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
        ++nb_visited;
  #endif
        vec=gt.construct_vector_3_object()(ORIGIN, converter(get(vpm, get(to_base_maps[level], v))));
        FT sp=gt.compute_scalar_product_3_object()(vec, dir.vector());
        if(compare(sp_max, sp)==SMALLER){
          sp_max=sp;
          argmax=v;
        }
      }

      // Go to level below
      if(level>0){
        argmax= get(next_in_hierarchy_maps[level], argmax);
        --level;
      } else {
        return converter(get(vpm, argmax));
      }
    }

    for(; true; --level){
      // Starting from the vertex of the previous level, we walk on the graph
      // along neighbors that increase the "score"
      const PolygonMesh &csm = hierarchy_sm[level];
      const V2VMap &cbase = to_base_maps[level];
      const V2VMap &cnext = next_in_hierarchy_maps[level];
      bool is_local_max;
      do{
        is_local_max=true;
        // CGAL_assertion(is_valid(argmax, csm));
        // CGAL_assertion(is_valid(halfedge(argmax, csm), csm));
        for(vertex_descriptor v: vertices_around_target(argmax ,csm)){
    #ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
          ++nb_visited;
    #endif
          vec=gt.construct_vector_3_object()(ORIGIN, converter(get(vpm, get(cbase, v))));
          FT sp=gt.compute_scalar_product_3_object()(vec,dir.vector());
          if(compare(sp_max, sp)==SMALLER){
            sp_max=sp;
            argmax=v;
            is_local_max=false; // repeat with the new vertex
            break;
          }
        }
      } while(!is_local_max);
      if(level>0){
        argmax= get(cnext, argmax);
      } else {
        return converter(get(vpm,argmax));
      }
    }
  }

private:
  template <typename Traits = Convex_hull_traits_3<GT , Mesh> >
  void init_hierarchy(const Traits &traits = Traits()){
    VPM vpm = get_const_property_map(vertex_point, hierarchy_sm[0]);

    size_t size=vertices(hierarchy_sm[0]).size();
    size_t level=0;

    to_base_maps.reserve(4);
    next_in_hierarchy_maps.reserve(4);

    V2VMap to_base_map = get(dynamic_vertex_property_t<vertex_descriptor>(), hierarchy_sm[0]);
    for(vertex_descriptor v : vertices(hierarchy_sm[0]))
      put(to_base_map, v, v);
    next_in_hierarchy_maps.push_back(V2VMap());
    to_base_maps.push_back(std::move(to_base_map));

    for(vertex_descriptor v: vertices(hierarchy_sm[0])){
      CGAL_assertion(get(to_base_maps[0], v) == v);
      CGAL_assertion(get(vpm, v) == get(make_compose_property_map(to_base_maps[0], vpm), v));
    }

    while(size>MINSIZE_FOR_NEXT_LEVEL){
      // Create a new level
      std::vector<Point_3> select_points;
      std::map<Point_3, vertex_descriptor> select_vertices;
      select_points.reserve(2*size/RATIO);
      PolygonMesh& above_sm=hierarchy_sm[level];
      ++level;

      // Select randomly vertices of the new level
      for(vertex_descriptor v: vertices(above_sm))
        if(rng.get_int(0,RATIO-1)==0){
          vertex_descriptor base = get(to_base_maps[level-1], v);
          select_points.emplace_back(get(vpm, base));
          select_vertices[select_points.back()] = v;
        }

      // Compute the graph of the new level
      PolygonMesh new_sm;
      convex_hull_3(select_points.begin(), select_points.end(), new_sm, traits);
      hierarchy_sm.push_back(std::move(new_sm));

      auto temp_vpm = get_property_map(vertex_point, hierarchy_sm.back());

      V2VMap to_base_map = get(CGAL::dynamic_vertex_property_t<vertex_descriptor>(),  hierarchy_sm.back());
      V2VMap next_in_hierarchy_map = get(CGAL::dynamic_vertex_property_t<vertex_descriptor>(),  hierarchy_sm.back());

      for(vertex_descriptor v : vertices( hierarchy_sm.back())){
        vertex_descriptor next = select_vertices[get(temp_vpm, v)];
        vertex_descriptor base = get(to_base_maps[level-1], next);
        put(next_in_hierarchy_map, v, next);
        put(to_base_map, v, base);
      }

      next_in_hierarchy_maps.push_back(std::move(next_in_hierarchy_map));
      to_base_maps.push_back(std::move(to_base_map));
      size=vertices( hierarchy_sm.back()).size();

      CGAL_assertion(size!=0);
    }
  }

  std::vector< PolygonMesh > hierarchy_sm;
  // Given a vertex of level n, give the correspondant vertex of level n-1
  std::vector<V2VMap> next_in_hierarchy_maps;
  // Given a vertex of level n, give the correspondant vertex of level 0
  std::vector<V2VMap> to_base_maps;
  Random rng;
};

}

#endif