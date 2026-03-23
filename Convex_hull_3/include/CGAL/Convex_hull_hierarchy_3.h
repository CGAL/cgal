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
#include <CGAL/extreme_point_3.h>
#include <vector>

namespace CGAL{

  /// \ingroup PkgConvexHull3Intersections
  /// \brief This class stores a convex hull in a data structure optimized for fast intersection tests.
  ///
  /// More specifically, the structure is optimized for `CGAL::extreme_vertex_3()`, which is used by `CGAL::Convex_hull_3::do_intersect()`.
  /// The computational complexities of `CGAL::extreme_point_3()` and consequently `CGAL::Convex_hull_3::do_intersect()` is \f$O(n)\f$ for a range of points,
  /// \f$O(\sqrt{n})\f$ on average for a graph, but only \f$O(\log{n})\f$ on average for a `Convex_hull_hierarchy_3`, where \f$n\f$ is the number of points of the object.
  ///
  /// Building this structure has linear complexity and is faster than computing a convex hull, but more costly than a single call to `CGAL::Convex_hull_3::do_intersect()`.
  /// It is therefore relevant when many intersection queries are performed, particularly when a convex hull has a large number of vertices.
  ///
  /// To evaluate this class, we generated convex hulls by sampling \f$n\f$ random points on the unit sphere, computing their convex hulls, and performing intersection tests between them.
  ///
  /// Average performance measurements (times in milliseconds).
  ///
  /// | Sphere size | Hull build | Hierarchy build | do_intersect (Range) | do_intersect (Graph) | do_intersect (Hull Hierarchy) |
  /// |-------------|-----------:|----------------:|---------------------:|---------------------:|------------------------------:|
  /// | 40          | 0.06431  | 0.006908 | 0.00260 | 0.00174 | 0.00170  |
  /// | 100         | 0.18069  | 0.016147 | 0.00594 | 0.00260 | 0.00272  |
  /// | 350         | 0.68038  | 0.076661 | 0.01976 | 0.00735 | 0.00697  |
  /// | 1,000       | 2.14852  | 0.207484 | 0.05664 | 0.02079 | 0.01292  |
  /// | 3,500       | 7.49052  | 0.740683 | 0.19527 | 0.04876 | 0.01842  |
  /// | 10,000      | 22.7057  | 2.22908  | 0.54977 | 0.07595 | 0.02789  |
  /// | 35,000      | 93.7957  | 8.512    | 2.00917 | 0.13722 | 0.03730  |
  /// | 100,000     | 399.377  | 32.0424  | 5.50728 | 0.39127 | 0.05841  |
  ///
  /// @tparam PolygonMesh The polygon mesh used to construct each level of the hierarchy. Must be a model of `VertexListGraph` and `MutableFaceGraph`.
  ///         An internal property map for `CGAL::vertex_point_t` must be available, with a value type that is a `Point_3` typedef of a \cgal kernel.
#if DOXYGEN_RUNNING
template < class PolygonMesh>
#else
template < class PolygonMesh,
           int MAX_HIERARCHY_DEPTH = 5,
           int RATIO_BETWEEN_LEVELS = 24,
           int MINSIZE_FOR_NEXT_LEVEL = RATIO_BETWEEN_LEVELS*6>
#endif
struct Convex_hull_hierarchy_3{
  // parameterization of the hierarchy
  /// @private
  constexpr static size_t RATIO = RATIO_BETWEEN_LEVELS;
  /// @private
  constexpr static size_t MAXSIZE_FOR_NAIVE_SEARCH = RATIO;

  /// @private
  using VPM = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type;

  /// @private
  using Point_3 = typename boost::property_traits<VPM>::value_type;

  /// @private
  using GT = typename Kernel_traits<Point_3>::Kernel;

  /// @private
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  /// @private
  using V2VMap = typename boost::property_map<PolygonMesh, dynamic_vertex_property_t<vertex_descriptor> >::type;

  /**
  * Constructor taking a graph `g`. If `compute_convex_hull` is set to `true`, compute and use the convex hull of `g` instead.
  *
  * @tparam Graph a model of `VertexListGraph`
  * @tparam NamedParameters a sequence of named parameters
  *
  * @param g the graph
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `g`}
  *     \cgalParamType{a model of `ReadablePropertyMap` whose value type is a point type}
  *     \cgalParamDefault{If this parameter is omitted, an internal property map for ` CGAL::vertex_point_t`  must be available in `VertexListGraph`}
  *   \cgalParamNEnd
  *   \cgalParamNBegin{compute_convex_hull}
  *     \cgalParamDescription{If set to `true`, the convex hull of `g` is computed. Otherwise, `g` must be convex.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *   \cgalParamNBegin{random_seed}
  *     \cgalParamDescription{Define the seed used by the hierarchy.}
  *     \cgalParamType{unsigned int}
  *     \cgalParamDefault{The seed used by `CGAL::get_default_random()`.}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <typename Graph,
            typename NamedParameters=parameters::Default_named_parameters>
  Convex_hull_hierarchy_3(const Graph &g, const NamedParameters& np = parameters::default_values()){
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    PolygonMesh ch;
    bool compute_convex_hull = choose_parameter(get_parameter(np, internal_np::compute_convex_hull), true);
    size_t seed = choose_parameter(get_parameter(np, internal_np::random_seed), 0);
    Random rng = (is_default_parameter<NamedParameters, internal_np::random_seed_t>::value)? get_default_random(): Random(seed);
    if(compute_convex_hull)
      convex_hull_3(g, ch, np);
    else
      copy_face_graph(g, ch, np);
    hierarchy_sm.reserve(MAX_HIERARCHY_DEPTH);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(rng);
  };

   /**
  * Constructor taking the points in the range `[first, last)`.
  *
  * @tparam PointIterator is an input iterator with a value type equivalent to `Traits::Point_3`
  * @tparam Traits is a model of the concept ConvexHullTraits_3. For the purposes of checking the postcondition that the convex hull is valid,
  *         `Traits` must also be a model of the concept `IsStronglyConvexTraits_3`.
  */
  template <typename PointIterator, typename Traits = Convex_hull_traits_3<GT , PolygonMesh> >
  Convex_hull_hierarchy_3(PointIterator begin, PointIterator end, const Traits &traits = Traits()){
    PolygonMesh ch;
    convex_hull_3(begin, end, ch, traits);
    hierarchy_sm.reserve(MAX_HIERARCHY_DEPTH);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(get_default_random(), traits);
  };

  // /// Returns the index of the deepest level of the hierarchy (the one with the fewest vertices).
  // /// This corresponds to the number of levels minus one.
  /// @private
  std::size_t maxlevel() const{
    return hierarchy_sm.size()-1;
  }

  // /// returns the mesh of the corresponding level of the hierarchy.
  /// @private
  const PolygonMesh& mesh(std::size_t level) const{
    return hierarchy_sm[level];
  }

  /// returns a const reference to the mesh stored by the class.
  const PolygonMesh& mesh() const{
    return hierarchy_sm[0];
  }

  /// @private
  template <typename Direction>
  typename Kernel_traits<Direction>::Kernel::Point_3 extreme_point_3(const Direction &dir) const {
    using Dir_GT = typename Kernel_traits<Direction>::Kernel;
    using GTC = Cartesian_converter<GT, Dir_GT>;
    GTC converter;
    VPM vpm = get_const_property_map(vertex_point, hierarchy_sm[0]);
    return converter(get(vpm, extreme_vertex_3(dir)));
  }

  /**
  * returns the vertex of the convex hull whose geometric position projected onto that direction is maximal.
  * If not unique, a single vertex is returned.
  *
  * @tparam Direction must be a `Direction_3` typedef of a \cgal kernel
  *
  * @param dir the direction
  */
  template <typename Direction>
  vertex_descriptor extreme_vertex_3(const Direction &dir) const
  {
    using Dir_GT = typename Kernel_traits<Direction>::Kernel;
    using GTC = Cartesian_converter<GT, Dir_GT>;
    Dir_GT gt;
    GTC converter;

    auto vector_3 = gt.construct_vector_3_object();
    auto csp = gt.compare_scalar_product_3_object();

    VPM vpm = get_const_property_map(vertex_point, hierarchy_sm[0]);
    const typename Dir_GT::Vector_3 &vdir = dir.vector();

    size_t level=maxlevel();

    const PolygonMesh &sm = hierarchy_sm[level];
    if(level == 0)
      return CGAL::extreme_vertex_3(sm, dir);

    vertex_descriptor argmax = *vertices(sm).begin();
    Vector_3 vec_max = vector_3(ORIGIN, converter(get(vpm, get(to_base_maps[level], argmax))));
    if(vertices(sm).size() <= MAXSIZE_FOR_NAIVE_SEARCH){
      //If maxlevel is small, we simply go through all its vertices
      for(auto vh=++(vertices(sm).begin()); vh!=vertices(sm).end(); ++vh){
        vertex_descriptor v = *vh;
        Vector_3 vec = vector_3(ORIGIN, converter(get(vpm, get(to_base_maps[level], v))));
        if(csp(vec_max, vdir, vec, vdir)==SMALLER){
          vec_max=vec;
          argmax=v;
        }
      }

      // Go to level below
      if(level>0){
        argmax= get(next_in_hierarchy_maps[level], argmax);
        --level;
      } else {
        return argmax;
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
        for(vertex_descriptor v: vertices_around_target(argmax ,csm)){
          Vector_3 vec = vector_3(ORIGIN, converter(get(vpm, get(cbase, v))));
          if(csp(vec_max, vdir, vec, vdir)==SMALLER){
            vec_max = vec;
            argmax = v;
            is_local_max = false; // repeat with the new vertex
            break;
          }
        }
      } while(!is_local_max);
      if(level>0){
        argmax= get(cnext, argmax);
      } else {
        return argmax;
      }
    }
  }

private:
  template <typename Random, typename Traits = Convex_hull_traits_3<GT , PolygonMesh> >
  void init_hierarchy(Random &rng = get_default_random(), const Traits &traits = Traits()){
    VPM vpm = get_const_property_map(vertex_point, hierarchy_sm[0]);

    size_t size=vertices(hierarchy_sm[0]).size();
    size_t level=0;

    to_base_maps.reserve(MAX_HIERARCHY_DEPTH);
    next_in_hierarchy_maps.reserve(MAX_HIERARCHY_DEPTH);

    V2VMap to_base_map = get(dynamic_vertex_property_t<vertex_descriptor>(), hierarchy_sm[0]);
    for(vertex_descriptor v : vertices(hierarchy_sm[0]))
      put(to_base_map, v, v);
    next_in_hierarchy_maps.push_back(V2VMap());
    to_base_maps.push_back(std::move(to_base_map));

    while(size>MINSIZE_FOR_NEXT_LEVEL && level<MAX_HIERARCHY_DEPTH-1){
      // Create a new level
      std::vector<Point_3> select_points;
      std::map<Point_3, vertex_descriptor> select_vertices;
      select_points.reserve(2*size/RATIO);
      PolygonMesh& above_sm=hierarchy_sm[level];
      ++level;

      // Select randomly vertices of the new level
      while(select_points.size()==0)
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
  // Given a vertex of level n, give the correspondent vertex of level n-1
  std::vector<V2VMap> next_in_hierarchy_maps;
  // Given a vertex of level n, give the correspondent vertex of level 0
  std::vector<V2VMap> to_base_maps;
};

/**
* \ingroup PkgConvexHull3Queries
*
* returns the vertex of the convex hull whose geometric position projected onto that direction is maximal.
* If not unique, a single vertex is returned.
*
* @tparam Mesh a model of `VertexListGraph` and `MutableFaceGraph`
* @tparam Direction is a model of `Kernel::Direction_3`.
*
* @param ch the convex hull
* @param dir the direction
*
*/
template <class PolygonMesh, class Direction>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
extreme_vertex_3(const Convex_hull_hierarchy_3<PolygonMesh> &ch, const Direction &dir)
{
  return ch.extreme_vertex_3(dir);
}

}

#endif // CGAL_CONVEX_HULL_HIERARCHY_3_H
