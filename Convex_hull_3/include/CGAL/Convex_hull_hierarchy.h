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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/license/Convex_hull_3.h>
#include <vector>

#ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
      size_t nb_visited=0;
#endif

template <typename, typename = void>
inline constexpr bool has_property_map = false;

template <typename T>
inline constexpr bool has_property_map<T, std::void_t<typename T::Property_map> > = true;

namespace CGAL{

  /// \ingroup PkgConvexHull3Ref
  /// This class implements a convex hull with a data structure optimized for finding the extreme point of the convex
  /// hull in a given direction. In particular, this operation is called by `CGAL::Convex_hull_3::do_intersect()` and therefore, this class
  ///  is optimized for very fast intersection tests.
  ///
  /// @tparam PolygonMesh The polygon mesh structure used to construct each level of the hierarchy. Must be a model of ` MutableFaceGraph`.
  ///         An internal property map for  ` CGAL::vertex_point_t` must be available, from which the point type `Point` is deduced.
  ///         There is no requirement on `Point`, besides being default constructible and assignable.
  ///         In typical use cases it will be a 3D point type.
template < class PolygonMesh>
struct Convex_hull_hierarchy{
  // parameterization of the hierarchy
  /// @private
  constexpr static size_t RATIO = 32;
  /// @private
  constexpr static size_t MINSIZE_FOR_NEXT_LEVEL = RATIO*4;
  /// @private
  constexpr static size_t MAXSIZE_FOR_NAIVE_SEARCH = RATIO;

  /// The point type.
  typedef typename PolygonMesh::Point Point;

  /// The mesh type.
  typedef PolygonMesh Mesh;

  /// @private
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  /// @private
  typedef typename PolygonMesh::Property_map<vertex_descriptor, vertex_descriptor> Property_map;

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
  * @tparam Traits: model of the concept ConvexHullTraits_3. For the purposes of checking the postcondition that the convex hull is valid, Traits must also be a model of the concept IsStronglyConvexTraits_3.
  *
  * @param g the graph
  * @param np an optional sequence of `Named Parameters` among the ones listed below
  * @param traits an instance of Traits
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of  `g`}
  *     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
  *     \cgalParamDefault{If this parameter is omitted, an internal property map for ` CGAL::vertex_point_t`  must be available in `VertexListGraph`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <typename Graph,
            typename NamedParameters=parameters::Default_named_parameters ,
            typename Traits = typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point>::type>
  Convex_hull_hierarchy(const Graph &g, const NamedParameters& np = parameters::default_values(), const Traits &traits=Traits()){
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    typedef typename GetVertexPointMap<Graph, NamedParameters>::const_type Vpmap;
    typedef CGAL::Property_map_to_unary_function<Vpmap> Vpmap_fct;
    Vpmap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_const_property_map(boost::vertex_point, g));
    Vpmap_fct v2p(vpm);
    // Convex_hull_hierarchy(boost::make_transform_iterator(vertices(g).begin(), v2p),
    //                       boost::make_transform_iterator(vertices(g).end(), v2p), traits);
    //                       PolygonMesh ch;
    PolygonMesh ch;
    convex_hull_3(boost::make_transform_iterator(vertices(g).begin(), v2p),
                  boost::make_transform_iterator(vertices(g).end(), v2p), ch, traits);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(traits);
  };

  /// @private
  template <typename Graph,
            typename Traits = typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point>::type>
  Convex_hull_hierarchy(const Graph &g, const Traits &traits): Convex_hull_hierarchy(g, parameters::default_values(), traits){};

  /**
  * constructor taking the points in the range `[first, last)`.
  *
  * @tparam RangeIterator must be an input iterator with a value type equivalent to `Traits::Point_3`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam Traits model of the concept ConvexHullTraits_3. For the purposes of checking the postcondition that the convex hull is valid, Traits must also be a model of the concept IsStronglyConvexTraits_3.
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{point_map}
  *     \cgalParamDescription{a property map associating points to the elements of `r`}
  *     \cgalParamType{a model of `ReadablePropertyMap`}
  *     \cgalParamDefault{`CGAL::Identity_property_map`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * If the kernel R of the points determined by the value type of `RangeIterator` is a kernel with exact predicates but inexact constructions (in practice we check R::Has_filtered_predicates_tag is Tag_true and R::FT is a floating point type), then the default traits class of `convex_hull_3()` is `Convex_hull_traits_3<R>`, and R otherwise.
  *
  */
   /**
  * constructor taking the points in the range `[first, last)`.
  *
  * @tparam RangeIterator must be an input iterator with a value type equivalent to `Traits::Point_3`
  * @tparam Traits must be a model of the concept ConvexHullTraits_3. For the purposes of checking the postcondition that the convex hull is valid,
  *         `Traits` must also be a model of the concept `IsStronglyConvexTraits_3`. Furthermore, `Traits` must define a type `PolygonMesh` that is a model of `MutableFaceGraph`.
  */
  template <typename RangeIterator,
            typename Traits = typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point>::type>
  Convex_hull_hierarchy(RangeIterator begin, RangeIterator end, const Traits &traits=Traits()){
    PolygonMesh ch;
    convex_hull_3(begin, end, ch, traits);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(traits);
    CGAL_assertion(hierarchy_sm.size()!=0);
  };

  /**
  * computes the furthest point of the convex hull along the direction.
  *
  * @tparam NamedParameters: a sequence of named parameters
  *
  * @param dir the direction
  * @param np an optional sequence of `Named Parameters` among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of the convex hull}
  *     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
  *     \cgalParamDefault{If this parameter is omitted, `Point` must be a 3D point type of the same kernel than `dir`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  */
  template <class Direction_3, typename NamedParameters=parameters::Default_named_parameters>
  typename Kernel_traits<Direction_3>::Kernel::Point_3 extreme_point_3(const Direction_3 &dir, const NamedParameters & np=parameters::default_values) const {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    using vertex_descriptor= typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

    using GetVertexPointMap = GetVertexPointMap<PolygonMesh, NamedParameters>;
    using VPM = typename GetVertexPointMap::const_type;
    // VPM pm = GetVertexPointMap::get_const_map(np, hierarchy_sm[0]);

    using GetGeomTraits = GetGeomTraits<PolygonMesh, NamedParameters>;
    using Mesh_GT= typename GetGeomTraits::type;
    using Default_GT = typename Kernel_traits<Direction_3>::Kernel;
    using GT=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_t,
        NamedParameters,
        Default_GT
      > ::type;
    GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
    using FT= typename GT::FT;

    using Default_geom_traits_converter = Cartesian_converter<Mesh_GT, GT>;
    using GTC=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_converter_t,
        NamedParameters,
        Default_geom_traits_converter
      > ::type;
    GTC converter = choose_parameter<GTC>(get_parameter(np, internal_np::geom_traits_converter));

    size_t level=maxlevel();

    const PolygonMesh &sm = hierarchy_sm[maxlevel()];
    const auto &init_pm = sm.points();
    // const auto &pm = make_compose_property_map(sm.points(), point_map);

    vertex_descriptor argmax=*vertices(sm).begin();
    Vector_3 vec=gt.construct_vector_3_object()(ORIGIN, converter(get(init_pm, argmax)));
    FT sp_max=gt.compute_scalar_product_3_object()(vec,dir.vector());
    if(vertices(sm).size() <= MAXSIZE_FOR_NAIVE_SEARCH){
      const auto &pm = sm.points();
      //If maxlevel is small, we simply go through all its vertices
      for(auto vh=++(vertices(sm).begin()); vh!=vertices(sm).end(); ++vh){
        vertex_descriptor v=*vh;
  #ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
        ++nb_visited;
  #endif
        vec=gt.construct_vector_3_object()(ORIGIN, converter(get(pm, v)));
        FT sp=gt.compute_scalar_product_3_object()(vec,dir.vector());
        if(compare(sp_max, sp)==SMALLER){
          sp_max=sp;
          argmax=v;
        }
      }

      // Go to under level
      if(level>0){
        argmax=(*(sm.template property_map<vertex_descriptor, vertex_descriptor>("v:next_in_hierarchy")))[argmax];
        --level;
      } else {
        return converter(get(pm,argmax));
      }
    }

    for(; level>=0; --level){
      // Starting from the vertex of the previous level, we walk on the graph
      // along neighbors that increase the "score"
      const PolygonMesh &csm = mesh(level);
      const auto &pm = csm.points();//make_compose_property_map(csm.points(), point_map);
      bool is_local_max;
      do{
        is_local_max=true;
        for(vertex_descriptor v: vertices_around_target(argmax ,csm)){
    #ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
          ++nb_visited;
    #endif
          vec=gt.construct_vector_3_object()(ORIGIN, converter(get(pm,v)));
          FT sp=gt.compute_scalar_product_3_object()(vec,dir.vector());
          if(compare(sp_max, sp)==SMALLER){
            sp_max=sp;
            argmax=v;
            is_local_max=false; // repeat with the new vertex
            break;
          }
        }
      } while(!is_local_max);
      if(level>0)
        argmax=(*(csm.template property_map<vertex_descriptor, vertex_descriptor>("v:next_in_hierarchy")))[argmax];
      else
        return converter(get(pm,argmax));
    }
  }

private:
  std::vector< PolygonMesh > hierarchy_sm;

  template <typename Traits>
  void init_hierarchy(const Traits &traits){
    size_t size=vertices(hierarchy_sm[0]).size();
    size_t level=0;
    if(size<=MINSIZE_FOR_NEXT_LEVEL)
      return;

    Random rng;
    while(size>MINSIZE_FOR_NEXT_LEVEL){
      std::vector<Point> select_points;
      std::map<Point, vertex_descriptor> select_vertices;
      select_points.reserve(2*size/RATIO);
      PolygonMesh& above_sm=hierarchy_sm[level];
      ++level;

      for(vertex_descriptor v: vertices(above_sm))
        if(rng.get_int(0,RATIO-1)==0){
          select_points.push_back(above_sm.point(v));
          select_vertices[above_sm.point(v)]=v;
        }

      PolygonMesh new_sm;
      convex_hull_3(select_points.begin(), select_points.end(), new_sm, traits);
      Property_map pm=new_sm.template add_property_map<vertex_descriptor>("v:next_in_hierarchy", *(vertices(above_sm).begin())).first;

      for(vertex_descriptor v : vertices(new_sm))
        put(pm,v,select_vertices[new_sm.point(v)]);
      size=vertices(new_sm).size();

      CGAL_assertion(size!=0);
      hierarchy_sm.push_back(std::move(new_sm));
    }
  }
};

}

#endif