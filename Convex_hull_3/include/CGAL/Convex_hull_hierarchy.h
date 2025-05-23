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

#include <CGAL/Surface_mesh.h>

#include <vector>

#ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
      size_t nb_visited=0;
#endif

namespace CGAL{

  /// \ingroup PkgConvexHull3Ref
  /// This class wraps a convex hull with a data structure optimized for finding the extreme point of the convex
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

  /// @private
  using Vertex_index = typename PolygonMesh::Vertex_index;

  /// The point type.
  typedef typename PolygonMesh::Point Point;

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
  */
  template <typename VertexListGraph, typename NamedParameters=parameters::Default_named_parameters>
  Convex_hull_hierarchy(const VertexListGraph &g, const NamedParameters & np=parameters::default_values){
    PolygonMesh ch;
    convex_hull_3(g, ch, np);
    hierarchy_sm.push_back(std::move(ch));
    // init_hierarchy(); //Hierarchy just need a property map
  }

  /// @private
  Convex_hull_hierarchy(const PolygonMesh &sm){
    using Traits=typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point>::type;
    PolygonMesh ch;
    convex_hull_3(sm, ch);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(Traits());
  };

  /// @private
  template <typename PointInput, typename Traits>
  Convex_hull_hierarchy(const Surface_mesh<PointInput> &sm, const Traits &traits){
    PolygonMesh ch;
    convex_hull_3(vertices(sm).begin(), vertices(sm).end(), ch, traits);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(traits);
  };

   /**
  * constructor taking the points in the range `[first, last)`.
  *
  * @tparam RangeIterator must be an input iterator with a value type equivalent to `Traits::Point_3`
  * @tparam Traits must be a model of the concept ConvexHullTraits_3. For the purposes of checking the postcondition that the convex hull is valid,
  *         `Traits` must also be a model of the concept `IsStronglyConvexTraits_3`. Furthermore, `Traits` must define a type `PolygonMesh` that is a model of `MutableFaceGraph`.
  */
  template <typename RangeIterator, typename Traits = typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point>::type>
  Convex_hull_hierarchy(RangeIterator begin, RangeIterator end, const Traits &traits=Traits()){
    PolygonMesh ch;
    convex_hull_3(begin, end, ch, traits);
    hierarchy_sm.push_back(std::move(ch));
    init_hierarchy(traits);
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
  const Point extreme_point(const Direction_3 &dir, const NamedParameters & np=parameters::default_values) const {
    // return extreme_point(dir.vector());
  }

  // /**
  // * Compute and return the extreme point of the convex hull along the provide direction.
  // *
  // * Using a converter, this function is typically used in filtered_predicate to compute the extreme_point with a different kernel than the one of `Point`.
  // * (TODO Maybe not documented this function)
  // *
  // * @tparam NamedParameters: a sequence of named parameters
  // * @tparam Traits: A traits class providing ` Vector_3`, `Point_3`, `ConstructVector_3` and `CrossProductVector_3`
  // * @tparam Converter: A converter that converts the type `Point` to `Traits::Point_3`
  // *
  // */
  /// @private
  template <class Converter, class Vector_3>
  const Point extreme_point(const Vector_3 &dir, const Converter& converter) const {
    using Point_3 = typename Kernel_traits<Vector_3>::Kernel::Point_3;
    using K= typename Kernel_traits<Vector_3>::Kernel;
    using FT= typename K::FT;

    size_t level=maxlevel();

    const PolygonMesh &sm = hierarchy_sm[maxlevel()];

    Vertex_index argmax=*sm.vertices().begin();
    FT tmax=K().construct_vector_3_object()(ORIGIN, converter(sm.point(argmax)))*dir;

    if(sm.vertices().size() <= MAXSIZE_FOR_NAIVE_SEARCH){
      //If maxlevel is small, we simply go through all its vertices
      for(auto vh=++(sm.vertices().begin()); vh!=sm.vertices().end(); ++vh){
        Vertex_index v=*vh;
  #ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
        ++nb_visited;
  #endif
        FT p=K().construct_vector_3_object()(ORIGIN, converter(sm.point(v)))*dir;
        if(compare(tmax, p)==SMALLER){
          tmax=p;
          argmax=v;
        }
      }

      // Go to under level
      if(level>0){
        argmax=(*(sm.template property_map<Vertex_index, Vertex_index>("v:next_in_hierarchy")))[argmax];
        --level;
      } else {
        return sm.point(argmax);
      }
    }

    for(; level>=0; --level){
      // Starting from the vertex of the previous level, we walk on the graph
      // along neighbors that increase the "score"
      const PolygonMesh &csm = mesh(level);
      bool is_local_max;
      do{
        is_local_max=true;
        for(Vertex_index v: vertices_around_target(argmax ,csm)){
    #ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
          ++nb_visited;
    #endif
          FT p=K().construct_vector_3_object()(ORIGIN, converter(csm.point(v)))*dir;
          if(compare(tmax, p)==SMALLER){
            tmax=p;
            argmax=v;
            is_local_max=false; // repeat with the new vertex
            break;
          }
        }
      } while(!is_local_max);
      if(level>0)
        argmax=(*(csm.template property_map<Vertex_index, Vertex_index>("v:next_in_hierarchy")))[argmax];
      else
        return csm.point(argmax);
    }
  }

private:
  std::vector< PolygonMesh > hierarchy_sm;

  template <typename Traits>
  void init_hierarchy(const Traits &traits){
    size_t size=hierarchy_sm[0].vertices().size();
    size_t level=0;

    if(size<=MINSIZE_FOR_NEXT_LEVEL)
      return;

    Random rng;
    while(size>MINSIZE_FOR_NEXT_LEVEL){
      std::vector<Point> select_points;
      std::map<Point, Vertex_index> select_vertices;
      select_points.reserve(2*size/RATIO);
      PolygonMesh& above_sm=hierarchy_sm[level];
      ++level;

      for(Vertex_index v: above_sm.vertices())
        if(rng.get_int(0,RATIO-1)==0){
          select_points.push_back(above_sm.point(v));
          select_vertices[above_sm.point(v)]=v;
        }

      PolygonMesh new_sm;
      convex_hull_3(select_points.begin(), select_points.end(), new_sm, traits);
      typename PolygonMesh::Property_map<Vertex_index, Vertex_index> pm=new_sm.template add_property_map<Vertex_index>("v:next_in_hierarchy", *(above_sm.vertices().begin())).first;
      for(Vertex_index v : new_sm.vertices())
        pm[v]=select_vertices[new_sm.point(v)];
      size=new_sm.vertices().size();

      hierarchy_sm.push_back(std::move(new_sm));
    }
  }
};

// Point can be different from K::Point_3
template <class Point, class Converter, class Vector_3>
const Point extreme_point(const Convex_hull_hierarchy<Point> &C, const Vector_3 &dir, const Converter &c){
  return C.template extreme_point(dir, c);
}

}

#endif