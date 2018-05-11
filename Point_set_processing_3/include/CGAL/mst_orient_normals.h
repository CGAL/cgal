// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Laurent Saboret and Andreas Fabri

#ifndef CGAL_MST_ORIENT_NORMALS_H
#define CGAL_MST_ORIENT_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/trace.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_vertex_handle_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Index_property_map.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/use.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <list>
#include <climits>
#include <math.h>

#include <CGAL/property_map.h>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h> // work around a bug in boost 1.54
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {

/// Generalization of std::distance() to compute the distance between 2 integers
inline std::size_t
distance(std::size_t _First, std::size_t _Last)
{
  // return int difference
  return _Last - _First;
}

// Bring std::distance() to scope
using std::distance;


/// Helper class: Riemannian graph.
///
/// This class is used internally by mst_orient_normals()
/// to encode:
/// - the adjacency relations of vertices in a K-neighboring.
/// - vertices contain the corresponding input point iterator.
/// - the edge weight = edge weight = 1 - | normal1 * normal2 |
///   where normal1 and normal2 are the normal at the edge extremities.

template <typename ForwardIterator> ///< Input point iterator
struct Riemannian_graph_vertex_properties {
    ForwardIterator input_point; ///< Input point
};
template <typename ForwardIterator> ///< Input point iterator
class Riemannian_graph
  : public boost::adjacency_list< boost::vecS, boost::vecS,
                                  boost::undirectedS,
                                  Riemannian_graph_vertex_properties<ForwardIterator>,
                                  boost::property<boost::edge_weight_t, float> >
{
};


/// Helper class: MST graph
///
/// This class is used internally by mst_orient_normals()
/// to encode:
/// - the adjacency relations of vertices in a Minimum Spanning Tree.
/// - vertices contain the corresponding input point iterator
//    and a boolean indicating if the normal is oriented.

template <typename ForwardIterator> ///< Input point iterator
struct MST_graph_vertex_properties {
    ForwardIterator input_point; ///< Input point
    bool is_oriented; ///< Is input point's normal oriented?
};
template <typename ForwardIterator, ///< Input point iterator
          typename NormalMap, ///< property map: value_type of ForwardIterator -> Normal
          typename Kernel ///< Geometric traits class
>
class MST_graph
  : public boost::adjacency_list< boost::vecS, boost::vecS,
                                  boost::directedS,
                                  MST_graph_vertex_properties<ForwardIterator> >
{
public:
    MST_graph(NormalMap normal_map) : m_normal_map(normal_map) {}

// Public data
    const NormalMap m_normal_map;
};


/// Helper class: Propagate_normal_orientation
///
/// This class is used internally by mst_orient_normals()
/// to propage the normal orientation, starting from a source point
/// and following the adjacency relations of vertices in a Minimum Spanning Tree.
/// It does not orient normals that are already oriented.
/// It does not propagate the orientation if the angle between 2 normals > angle_max.
///
/// \pre Normals must be unit vectors
/// \pre `0 < angle_max <= PI/2`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam NormalMap is a model of `ReadWritePropertyMap`.
/// @tparam Kernel Geometric traits class.

template <typename ForwardIterator, ///< Input point iterator
          typename NormalMap, ///< property map: value_type of ForwardIterator -> Normal
          typename Kernel
>
struct Propagate_normal_orientation
  : public boost::base_visitor< Propagate_normal_orientation<ForwardIterator, NormalMap, Kernel> >
{
    typedef internal::MST_graph<ForwardIterator, NormalMap, Kernel> MST_graph;
    typedef boost::on_examine_edge event_filter;

    Propagate_normal_orientation(double angle_max = CGAL_PI/2.) ///< max angle to propagate the normal orientation (radians)
    : m_angle_max(angle_max)
    {
        // Precondition: 0 < angle_max <= PI/2
        CGAL_point_set_processing_precondition(0 < angle_max && angle_max <= CGAL_PI/2.);
    }

    template <class Edge>
    void operator()(Edge& edge, const MST_graph& mst_graph)
    {
        typedef typename boost::property_traits<NormalMap>::reference Vector_ref;
        typedef typename MST_graph::vertex_descriptor vertex_descriptor;

        // Gets source normal
        vertex_descriptor source_vertex = source(edge, mst_graph);
        Vector_ref source_normal = get(mst_graph.m_normal_map, *(mst_graph[source_vertex].input_point) );
        const bool source_normal_is_oriented = mst_graph[source_vertex].is_oriented;
        // Gets target normal
        vertex_descriptor target_vertex = target(edge, mst_graph);
        Vector_ref target_normal = get( mst_graph.m_normal_map, *(mst_graph[target_vertex].input_point) );
        bool& target_normal_is_oriented = ((MST_graph&)mst_graph)[target_vertex].is_oriented;
        if ( ! target_normal_is_oriented )
        {
          //             ->                        ->
          // Orients target_normal parallel to source_normal
          double normals_dot = source_normal * target_normal;
          if (normals_dot < 0)
          {
            put( mst_graph.m_normal_map, *(mst_graph[target_vertex].input_point), -target_normal );
          }

          // Is orientation robust?
          target_normal_is_oriented
            = source_normal_is_oriented &&
              (std::abs(normals_dot) >= std::cos(m_angle_max)); // oriented iff angle <= m_angle_max
        }
    }

// Data
// Implementation note: boost::breadth_first_search() makes copies of this object => data must be constant or shared.
private:
    const double m_angle_max; ///< max angle to propagate the normal orientation (radians).
};

/// Orients the normal of the point with maximum Z towards +Z axis.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
/// @tparam NormalMap is a model of `ReadWritePropertyMap` with a value_type = Vector_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///
/// @return iterator over the top point.
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename Kernel
>
ForwardIterator
mst_find_source(
    ForwardIterator first,   ///< iterator over the first input point.
    ForwardIterator beyond,  ///< past-the-end iterator over the input points.
    PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3
    NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3
    const Kernel& /*kernel*/)    ///< geometric traits.
{
    CGAL_TRACE("  mst_find_source()\n");

    // Input points types
    typedef typename boost::property_traits<NormalMap>::value_type Vector;
    typedef typename boost::property_traits<NormalMap>::reference Vector_ref;

    // Precondition: at least one element in the container
    CGAL_point_set_processing_precondition(first != beyond);

    // Find top point
    ForwardIterator top_point = first;
    for (ForwardIterator v = ++first; v != beyond; v++)
    {
      
      double top_z = get(point_map,*top_point).z(); // top_point's Z coordinate
      double z = get(point_map,*v).z();
      
      if (top_z < z)
        top_point = v;
    }

    // Orients its normal towards +Z axis
    Vector_ref normal = get(normal_map,*top_point);
    const Vector Z(0, 0, 1);
    if (Z * normal < 0) {
      CGAL_TRACE("  Flip top point normal\n");
    put(normal_map,*top_point, -normal);
    }

    return top_point;
}

/// Iterates over input points and creates Riemannian Graph:
/// - vertices are numbered like the input points index.
/// - vertices contain the corresponding input point iterator.
/// - we add the edge (i, j) if either vertex i is in the k-neighborhood of vertex j,
///   or vertex j is in the k-neighborhood of vertex i.
///
/// \pre Normals must be unit vectors.
/// \pre `k >= 2`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam IndexMap is a model of `ReadablePropertyMap` with an integral value_type.
/// @tparam PointMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
/// @tparam NormalMap is a model of `ReadWritePropertyMap` with a value_type = Vector_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///
/// @return the Riemannian graph
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename IndexMap,
          typename Kernel
>
Riemannian_graph<ForwardIterator>
create_riemannian_graph(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3
    NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3
    IndexMap index_map, ///< property map ForwardIterator -> index
    unsigned int k, ///< number of neighbors
    const Kernel& /*kernel*/) ///< geometric traits.
{
    // Input points types
    typedef typename boost::property_traits<PointMap>::reference Point_ref;
    typedef typename boost::property_traits<NormalMap>::reference Vector_ref;

    // Types for K nearest neighbors search structure
    typedef Point_vertex_handle_3<ForwardIterator> Point_vertex_handle_3;
    typedef Search_traits_vertex_handle_3<ForwardIterator> Traits;
    typedef Euclidean_distance_vertex_handle_3<ForwardIterator> KDistance;
    typedef Orthogonal_k_neighbor_search<Traits,KDistance> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Neighbor_search::iterator Search_iterator;

    // Riemannian_graph types
    typedef internal::Riemannian_graph<ForwardIterator> Riemannian_graph;
    typedef typename boost::property_map<Riemannian_graph, boost::edge_weight_t>::type Riemannian_graph_weight_map;

    // Precondition: at least one element in the container.
    CGAL_point_set_processing_precondition(first != beyond);

    // Precondition: at least 2 nearest neighbors
    CGAL_point_set_processing_precondition(k >= 2);

    // Number of input points
    const std::size_t num_input_points = distance(first, beyond);

    std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Creates KD-tree\n");

    // Instanciate a KD-tree search.
    // Notes: We have to wrap each input point by a Point_vertex_handle_3.
    //        The KD-tree is allocated dynamically to recover RAM as soon as possible.
    std::vector<Point_vertex_handle_3> kd_tree_points; kd_tree_points.reserve(num_input_points);
    for (ForwardIterator it = first; it != beyond; it++)
    {
        
        Point_ref point = get(point_map, *it);
        Point_vertex_handle_3 point_wrapper(point.x(), point.y(), point.z(), it);
        kd_tree_points.push_back(point_wrapper);
    }
    boost::shared_ptr<Tree> tree( new Tree(kd_tree_points.begin(), kd_tree_points.end()) );

    // Recover RAM
    kd_tree_points.clear();

    memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Creates Riemannian Graph\n");

    // Iterates over input points and creates Riemannian Graph:
    // - vertices are numbered like the input points index.
    // - vertices contain the corresponding input point iterator.
    // - we add the edge (i, j) if either vertex i is in the k-neighborhood of vertex j,
    //   or vertex j is in the k-neighborhood of vertex i.
    Riemannian_graph riemannian_graph;
    //
    // add vertices
    for (ForwardIterator it = first; it != beyond; it++)
    {
        typename Riemannian_graph::vertex_descriptor v = add_vertex(riemannian_graph);
        CGAL_point_set_processing_assertion(v == get(index_map,it));
        riemannian_graph[v].input_point = it;
    }
    //
    // add edges
    Riemannian_graph_weight_map riemannian_graph_weight_map = get(boost::edge_weight, riemannian_graph);
    for (ForwardIterator it = first; it != beyond; it++)
    {
        std::size_t it_index = get(index_map,it);
        Vector_ref it_normal_vector = get(normal_map,*it);
        
        // Gather set of (k+1) neighboring points.
        // Perform k+1 queries (as in point set, the query point is
        // output first). Search may be aborted if k is greater
        // than number of input points.
        
        Point_ref point = get(point_map, *it);
        Point_vertex_handle_3 point_wrapper(point.x(), point.y(), point.z(), it);
        Neighbor_search search(*tree, point_wrapper, k+1);
        Search_iterator search_iterator = search.begin();
        for(std::size_t i=0;i<(k+1);i++)
        {
            if(search_iterator == search.end())
                break; // premature ending

            ForwardIterator neighbor = search_iterator->first;
            std::size_t neighbor_index = get(index_map,neighbor);
            if (neighbor_index > it_index) // undirected graph
            {
                // Add edge
                typename boost::graph_traits<Riemannian_graph>::edge_descriptor e;
                bool inserted;
                boost::tie(e, inserted) = add_edge(vertex(it_index, riemannian_graph),
                                                   vertex(neighbor_index, riemannian_graph),
                                                   riemannian_graph);
                CGAL_point_set_processing_assertion(inserted);

                //                               ->        ->
                // Computes edge weight = 1 - | normal1 * normal2 |
                // where normal1 and normal2 are the normal at the edge extremities.
                
                Vector_ref neighbor_normal_vector = get(normal_map,*neighbor);
                double weight = 1.0 - std::abs(it_normal_vector * neighbor_normal_vector);
                if (weight < 0)
                    weight = 0; // safety check
                riemannian_graph_weight_map[e] = (float)weight;
            }

            search_iterator++;
        }
    }

    return riemannian_graph;
}

/// Computes Minimum Spanning Tree and store it in a Boost graph:
/// - vertices are numbered like the input points index.
/// - vertices contain the corresponding input point iterator.
/// - we add the edge (predecessor[i], i) for each element of the MST.
///
/// \pre Normals must be unit vectors.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam IndexMap is a model of `ReadablePropertyMap` with an integral value_type.
/// @tparam PointMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
/// @tparam NormalMap is a model of `ReadWritePropertyMap` with a value_type = Vector_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///
/// @return the MST graph.
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename IndexMap,
          typename Kernel
>
MST_graph<ForwardIterator, NormalMap, Kernel>
create_mst_graph(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3
    NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3
    IndexMap index_map, ///< property map ForwardIterator -> index
    unsigned int k, ///< number of neighbors
    const Kernel& kernel, ///< geometric traits.
    const Riemannian_graph<ForwardIterator>& riemannian_graph, ///< graph connecting each vertex to its knn
    ForwardIterator source_point) ///< source point (with an oriented normal)
{
    // prevents warnings
    CGAL_USE(point_map);
    CGAL_USE(k);
    CGAL_USE(kernel);

    // Bring private stuff to scope
    using namespace internal;

    // Riemannian_graph types
    typedef internal::Riemannian_graph<ForwardIterator> Riemannian_graph;
    typedef typename boost::property_map<Riemannian_graph, boost::edge_weight_t>::const_type Riemannian_graph_weight_map;

    // MST_graph types
    typedef internal::MST_graph<ForwardIterator, NormalMap, Kernel> MST_graph;

    // Precondition: at least one element in the container.
    CGAL_point_set_processing_precondition(first != beyond);

    // Number of input points
    const std::size_t num_input_points = num_vertices(riemannian_graph);

    std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Calls boost::prim_minimum_spanning_tree()\n");

    // Computes Minimum Spanning Tree.
    std::size_t source_point_index = get(index_map, source_point);
    Riemannian_graph_weight_map riemannian_graph_weight_map = get(boost::edge_weight, riemannian_graph);
    typedef std::vector<typename Riemannian_graph::vertex_descriptor> PredecessorMap;
    PredecessorMap predecessor(num_input_points);
    boost::prim_minimum_spanning_tree(riemannian_graph, &predecessor[0],
                                      weight_map( riemannian_graph_weight_map )
                                     .root_vertex( vertex(source_point_index, riemannian_graph) ));

    memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Creates MST Graph\n");

    // Converts predecessor map to a MST graph:
    // - vertices are numbered like the input points index.
    // - vertices contain the corresponding input point iterator.
    // - we add the edge (predecessor[i], i) for each element of the predecessor map.
    MST_graph mst_graph(normal_map);
    //
    // Add vertices. source_point is the unique point marked "oriented".
    for (ForwardIterator it = first; it != beyond; it++)
    {
        // With C++11, the following line triggers a bug in Boost versions
        // 1.56 and 1.57:
        //   https://svn.boost.org/trac/boost/ticket/10382
        typename MST_graph::vertex_descriptor v = add_vertex(mst_graph);
        CGAL_point_set_processing_assertion(v == get(index_map,it));
        mst_graph[v].input_point = it;
        mst_graph[v].is_oriented = (it == source_point);
    }
    // add edges
    for (std::size_t i=0; i < predecessor.size(); i++) // add edges
    {
        if (i != predecessor[i])
        {
            // check that bi-directed graph is useless
            CGAL_point_set_processing_assertion(predecessor[predecessor[i]] != i);

            add_edge(vertex(predecessor[i], mst_graph),
                     vertex(i,     mst_graph),
                     mst_graph);
        }
    }

    return mst_graph;
}

} /* namespace internal */
/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**  
   \ingroup PkgPointSetProcessingAlgorithms
   Orients the normals of the range of `points` using the propagation
   of a seed orientation through a minimum spanning tree of the Riemannian graph.
   This method modifies the order of input points so as to pack all sucessfully oriented points first,
   and returns an iterator over the first point with an unoriented normal (see erase-remove idiom).
   For this reason it should not be called on sorted containers.
   It is based on \cgalCite{cgal:hddms-srup-92}.

   \warning This function may fail when Boost version 1.54 is used,
   because of the following bug: https://svn.boost.org/trac/boost/ticket/9012

   \pre Normals must be unit vectors
   \pre `k >= 2`

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return iterator over the first point with an unoriented normal.
*/
template <typename PointRange,
          typename NamedParameters
>
typename PointRange::iterator
mst_orient_normals(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
    using boost::choose_param;
    CGAL_TRACE("Calls mst_orient_normals()\n");

    typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
    typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
    typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

    CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                                typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                              "Error: no normal map");

    PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
    NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());
    Kernel kernel;

  // Bring private stuff to scope
    using namespace internal;

    // Input points types
    typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Enriched_point; // actual type of input points
    // Property map typename PointRange::iterator -> index
    typedef Index_property_map<typename PointRange::iterator> IndexMap;

    // Riemannian_graph types
    typedef Riemannian_graph<typename PointRange::iterator> Riemannian_graph;

    // MST_graph types
    typedef MST_graph<typename PointRange::iterator, NormalMap, Kernel> MST_graph;

    // Precondition: at least one element in the container.
    CGAL_point_set_processing_precondition(points.begin() != points.end());

    // Precondition: at least 2 nearest neighbors
    CGAL_point_set_processing_precondition(k >= 2);

    std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Create Index_property_map\n");

    // Create a property map Iterator -> index.
    // - if typename PointRange::iterator is a random access iterator (typically vector and deque),
    // get() just calls std::distance() and is very efficient;
    // - else, the property map allocates a std::map to store indices
    // and get() requires a lookup in the map.
    IndexMap index_map(points.begin(), points.end());

    // Orients the normal of the point with maximum Z towards +Z axis.
    typename PointRange::iterator source_point
      = mst_find_source(points.begin(), points.end(),
                        point_map, normal_map,
                        kernel);

    // Iterates over input points and creates Riemannian Graph:
    // - vertices are numbered like the input points index.
    // - vertices are empty.
    // - we add the edge (i, j) if either vertex i is in the k-neighborhood of vertex j,
    //   or vertex j is in the k-neighborhood of vertex i.
    Riemannian_graph riemannian_graph
      = create_riemannian_graph(points.begin(), points.end(),
                                point_map, normal_map, index_map,
                                k,
                                kernel);

    // Creates a Minimum Spanning Tree starting at source_point
    MST_graph mst_graph = create_mst_graph(points.begin(), points.end(),
                                           point_map, normal_map, index_map,
                                           k,
                                           kernel,
                                           riemannian_graph,
                                           source_point);

    memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Calls boost::breadth_first_search()\n");

    // Traverse the point set along the MST to propagate source_point's orientation
    Propagate_normal_orientation<typename PointRange::iterator, NormalMap, Kernel> orienter;
    std::size_t source_point_index = get(index_map, source_point);
    boost::breadth_first_search(mst_graph,
                                vertex(source_point_index, mst_graph), // source
                                visitor(boost::make_bfs_visitor(orienter)));

    // Copy points with robust normal orientation to oriented_points[], the others to unoriented_points[].
    std::deque<Enriched_point> oriented_points, unoriented_points;
    for (typename PointRange::iterator it = points.begin(); it != points.end(); it++)
    {
        std::size_t it_index = get(index_map,it);
        typename MST_graph::vertex_descriptor v = vertex(it_index, mst_graph);
        if (mst_graph[v].is_oriented)
          oriented_points.push_back(*it);
        else
          unoriented_points.push_back(*it);
    }

    // Replaces [points.begin(), points.end()) range by the content of oriented_points[], then unoriented_points[].
    typename PointRange::iterator first_unoriented_point =
      std::copy(oriented_points.begin(), oriented_points.end(), points.begin());
    std::copy(unoriented_points.begin(), unoriented_points.end(), first_unoriented_point);

    // At this stage, we have typically 0 unoriented normals if k is large enough
    CGAL_TRACE("  => %u normals are unoriented\n", unoriented_points.size());

    memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("End of mst_orient_normals()\n");

    return first_unoriented_point;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
typename PointRange::iterator
mst_orient_normals(
  PointRange& points,
  unsigned int k) ///< number of neighbors

{
  return mst_orient_normals (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::mst_orient_normals(), please update your code")
ForwardIterator
mst_orient_normals(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
    NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
    unsigned int k, ///< number of neighbors
    const Kernel& kernel) ///< geometric traits.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return mst_orient_normals
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     geom_traits(kernel));
}
  
// deprecated API
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::mst_orient_normals(), please update your code")
ForwardIterator
mst_orient_normals(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
    NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
    unsigned int k) ///< number of neighbors
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return mst_orient_normals
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template <typename ForwardIterator,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::mst_orient_normals(), please update your code")
ForwardIterator
mst_orient_normals(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
    unsigned int k) ///< number of neighbors
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return mst_orient_normals
    (points,
     k,
     CGAL::parameters::normal_map (normal_map));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MST_ORIENT_NORMALS_H
