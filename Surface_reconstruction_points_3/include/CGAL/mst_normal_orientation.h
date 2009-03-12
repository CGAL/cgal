// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s) : Pierre Alliez and Laurent Saboret and Andreas Fabri

#ifndef CGAL_MST_NORMAL_ORIENTATION_H
#define CGAL_MST_NORMAL_ORIENTATION_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_vertex_handle_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include <iterator>
#include <list>
#include <climits>
#include <math.h>
#ifndef M_PI
  #define M_PI       3.14159265358979323846
#endif

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi { 


/// Generalization of std::distance() to compute the distance between 2 integers
inline int
distance(std::size_t _First, std::size_t _Last)
{
  // return int difference
  return _Last - _First;
}


/// Helper class: Riemannian graph.
///
/// This class is used internally by mst_normal_orientation()
/// to encode:
/// - the adjacency relations of vertices in a K-neighboring
/// - vertices contain the corresponding input vertex handle
/// - the edge weight = edge weight = 1 - | normal1 * normal2 |
///   where normal1 and normal2 are the normal at the edge extremities

template <class VertexIterator>
struct Riemannian_graph_vertex_properties {
    VertexIterator vertex; ///< Input vertex
};
template <class VertexIterator>     ///< Input vertex iterator
class Riemannian_graph
  : public boost::adjacency_list< boost::vecS, boost::vecS,
                                  boost::undirectedS,
                                  Riemannian_graph_vertex_properties<VertexIterator>,
                                  boost::property<boost::edge_weight_t, float> >
{
};


/// Helper class: MST graph
///
/// This class is used internally by mst_normal_orientation()
/// to encode:
/// - the adjacency relations of vertices in a Minimum Spanning Tree
/// - vertices contain the corresponding input vertex handle

template <class VertexIterator>
struct MST_graph_vertex_properties {
    VertexIterator vertex; ///< Input vertex
};
template <class VertexIterator,     ///< Input vertex iterator
          class VertexNormalMap>    ///< property map VertexIterator -> Normal (in and out)
class MST_graph
  : public boost::adjacency_list< boost::vecS, boost::vecS,
                                  boost::directedS,
                                  MST_graph_vertex_properties<VertexIterator> >
{
public:
    MST_graph(VertexNormalMap vertex_normal_map) : m_vertex_normal_map(vertex_normal_map) {}

// Public data
    const VertexNormalMap m_vertex_normal_map;
};


/// Helper class: Propagate_normal_orientation
///
/// This class is used internally by mst_normal_orientation()
/// to propage the normal orientation, starting from a source vertex
/// and following the adjacency relations of vertices in a Minimum Spanning Tree.
///
/// This variant:
/// - does not orient normals that are already oriented.
/// - does not propagate the orientation if the angle between 2 normals > angle_max.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - Normals must be unit vectors.
/// - 0 < angle_max <= PI/2.

template <class VertexIterator,     ///< Input vertex iterator
          class VertexNormalMap>    ///< property map VertexIterator -> Normal (in and out)
struct Propagate_normal_orientation
  : public boost::base_visitor< Propagate_normal_orientation<VertexIterator, VertexNormalMap> >
{
    typedef CGALi::MST_graph<VertexIterator, VertexNormalMap> MST_graph;
    typedef boost::on_examine_edge event_filter;

    Propagate_normal_orientation(double angle_max) ///< max angle to propagate the normal orientation (radians)
    : m_angle_max(angle_max)
    {
        // Precondition: 0 < angle_max <= PI/2
        CGAL_surface_reconstruction_precondition(0 < angle_max && angle_max <= M_PI/2.);
    }

    template <class Edge>
    void operator()(const Edge& edge, const MST_graph& graph)
    {
        typedef typename boost::property_traits<VertexNormalMap>::value_type Normal;
        typedef typename Normal::Vector Vector;
        typedef typename MST_graph::vertex_descriptor vertex_descriptor;

        // Get source normal
        vertex_descriptor source_vertex = boost::source(edge, graph);
        Normal& source_normal = graph.m_vertex_normal_map[graph[source_vertex].vertex];
        Vector source_vector = source_normal;

        // Get target normal
        vertex_descriptor target_vertex = boost::target(edge, graph);
        Normal& target_normal = graph.m_vertex_normal_map[graph[target_vertex].vertex];
        Vector target_vector = target_normal;

        if ( ! target_normal.is_oriented() )
        {
          //             ->                        ->
          // Orient target_normal parallel to source_normal
          double normals_dot = source_vector * target_vector;
          if (normals_dot < 0) {
            //CGAL_TRACE("    flip %d\n", (int)target_vertex);
            target_vector = -target_vector;
          }

          // Is orientation robust?
          bool oriented = source_normal.is_oriented() &&
                          (std::abs(normals_dot) >= std::cos(m_angle_max)); // oriented iff angle <= m_angle_max
          target_normal = Normal(target_vector, oriented);
        }
    }

// Data
// Implementation note: boost::breadth_first_search() makes copies of this object => data must be constant or shared.
private:
    const double m_angle_max; ///< max angle to propagate the normal orientation (radians).
};

/// If the point set contains points with an oriented normal, find one.
/// Else, orient the normal of the vertex with maximum Z towards +Z axis.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexIndexMap is a model of boost::readable_property_map.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
///
/// @return the top point.

template<class VertexIterator, class VertexPointMap, class VertexIndexMap, class VertexNormalMap>
VertexIterator
find_source_mst_3(
    VertexIterator first, ///< first input vertex
    VertexIterator beyond, ///< past-the-end input vertex
    VertexIndexMap vertex_index_map, ///< property map VertexIterator -> index
    VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
    VertexNormalMap vertex_normal_map) ///< property map VertexIterator -> Normal (in and out)
{
    CGAL_TRACE("  orient_highest_point_normal_3()\n");

    // Input mesh's types
    typedef typename boost::property_traits<VertexPointMap>::value_type Point;
    typedef typename boost::property_traits<VertexNormalMap>::value_type Normal;
    typedef typename Normal::Vector Vector;

    // Precondition: at least one element in the container.
    CGAL_surface_reconstruction_precondition(first != beyond);

    // If the point set contains points with an oriented normal, find top one.
    // Else, find top point.
    //
    // Invariant: among traversed vertices, top_vertex is
    //            the top vertex with an oriented normal if at least one was traversed,
    //            else the top vertex.
    VertexIterator top_vertex = first;
    double top_z = get(vertex_point_map,first).z(); // top_vertex's Z coordinate
    bool top_normal_is_oriented = vertex_normal_map[top_vertex].is_oriented();
    for (VertexIterator v = ++first; v != beyond; v++)
    {
      double z = get(vertex_point_map,v).z();
      bool normal_is_oriented = vertex_normal_map[v].is_oriented();
      if (top_normal_is_oriented)
      {
        if (normal_is_oriented && top_z < z) {
            top_vertex = v;
            top_z = z;
            top_normal_is_oriented = normal_is_oriented;
        }
      }
      else
      {
        if (normal_is_oriented || top_z < z) {
            top_vertex = v;
            top_z = z;
            top_normal_is_oriented = normal_is_oriented;
        }
      }
    }
    //CGAL_TRACE_STREAM << "  Top vertex index = " << get(vertex_index_map,top_vertex) << std::endl;

    // Orient its normal towards +Z axis
    Normal& normal = vertex_normal_map[top_vertex];
    if ( ! normal.is_oriented() )
    {
      const Vector Z(0, 0, 1);
      Vector vec = normal;
      if (Z * vec < 0) {
          CGAL_TRACE("  Flip top vertex normal\n");
          vec = -vec;
      }
      normal = Normal(vec, true /* oriented */);
    }

    return top_vertex;
}

/// Iterate over input points and create Riemannian Graph:
/// - vertices are numbered like the input vertices' index.
/// - vertices are empty.
/// - we add the edge (i, j) if either vertex i is in the k-neighborhood of vertex j,
///   or vertex j is in the k-neighborhood of vertex i.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexIndexMap is a model of boost::readable_property_map.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - Normals must be unit vectors.
/// - k >= 2.

template<class VertexIterator, class VertexPointMap, class VertexIndexMap, class VertexNormalMap>
Riemannian_graph<VertexIterator>
create_riemannian_graph(
    VertexIterator first, ///< first input vertex
    VertexIterator beyond, ///< past-the-end input vertex
    VertexIndexMap vertex_index_map, ///< property map VertexIterator -> index
    VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
    VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
    unsigned int k) ///< number of neighbors
{
    // Input mesh's types
    typedef typename boost::property_traits<VertexPointMap>::value_type Point;
    typedef typename boost::property_traits<VertexNormalMap>::value_type Normal;
    typedef typename Normal::Vector Vector;

    // Types for K nearest neighbors search structure
    typedef Point_vertex_handle_3<VertexIterator> Point_vertex_handle_3;
    typedef Search_traits_vertex_handle_3<VertexIterator> Traits;
    typedef Euclidean_distance_vertex_handle_3<VertexIterator> KDistance;
    typedef Orthogonal_k_neighbor_search<Traits,KDistance> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Neighbor_search::iterator Search_iterator;

    // Riemannian_graph types
    typedef CGALi::Riemannian_graph<VertexIterator> Riemannian_graph;
    typedef typename boost::property_map<Riemannian_graph, boost::edge_weight_t>::type Riemannian_graph_weight_map;

    // Precondition: at least one element in the container.
    CGAL_surface_reconstruction_precondition(first != beyond);

    // Precondition: at least 2 nearest neighbors
    CGAL_surface_reconstruction_precondition(k >= 2);

    // Number of input vertices
    const int num_input_vertices = distance(first, beyond);

    long memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Create KD-tree\n");

    // Instanciate a KD-tree search.
    // Notes: We have to wrap each input vertex by a Point_vertex_handle_3.
    //        The KD-tree is allocated dynamically to recover RAM as soon as possible.
    std::vector<Point_vertex_handle_3> kd_tree_points; kd_tree_points.reserve(num_input_vertices);
    for (VertexIterator it = first; it != beyond; it++)
    {
        Point point = get(vertex_point_map, it);
        Point_vertex_handle_3 point_wrapper(point.x(), point.y(), point.z(), it);
        kd_tree_points.push_back(point_wrapper);
    }
    std::auto_ptr<Tree> tree( new Tree(kd_tree_points.begin(), kd_tree_points.end()) );

    // Recover RAM
    kd_tree_points.clear();

    /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Create Riemannian Graph\n");

    // Iterate over input points and create Riemannian Graph:
    // - vertices are numbered like the input vertices' index.
    // - vertices contain the corresponding input vertex handle.
    // - we add the edge (i, j) if either vertex i is in the k-neighborhood of vertex j,
    //   or vertex j is in the k-neighborhood of vertex i.
    Riemannian_graph riemannian_graph;
    //
    // add vertices
    for (VertexIterator it = first; it != beyond; it++)
    {
        typename Riemannian_graph::vertex_descriptor v = add_vertex(riemannian_graph);
        CGAL_surface_reconstruction_assertion(v == get(vertex_index_map,it));
        riemannian_graph[v].vertex = it;
    }
    //
    // add edges
    Riemannian_graph_weight_map riemannian_graph_weight_map = get(boost::edge_weight, riemannian_graph);
    for (VertexIterator it = first; it != beyond; it++)
    {
        unsigned int it_index = get(vertex_index_map,it);
        Vector it_normal_vector = vertex_normal_map[it];

        // Gather set of (k+1) neighboring points.
        // Perform k+1 queries (as in point set, the query point is
        // output first). Search may be aborted when k is greater
        // than number of input points.
        Point point = get(vertex_point_map, it);
        Point_vertex_handle_3 point_wrapper(point.x(), point.y(), point.z(), it);
        Neighbor_search search(*tree, point_wrapper, k+1);
        Search_iterator search_iterator = search.begin();
        for(unsigned int i=0;i<(k+1);i++)
        {
            if(search_iterator == search.end())
                break; // premature ending

            VertexIterator neighbor = search_iterator->first;
            unsigned int neighbor_index = get(vertex_index_map,neighbor);
            if (neighbor_index > it_index) // undirected graph
            {
                // Add edge
                typename boost::graph_traits<Riemannian_graph>::edge_descriptor e;
                bool inserted;
                boost::tie(e, inserted) = boost::add_edge(boost::vertex(it_index, riemannian_graph),
                                                          boost::vertex(neighbor_index, riemannian_graph),
                                                          riemannian_graph);
                CGAL_surface_reconstruction_assertion(inserted);

                //                               ->        ->
                // Compute edge weight = 1 - | normal1 * normal2 |
                // where normal1 and normal2 are the normal at the edge extremities.
                Vector neighbor_normal_vector = vertex_normal_map[neighbor];
                double weight = 1.0 - std::abs(it_normal_vector * neighbor_normal_vector);
                if (weight < 0)
                    weight = 0; // safety check
                //CGAL_TRACE("    %d (%1.3lf,%1.3lf,%1.3lf) -> %d (%1.3lf,%1.3lf,%1.3lf): weight=%1.3lf\n",
                //           (int)it_index, it_normal_vector.x(),it_normal_vector.y(),it_normal_vector.z(),
                //           (int)neighbor_index, neighbor_normal_vector.x(),neighbor_normal_vector.y(),neighbor_normal_vector.z(),
                //           weight);
                riemannian_graph_weight_map[e] = (float)weight;
            }

            search_iterator++;
        }
    }

    return riemannian_graph;
}

/// Compute Minimum Spanning Tree and store it in a Boost graph:
/// - vertices are numbered like the input vertices' index.
/// - vertices contain the corresponding input vertex handle.
/// - we add the edge (predecessor[i], i) for each element of the MST.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexIndexMap is a model of boost::readable_property_map.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - Normals must be unit vectors.
///
/// @return the number of un-oriented normals.

template<class VertexIterator, class VertexPointMap, class VertexIndexMap, class VertexNormalMap>
MST_graph<VertexIterator, VertexNormalMap>
create_mst_graph(
    VertexIterator first, ///< first input vertex
    VertexIterator beyond, ///< past-the-end input vertex
    VertexIndexMap vertex_index_map, ///< property map VertexIterator -> index
    VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
    VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
    const Riemannian_graph<VertexIterator>& riemannian_graph, ///< graph connecting each vertex to its k
    VertexIterator source_vertex) ///< source vertex (with an oriented normal)
{
    // Bring private stuff to scope
    using namespace CGALi;

    // Input mesh's types
    typedef typename boost::property_traits<VertexPointMap>::value_type Point;
    typedef typename boost::property_traits<VertexNormalMap>::value_type Normal;
    typedef typename Normal::Vector Vector;

    // Riemannian_graph types
    typedef CGALi::Riemannian_graph<VertexIterator> Riemannian_graph;
    typedef typename boost::property_map<Riemannian_graph, boost::edge_weight_t>::const_type Riemannian_graph_weight_map;

    // MST_graph types
    typedef CGALi::MST_graph<VertexIterator, VertexNormalMap> MST_graph;

    // Precondition: at least one element in the container.
    CGAL_surface_reconstruction_precondition(first != beyond);

    // Precondition: the source vertex's normal must be oriented.
    CGAL_surface_reconstruction_precondition(vertex_normal_map[source_vertex].is_oriented());

    // Number of input vertices
    const int num_input_vertices = boost::num_vertices(riemannian_graph);

    long memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Call boost::prim_minimum_spanning_tree()\n");

    // Compute Minimum Spanning Tree.
    unsigned int source_vertex_index = get(vertex_index_map, source_vertex);
    Riemannian_graph_weight_map riemannian_graph_weight_map = get(boost::edge_weight, riemannian_graph);
    typedef std::vector<typename Riemannian_graph::vertex_descriptor> PredecessorMap;
    PredecessorMap predecessor(num_input_vertices);
    boost::prim_minimum_spanning_tree(riemannian_graph, &predecessor[0],
                                      weight_map( riemannian_graph_weight_map )
                                     .root_vertex( boost::vertex(source_vertex_index, riemannian_graph) ));

    /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Create MST Graph\n");

    // Convert predecessor map to a MST graph:
    // - vertices are numbered like the input vertices' index.
    // - vertices contain the corresponding input vertex handle.
    // - we add the edge (predecessor[i], i) for each element of the predecessor map.
    MST_graph mst_graph(vertex_normal_map);
    //
    // add vertices
    for (VertexIterator it = first; it != beyond; it++)
    {
        typename MST_graph::vertex_descriptor v = add_vertex(mst_graph);
        CGAL_surface_reconstruction_assertion(v == get(vertex_index_map,it));
        mst_graph[v].vertex = it;
    }
    // add edges
    for (unsigned int i=0; i < predecessor.size(); i++) // add edges
    {
        if (i != predecessor[i])
        {
            // check that bi-directed graph is useless
            CGAL_surface_reconstruction_assertion(predecessor[predecessor[i]] != i);

            boost::add_edge(boost::vertex(predecessor[i], mst_graph),
                            boost::vertex(i,     mst_graph),
                            mst_graph);
        }
    }

    return mst_graph;
}


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Orient the normals of a point set using the method described by
/// Hoppe, DeRose, Duchamp, McDonald and Stuetzle in
/// "Surface reconstruction from unorganized points" [Hoppe92].
///
/// This variant implements the original algorithm.
/// Note that it does not orient normals that are already oriented.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexIndexMap is a model of boost::readable_property_map.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - Normals must be unit vectors.
/// - k >= 2.
///
/// @return the number of un-oriented normals.

template<class VertexIterator, class VertexPointMap, class VertexIndexMap, class VertexNormalMap>
unsigned int
mst_normal_orientation(
    VertexIterator first, ///< first input vertex
    VertexIterator beyond, ///< past-the-end input vertex
    VertexIndexMap vertex_index_map, ///< property map VertexIterator -> index
    VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
    VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
    unsigned int k) ///< number of neighbors
{
    return mst_normal_orientation(first, beyond,
                                  vertex_index_map, vertex_point_map, vertex_normal_map,
                                  k,
                                  M_PI/2.); // always propagate normal orientation
}

/// Orient the normals of a point set using the method described by
/// Hoppe et al. in
/// "Surface reconstruction from unorganized points" [Hoppe92].
///
/// This is a variant of the original algorithm. It:
/// - orients the top point towards +Z axis.
/// - does not orient normals that are already oriented.
/// - does not propagate the orientation if the angle between 2 normals > angle_max.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexIndexMap is a model of boost::readable_property_map.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - Normals must be unit vectors.
/// - k >= 2.
/// - 0 < angle_max <= PI/2.
///
/// @return the number of un-oriented normals.

template<class VertexIterator, class VertexPointMap, class VertexIndexMap, class VertexNormalMap>
unsigned int
mst_normal_orientation(
    VertexIterator first, ///< first input vertex
    VertexIterator beyond, ///< past-the-end input vertex
    VertexIndexMap vertex_index_map, ///< property map VertexIterator -> index
    VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
    VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
    unsigned int k, ///< number of neighbors
    double angle_max) ///< max angle to propagate the normal orientation (radians)
{
    CGAL_TRACE("Call mst_normal_orientation(angle_max=%lf degrees)\n", angle_max*180.0/M_PI);

    // Bring private stuff to scope
    using namespace CGALi;

    // Input mesh's types
    typedef typename boost::property_traits<VertexPointMap>::value_type Point;
    typedef typename boost::property_traits<VertexNormalMap>::value_type Normal;
    typedef typename Normal::Vector Vector;

    // Riemannian_graph types
    typedef Riemannian_graph<VertexIterator> Riemannian_graph;

    // MST_graph types
    typedef MST_graph<VertexIterator, VertexNormalMap> MST_graph;

    // Precondition: at least one element in the container.
    CGAL_surface_reconstruction_precondition(first != beyond);

    // Precondition: at least 2 nearest neighbors
    CGAL_surface_reconstruction_precondition(k >= 2);

    // Precondition: 0 < angle_max <= PI/2
    CGAL_surface_reconstruction_precondition(0 < angle_max && angle_max <= M_PI/2.);

    // If the point set contains points with an oriented normal, find one.
    // Else, orient the normal of the vertex with maximum Z towards +Z axis.
    VertexIterator source_vertex
      = find_source_mst_3(first, beyond,
                          vertex_index_map, vertex_point_map, vertex_normal_map);

    // Iterate over input points and create Riemannian Graph:
    // - vertices are numbered like the input vertices' index.
    // - vertices are empty.
    // - we add the edge (i, j) if either vertex i is in the k-neighborhood of vertex j,
    //   or vertex j is in the k-neighborhood of vertex i.
    Riemannian_graph riemannian_graph
      = create_riemannian_graph(first, beyond,
                                vertex_index_map, vertex_point_map, vertex_normal_map,
                                k);

    // Create a Minimum Spanning Tree starting at source_vertex
    MST_graph mst_graph = create_mst_graph(first, beyond,
                                           vertex_index_map, vertex_point_map, vertex_normal_map,
                                           riemannian_graph,
                                           source_vertex);

    long memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("  Call boost::breadth_first_search()\n");

    // Traverse the point set along the MST to propagate source_vertex's orientation
    Propagate_normal_orientation<VertexIterator, VertexNormalMap> orienter(angle_max);
    unsigned int source_vertex_index = get(vertex_index_map, source_vertex);
    boost::breadth_first_search(mst_graph,
                                boost::vertex(source_vertex_index, mst_graph), // source
                                visitor(boost::make_bfs_visitor(orienter)));

    // Count un-oriented normals
    unsigned int unoriented_normals = 0;
    for (VertexIterator it = first; it != beyond; it++)
        if ( ! vertex_normal_map[it].is_oriented() )
          unoriented_normals++;
    CGAL_TRACE("  => %u normals are unoriented\n", unoriented_normals);

    // At this stage, we have typically:
    // - 0 unoriented normals if angle_max = PI/2 and k is large enough
    // - <1% of unoriented normals if angle_max = PI/4

#if 0
    // Enhanced version of the algorithm: 2nd pass
    int pass = 1;
    unsigned int old_unoriented_normals = -1;
    while (unoriented_normals != 0 && unoriented_normals != old_unoriented_normals && angle_max < M_PI/2.)
    {
      /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
      CGAL_TRACE_STREAM << "  pass " << ++pass << "\n";

      old_unoriented_normals = unoriented_normals;

      // For each unoriented normal
      for (VertexIterator target_it = first; (target_it != beyond) && (unoriented_normals != 0); target_it++)
      {
          Normal& target_normal = vertex_normal_map[target_it];
          Vector target_vector = target_normal;
          if ( ! target_normal.is_oriented() )
          {
              typedef typename Riemannian_graph::vertex_descriptor vertex_descriptor;

              // Convert target_it to a riemannian_graph vertex target_vertex
              unsigned int index = get(vertex_index_map, target_it);
              vertex_descriptor target_vertex = boost::vertex(index, riemannian_graph);

              // Find neighbor source_vertex which has an oriented normal
              // and such that the normals angle is minimum.
              vertex_descriptor source_vertex = -1;
              double normals_dot_max = 0;
              boost::graph_traits<Riemannian_graph>::out_edge_iterator e, e_end;
              for (boost::tie(e, e_end) = boost::out_edges(target_vertex, riemannian_graph); e != e_end; ++e)
              {
                // Get neighbor's normal
                vertex_descriptor s = boost::target(*e, riemannian_graph);
                Normal& n = vertex_normal_map[riemannian_graph[s].vertex];
                Vector v = n;
                double normals_dot = v * target_vector;
                if ( n.is_oriented() && (std::abs(normals_dot) > normals_dot_max) )
                {
                  source_vertex = s;
                  normals_dot_max = std::abs(normals_dot);
                }
              }

              if (source_vertex != -1)
              {
                //             ->                        ->
                // Orient target_normal parallel to source_normal
                Normal& source_normal = vertex_normal_map[riemannian_graph[source_vertex].vertex];
                Vector source_vector = source_normal;
                double normals_dot = source_vector * target_vector;
                if (normals_dot < 0) {
                  //CGAL_TRACE("    flip %d\n", (int)target_vertex);
                  target_vector = -target_vector;
                }

                // Is orientation robust?
                CGAL_surface_reconstruction_assertion(source_normal.is_oriented());
                //bool oriented = (std::abs(normals_dot) >= std::cos(angle_max)); // oriented iff angle <= m_angle_max
                bool oriented = true;
                target_normal = Normal(target_vector, oriented);

                // Update the number of unoriented normals
                if (oriented)
                {
                  //CGAL_TRACE("    orient %d\n", (int)target_vertex);
                  CGAL_surface_reconstruction_assertion(unoriented_normals > 0);
                  unoriented_normals--;
                }
              }
          }
      }
      CGAL_TRACE("  => %u normals are unoriented\n", unoriented_normals);
    }
#endif // 0

#if 0
    // Enhanced version of the algorithm: 2nd pass
    if (unoriented_normals != 0 && angle_max < M_PI/2.)
    {
      long memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
      CGAL_TRACE("  2nd pass\n");
      CGAL_TRACE("  Call boost::breadth_first_search()\n");

      // Traverse the point set along the MST to propagate source_vertex's orientation
      Propagate_normal_orientation<VertexIterator, VertexNormalMap> orienter(M_PI/2.);
      unsigned int source_vertex_index = get(vertex_index_map, source_vertex);
      boost::breadth_first_search(mst_graph,
                                  boost::vertex(source_vertex_index, mst_graph), // source
                                  visitor(boost::make_bfs_visitor(orienter)));

      // Count un-oriented normals
      unoriented_normals = 0;
      for (VertexIterator it = first; it != beyond; it++)
          if ( ! vertex_normal_map[it].is_oriented() )
            unoriented_normals++;
      CGAL_TRACE("  => %u normals are unoriented\n", unoriented_normals);
    }
#endif // 0

    /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
    CGAL_TRACE("End of mst_normal_orientation()\n");

    return unoriented_normals;
}


CGAL_END_NAMESPACE

#endif // CGAL_MST_NORMAL_ORIENTATION_H

