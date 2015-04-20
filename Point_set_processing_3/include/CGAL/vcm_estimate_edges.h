// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
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
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_VCM_ESTIMATE_EDGES_H
#define CGAL_VCM_ESTIMATE_EDGES_H

#include <Eigen/Dense>

#include <CGAL/vcm_estimate_normals.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <fstream>
#include <cmath>


namespace CGAL {

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// @cond SKIP_IN_MANUAL
// Determine if a point is on an edge using the VCM of the point.
// A point will be considered as an edge point iff it satisfies a criteria
// relating the eigenvalues of its VCM.
template <class Covariance>
bool
is_on_edge (Covariance &cov,
            double threshold,
            Eigen::Vector3f &dir) {
    // Construct covariance matrix
    Eigen::Matrix3f m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    Eigen::Vector3f eigenvalues;
    Eigen::Matrix3f eigenvectors;
    if (! diagonalize_selfadjoint_matrix(m, eigenvectors, eigenvalues)) {
        return false;
    }

    // Compute the ratio
    float r = eigenvalues(1) / (eigenvalues(0) + eigenvalues(1) + eigenvalues(2));
    if (r >= threshold) {
        dir = eigenvectors.col(1);
        return true;
    }

    return false;
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// Computes the Rips graph of a set of points.
// The Rips graph is a graph where the vertices are the points and
// there is an edge between each vertex contained in a ball of a given radius.
// There is also a cost at each edge.
// The cost between the points pi and pj is given by: c(pi, pj) = norm(p, pi, pj) ^ p
// where p is a given exponent (whose default value is 2 so the cost is the euclidean
// squared distance).
template <typename Undirected_Graph,
          typename K>
void
compute_rips_graph (Undirected_Graph& g, ///< constructed graph.
                    std::vector<typename K::Point_3> points_on_edges, ///< given points.
                    std::map<typename K::Point_3, size_t> indices, ///< map between each point and its corresponding index in the array.
                    double rips_radius, ///< radius of the sphere used to add edges.
                    float exponent, ///< exponent used in the computation of the cost.
                    const K & /* kernel */) ///< geometric traits
{
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::FT FT;

    // KD-Tree
    typedef typename CGAL::Search_traits_3<K> Traits;
    typedef typename CGAL::Kd_tree<Traits> Tree;
    typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
    Tree tree(points_on_edges.begin(), points_on_edges.end());

    // Vertices
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind)
        boost::add_vertex(g);

    // Edges
    // We put edges between two vertices which are in a given ball
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind) {
        std::vector<Point> nn;
        Point p = points_on_edges[ind];
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere(p, rips_radius));

        for (unsigned int k = 0; k < nn.size(); k++) {
            Vector v = nn[k] - p;
            FT cost = pow(fabs(v.x()), exponent) +
                      pow(fabs(v.y()), exponent) +
                      pow(fabs(v.z()), exponent);
            boost::add_edge(ind, indices[nn[k]], cost, g);
        }
    }
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// Computes the nearest neighbors graph of a set of points.
// The nearest neighbors graph is a graph where the vertices are the points and
// there is an edge between each vertex and its k nearest neighbors where k is a given parameter.
// There is also a cost at each edge.
// The cost between the points pi and pj is given by: c(pi, pj) = norm(p, pi, pj) ^ p
// where p is a given exponent (whose default value is 2 so the cost is the euclidean
// squared distance).
template <typename Undirected_Graph,
          typename K>
void
compute_nearest_neighbors_graph (Undirected_Graph& g, ///< constructed graph.
                                 std::vector<typename K::Point_3> points_on_edges, ///< given points.
                                 std::map<typename K::Point_3, size_t> indices, ///< map between each point and its corresponding index in the array.
                                 int nb_neighbors, ///< number of neighbors to consider when constructing the graph.
                                 float exponent, ///< exponent used in the computation of the cost.
                                 const K & /* kernel */) ///< geometric traits
{
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::FT FT;

    typedef typename CGAL::Search_traits_3<K> Traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;

    Tree tree(points_on_edges.begin(), points_on_edges.end());

    // Vertices
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind)
        boost::add_vertex(g);

    // Edges
    // We put edges between each vertex and its k nearest neighbors
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind) {
        std::vector<Point> nn;
        Point p = points_on_edges[ind];
        Neighbor_search search(tree, p, nb_neighbors);
        for (typename Neighbor_search::iterator nit = search.begin();
             nit != search.end();
             ++nit) {
            nn.push_back(nit->first);
        }

        for (unsigned int k = 0; k < nn.size(); k++) {
            Vector v = nn[k] - p;
            FT cost = pow(fabs(v.x()), exponent) +
                      pow(fabs(v.y()), exponent) +
                      pow(fabs(v.z()), exponent);
            boost::add_edge(ind, indices[nn[k]], cost, g);
        }
    }
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// Computes the Delaunay graph of a set of points.
// The Delaunay graph is a graph where the vertices are the points and
// there the edges are the Delaunay edges.
// There is also a cost at each edge.
// The cost between the points pi and pj is given by: c(pi, pj) = norm(p, pi, pj) ^ p
// where p is a given exponent (whose default value is 2 so the cost is the euclidean
// squared distance).
template <typename Undirected_Graph,
          typename K>
void
compute_delaunay_graph (Undirected_Graph& g, ///< constructed graph.
                        std::vector<typename K::Point_3> points_on_edges, ///< given points.
                        std::map<typename K::Point_3, size_t> indices, ///< map between each point and its corresponding index in the array.
                        float exponent, ///< exponent used in the computation of the cost.
                        const K & /* kernel */) ///< geometric traits
{
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::FT FT;

    typedef CGAL::Delaunay_triangulation_3<K> DT;
    typedef typename DT::Finite_edges_iterator Finite_edges_iterator;
    DT dt(points_on_edges.begin(), points_on_edges.end());

    // Vertices
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind)
        boost::add_vertex(g);

    // Edges
    // The edges are the Delaunay edges
    for (Finite_edges_iterator eit = dt.finite_edges_begin();
         eit != dt.finite_edges_end();
         ++eit) {
        typename DT::Segment seg = dt.segment(*eit);
        Point s = seg.source(),
              t = seg.target();
        Vector v = s - t;
        FT cost = pow(fabs(v.x()), exponent) +
                  pow(fabs(v.y()), exponent) +
                  pow(fabs(v.z()), exponent);
        boost::add_edge(indices[s], indices[t], cost, g);
    }
}
/// @endcond

} // namespace internal

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// Estimates the feature edges of the `[first, beyond)` range of points
/// using the Voronoi Covariance Measure.
/// It returns a vector of all the points that have been estimated as edge points.
/// It mainly consists in computing the VCM using `vcm_compute` and then
/// determining which point must be considered as an edge one or not (using a criterion
/// based on the eigenvalues of the covariance matrices).
///
/// See `vcm_compute()` for more details on the VCM.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam Kernel Geometric traits class.
/// @tparam Covariance Covariance matrix type..
template < typename ForwardIterator,
           typename PointPMap,
           typename Kernel
>
std::vector<typename Kernel::Point_3>
vcm_estimate_edges (ForwardIterator first, ///< iterator over the first input point.
                    ForwardIterator beyond, ///< past-the-end iterator over the input points.
                    PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                    double R, ///< offset radius: radius of the sphere to intersect the Voronoi cell with.
                    double r, ///< convolution radius: all points in a sphere with this radius will be convolved.
                    double threshold, ///< threshold used to determine if a point is an edge point or not.
                    const Kernel& k ///< geometric traits.
)
{
    typedef typename Kernel::FT FT;
    typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;

    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    vcm_compute(first, beyond,
                point_pmap,
                cov,
                R,
                r,
                k);

    // Find the potential points on the edges
    std::vector<Point> points_on_edges;
    int i = 0;
    for (ForwardIterator it = first; it != beyond; ++it) {
        Eigen::Vector3f dir;
        if (internal::is_on_edge(cov[i], threshold, dir)) {
            points_on_edges.push_back(get(point_pmap, *it));
        }
        i++;
    }

    return points_on_edges;
}

/// \ingroup PkgPointSetProcessing
/// Constructs a minimum spanning tree where the vertices are the points
/// given in argument.
/// First, it creates a graph where the vertices are the points and there is an edge
/// between each vertices contained in a sphere of a given radius (`rips_radius`).
/// The cost on the edges is the norm on the Lp space to the power p (where p = `exponent`).
/// Secondly, it constructs the MST using BGL and the Kruskal algorithm.
/// It returns a vector of segments where each segment represent a polyline which belongs to the estimated feature.
/// @tparam Kernel Geometric traits class.
template < typename Kernel >
std::vector<typename Kernel::Segment_3>
construct_mst (std::vector<typename Kernel::Point_3> points_on_edges, ///< estimated edges points (see `vcm_estimate_edges()`).
               const Kernel& k, ///< geometric traits.
               double rips_radius = 0.1, ///< radius used for the construction of the Rips graph.
               float exponent = 2 ///< exponent of the cost between edges.
)
{
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Segment_3 Segment;

    // Map between points and their corresponding indices
    std::map<Point, size_t> indices;
    for (size_t s = 0; s < points_on_edges.size(); ++s)
        indices[points_on_edges[s]] = s;

    // Compute the graph
    typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
    typedef boost::adjacency_list<boost::vecS,
                                  boost::vecS,
                                  boost::undirectedS,
                                  boost::no_property,
                                  EdgeWeightProperty > Undirected_Graph;
    typedef boost::graph_traits<Undirected_Graph>::edge_descriptor Edge_descriptor;
    Undirected_Graph g;
    compute_rips_graph(g, points_on_edges, indices, rips_radius, exponent, k);

    // Compute the MST
    boost::property_map<Undirected_Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, g);
    std::vector<Edge_descriptor> spanning_tree;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

    // Construct the polylines
    std::vector<Segment> polylines;
    for (std::vector<Edge_descriptor>::iterator ei = spanning_tree.begin();
         ei != spanning_tree.end();
         ++ei) {
        unsigned int si = boost::source(*ei, g);
        unsigned int ti = boost::target(*ei, g);
        Segment s(points_on_edges[si], points_on_edges[ti]);
        polylines.push_back(s);
    }

    return polylines;
}

} // namespace CGAL

#endif // CGAL_VCM_ESTIMATE_EDGES_H
