#ifndef CGAL_VCM_ESTIMATE_EDGES_H
#define CGAL_VCM_ESTIMATE_EDGES_H

#include <Eigen/Dense>

#include <CGAL/vcm_utilities.h>
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

// TODO:
// - remove traces
// - use average spacing to determine edge_radius

namespace CGAL {

namespace internal {

// Determine if a point is on an edge
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
    // TODO: remove trace
    /* std::cout << r << std::endl; */
    if (r >= threshold) {
        dir = eigenvectors.col(1);
        return true;
    }

    return false;
}

template <typename Undirected_Graph,
          typename K>
void
compute_rips_graph (Undirected_Graph& g,
                    std::vector<typename K::Point_3> points_on_edges,
                    std::map<typename K::Point_3, size_t> indices,
                    double rips_radius,
                    float exponent,
                    const K & /* kernel */) {
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

template <typename Undirected_Graph,
          typename K>
void
compute_nearest_neighbors_graph (Undirected_Graph& g,
                                 std::vector<typename K::Point_3> points_on_edges,
                                 std::map<typename K::Point_3, size_t> indices,
                                 int nb_neighbors,
                                 float exponent,
                                 const K & /* kernel */) {
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

template <typename Undirected_Graph,
          typename K>
void
compute_delaunay_graph (Undirected_Graph& g,
                        std::vector<typename K::Point_3> points_on_edges,
                        std::map<typename K::Point_3, size_t> indices,
                        float exponent,
                        const K & /* kernel */) {
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

} // namespace internal

// Estimate sharp edges using VCM
template < typename ForwardIterator,
           typename PointPMap,
           typename Kernel
>
std::vector<typename Kernel::Segment_3>
vcm_estimate_edges (ForwardIterator first,
                    ForwardIterator beyond,
                    PointPMap point_pmap,
                    double R,
                    double r,
                    double threshold,
                    const Kernel& k,
                    double rips_radius = 0.1,
                    float exponent = 2)
{
    typedef typename Kernel::FT FT;
    typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Segment_3 Segment;
    typedef typename Kernel::FT FT;

    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    internal::vcm_offset_and_convolve(first, beyond,
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

    // Map between points and their corresponding indices
    std::map<Point, size_t> indices;
    for (size_t s = 0; s < points_on_edges.size(); ++s)
        indices[points_on_edges[s]] = s;

    // TODO: debug
    std::ofstream file_edges("edges.xyz");
    for (typename std::vector<Point>::iterator it = points_on_edges.begin();
         it != points_on_edges.end();
         ++it) {
        file_edges << *it << "\n";
    }

    // Compute the graph
    typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
    typedef boost::adjacency_list<boost::vecS,
                                  boost::vecS,
                                  boost::undirectedS,
                                  boost::no_property,
                                  EdgeWeightProperty > Undirected_Graph;
    typedef typename boost::graph_traits<Undirected_Graph>::vertex_iterator Vertex_iterator;
    typedef boost::graph_traits<Undirected_Graph>::vertex_descriptor Vertex_descriptor;
    typedef boost::graph_traits<Undirected_Graph>::edge_descriptor Edge_descriptor;
    typedef boost::graph_traits<Undirected_Graph>::edge_iterator Edge_iterator;
    Undirected_Graph g;
    compute_rips_graph(g, points_on_edges, indices, rips_radius, exponent, k);
    /* compute_nearest_neighbors_graph(g, points_on_edges, indices, 50, exponent, k); */
    /* compute_delaunay_graph(g, points_on_edges, indices, exponent, k); */

    // Compute the MST
    boost::property_map<Undirected_Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, g);
    std::vector<Edge_descriptor> spanning_tree;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

    // Compute the distance function
    // We first create a new graph containing the tree
    Undirected_Graph mst;
    Vertex_descriptor v0 = boost::add_vertex(mst);
    for (unsigned int ind = 1; ind < points_on_edges.size(); ++ind)
        boost::add_vertex(mst);

    for (std::vector<Edge_descriptor>::iterator ei = spanning_tree.begin();
         ei != spanning_tree.end();
         ++ei) {
        unsigned int si = boost::source(*ei, g);
        unsigned int ti = boost::target(*ei, g);
        boost::add_edge(si, ti, weight[*ei], mst);
    }

    // Then, we compute the shortest paths using Dijkstra algorithm
    std::vector<Vertex_descriptor> parents(boost::num_vertices(mst));
    std::vector<float> distances(boost::num_vertices(mst));
    boost::dijkstra_shortest_paths(mst, v0, boost::predecessor_map(&parents[0]).distance_map(&distances[0]));

    // TODO: remove trace
    /* std::cout << "distances and parents:" << std::endl; */
    std::ofstream file_distances("distances.xyz");
    Vertex_iterator vertexIterator, vend;
    for (boost::tie(vertexIterator, vend) = boost::vertices(g); vertexIterator != vend; ++vertexIterator) {
        /* std::cout << "distance(" << *vertexIterator << ") = " << distances[*vertexIterator] << ", "; */
        /* std::cout << "parent(" << *vertexIterator << ") = " << parents[*vertexIterator] << std::endl; */
        file_distances << distances[*vertexIterator] << "\n";
    }
    /* std::cout << std::endl; */

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
