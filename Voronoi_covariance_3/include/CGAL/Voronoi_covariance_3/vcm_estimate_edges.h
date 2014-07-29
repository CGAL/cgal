#ifndef CGAL_VCM_ESTIMATE_EDGES_H
#define CGAL_VCM_ESTIMATE_EDGES_H

#include <Eigen/Dense>

#include <CGAL/vcm_utilities.h>
#include <CGAL/vcm_estimate_normals.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <fstream>

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
    std::cout << r << std::endl;
    if (r >= threshold) {
        dir = eigenvectors.col(1);
        return true;
    }

    return false;
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
                    const Kernel& k)
{
    typedef typename Kernel::FT FT;
    typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Segment_3 Segment;

    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    internal::vcm_offset_and_convolve(first, beyond,
                                      point_pmap,
                                      cov,
                                      R,
                                      r,
                                      k);

    // Find the potential points on the edges
    // TODO: direction?
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
                                  EdgeWeightProperty > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex_descriptor;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge_descriptor;
    Graph g;

    // KD-Tree
    typedef CGAL::Search_traits_3<Kernel> Traits;
    typedef CGAL::Kd_tree<Traits> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
    Tree tree(points_on_edges.begin(), points_on_edges.end());

    // Vertices
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind)
        Vertex_descriptor v = boost::add_vertex(g);

    // Edges
    // We put edges between two vertices which are in a given ball
    // The cost is the squared distance
    // TODO: radius of sphere
    std::cout << "tree:" << std::endl;
    for (unsigned int ind = 0; ind < points_on_edges.size(); ++ind) {
        std::vector<Point> nn;
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere (points_on_edges[ind], 0.1));
        std::cout << nn.size() << std::endl;
        for (unsigned int k = 0; k < nn.size(); k++) {
            Vector v = nn[k] - points_on_edges[ind];
            boost::add_edge(ind, indices[nn[k]], v.squared_length(), g);
        }
    }

    // Compute the MST
    std::cout << "mst:" << std::endl;
    boost::property_map<Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, g);
    std::vector<Edge_descriptor> spanning_tree;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

    // TODO: remove trace
    for (std::vector<Edge_descriptor>::iterator ei = spanning_tree.begin();
         ei != spanning_tree.end();
         ++ei) {
        std::cout << boost::source(*ei, g) << " <--> " << boost::target(*ei, g)
            << " with weight of " << weight[*ei]
            << std::endl;
    }

    // Compute the polylines as segments
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
