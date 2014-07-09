#ifndef CGAL_VCM_ESTIMATE_EDGES_H
#define CGAL_VCM_ESTIMATE_EDGES_H

#include <eigen3/Eigen/Dense>
#include <CGAL/vcm_utilities.h>
#include <CGAL/vcm_estimate_normals.h>

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
    // TODO
    float r = eigenvalues(0) / (eigenvalues(0) + eigenvalues(1) + eigenvalues(2));
    std::cout << r << std::endl;
    if (r <= threshold) {
        dir = eigenvectors.col(1);
        return true;
    }

    return false;
}

} // namespace internal

// Estimate sharp edges using VCM
template < typename ForwardIterator,
           typename PointPMap,
           typename Kernel,
           typename Covariance
>
void
vcm_estimate_edges (ForwardIterator first,
                    ForwardIterator beyond,
                    PointPMap point_pmap,
                    double R,
                    double r,
                    double threshold,
                    const Kernel &,
                    const Covariance &)
{
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;

    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    internal::vcm_offset_and_convolve(first, beyond,
                                      point_pmap,
                                      cov,
                                      R,
                                      r,
                                      Kernel());

    // Find the potential points on the edges
    std::map<Point, Vector> points_on_edges;
    int i = 0;
    for (ForwardIterator it = first; it != beyond; ++it) {
        Eigen::Vector3f dir;
        if (internal::is_on_edge(cov[i], threshold, dir)) {
            Vector cdir(dir[0], dir[1], dir[2]);
            points_on_edges[get(point_pmap, *it)] = cdir;
        }
        i++;
    }

    // Compute the graph
    // TODO

    // Compute the MST
    // TODO

    // Compute the Polylines
    // TODO

    // TODO: debug
    /* std::ofstream file_edges("edges.xyz"); */
    /* std::copy(points_on_edges.begin(), points_on_edges.end(), */
    /*           std::ostream_iterator<Covariance>(file_edges)); */
}

} // namespace CGAL

#endif // CGAL_VCM_ESTIMATE_EDGES_H
