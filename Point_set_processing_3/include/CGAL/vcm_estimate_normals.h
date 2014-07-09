#ifndef CGAL_VCM_ESTIMATE_NORMALS_H
#define CGAL_VCM_ESTIMATE_NORMALS_H

#include <CGAL/property_map.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/voronoi_covariance_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <iterator>
#include <vector>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// \cond SKIP_IN_MANUAL

/// Computes the VCM for each point in the property map.
/// The matrix is computed by intersecting the Voronoi cell
/// of a point and a sphere whose radius is `R` and discretized
/// by `N` planes.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
/// @tparam K Geometric traits class.
/// @tparam Covariance Covariance matrix type.
template < typename ForwardIterator,
           typename PointPMap,
           class K,
           class Covariance
>
void
vcm_offset (ForwardIterator first, ///< iterator over the first input point.
            ForwardIterator beyond, ///< past-the-end iterator over the input points.
            PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
            std::vector<Covariance> &cov, ///< vector of covariance matrices.
            double R, ///< radius of the sphere.
            size_t N, ///< number of planes used to discretize the sphere.
            const K & /*kernel*/) ///< geometric traits.
{
    // Sphere discretization
    typename CGAL::Voronoi_covariance_3::Sphere_discretization<K> sphere(R, N);

    // Compute the Delaunay Triangulation
    typedef CGAL::Delaunay_triangulation_3<K> DT;
    DT dt;
    ForwardIterator it;
    for (it = first; it != beyond; ++it) {
        dt.insert(get(point_pmap, *it));
    }

    cov.clear();

    // Compute the VCM
    for (it = first; it != beyond; ++it) {
        typename DT::Vertex_handle vh = dt.nearest_vertex(get(point_pmap, *it));
        Covariance c = Voronoi_covariance_3::voronoi_covariance_3(dt, vh, sphere);
        cov.push_back(c);
    }
}

// Convolve
template < class ForwardIterator,
           class PointPMap,
           class K,
           class Covariance
>
void
vcm_convolve (ForwardIterator first,
              ForwardIterator beyond,
              PointPMap point_pmap,
              const std::vector<Covariance> &cov,
              std::vector<Covariance> &ncov,
              double r,
              const K &)
{
    typedef typename CGAL::Point_3<K> Point;
    typedef typename CGAL::Search_traits_3<K> Traits;
    typedef typename CGAL::Kd_tree<Traits> Tree;
    typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

    ForwardIterator it;

    // Kd tree
    Tree tree;
    for (it = first; it != beyond; ++it) {
        tree.insert(get(point_pmap, *it));
    }

    std::map<Point, size_t> indices;
    size_t i = 0;
    for (it = first; it != beyond; ++it) {
        indices[get(point_pmap, *it)] = i;
        i++;
    }

    // Convolving
    ncov.clear();
    for (it = first; it != beyond; ++it) {
        std::vector<Point> nn;
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere (get(point_pmap, *it), r));

        Covariance m;
        for (size_t k = 0; k < nn.size(); ++k)
            m += cov[indices [nn[k]]];
        ncov.push_back(m);
    }
}

// Construct the covariance matrix
template <class Covariance>
Eigen::Matrix3f
construct_covariance_matrix (Covariance &cov) {
    Eigen::Matrix3f m;

    m(0,0) = cov[0]; m(0,1) = cov[1]; m(0,2) = cov[2];
    m(1,1) = cov[3]; m(1,2) = cov[4];
    m(2,2) = cov[5];

    m(1, 0) = m(0,1); m(2, 0) = m(0, 2); m(2, 1) = m(1, 2);

    return m;
}

// Diagonalize a selfadjoint matrix
bool
diagonalize_selfadjoint_matrix (Eigen::Matrix3f &m,
                                Eigen::Matrix3f &eigenvectors,
                                Eigen::Vector3f &eigenvalues) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(m);

    if (eigensolver.info() != Eigen::Success) {
        return false;
    }

    eigenvalues = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();

    return true;
}

// Determine if a point is on an edge
template <class Covariance>
bool
is_on_edge (Covariance &cov,
            float threshold,
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

// Extract the eigenvector associated to the greatest eigenvalue
template <class Covariance>
bool
extract_greater_eigenvector (Covariance &cov,
                             Eigen::Vector3f &normal) {
    // Construct covariance matrix
    Eigen::Matrix3f m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    Eigen::Vector3f eigenvalues;
    Eigen::Matrix3f eigenvectors;
    if (! diagonalize_selfadjoint_matrix(m, eigenvectors, eigenvalues)) {
        return false;
    }

    // Eigenvalues are already sorted by increasing order
    normal = eigenvectors.col(0);

    return true;
}

/// \endcond

} // namespace internal

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// Estimates normal directions of the `[first, beyond)` range of points
/// using the Voronoi Covariance Measure.
/// The output normals are randomly oriented.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = Vector_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
/// @tparam Covariance Covariance matrix type.

// This variant requires all parameters.
template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel,
           typename Covariance
>
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double R, ///< offset radius.
                      double r, ///< convolution radius.
                      const Kernel & /*kernel*/, ///< geometric traits.
                      const Covariance &) ///< covariance matrix type.
{
    // First, compute the VCM for each point
    std::vector<Covariance> cov;
    size_t N = 20;
    internal::vcm_offset (first, beyond,
                          point_pmap,
                          cov,
                          R,
                          N,
                          Kernel());

    // TODO: debug
    std::ofstream file_offset("offset.p");
    std::copy(cov.begin(), cov.end(),
              std::ostream_iterator<Covariance>(file_offset));

    // Then, convolve it (only when r != 0)
    std::vector<Covariance> ccov;
    if (r == 0) {
        std::copy(cov.begin(), cov.end(), std::back_inserter(ccov));
    } else {
        internal::vcm_convolve(first, beyond,
                               point_pmap,
                               cov,
                               ccov,
                               r,
                               Kernel());
    }

    // TODO: debug
    std::ofstream file_convolve("convolve.p");
    std::copy(ccov.begin(), ccov.end(),
              std::ostream_iterator<Covariance>(file_convolve));

    // Sharp edges
    float threshold = 0.5;
    std::vector<typename Kernel::Point_3> points_on_edges;
    // And finally, compute the normals
    int i = 0;
    for (ForwardIterator it = first; it != beyond; ++it) {
        Eigen::Vector3f enormal;
        Eigen::Vector3f dir;
        if (! internal::extract_greater_eigenvector(ccov[i], enormal)) {
            std::cerr << "Error during extraction of normal: " <<
                "the covariance matrix is not diagonalizable!\n";
            exit(1);
        }

        // TODO
        if (internal::is_on_edge(ccov[i], threshold, dir)) {
            points_on_edges.push_back(get(point_pmap, *it));
        }

        CGAL::Vector_3<Kernel> normal(enormal[0],
                                      enormal[1],
                                      enormal[2]);
        put(normal_pmap, *it, normal);
        i++;
    }

    // TODO: debug
    std::ofstream file_edges("edges.xyz");
    std::copy(points_on_edges.begin(), points_on_edges.end(),
              std::ostream_iterator<typename Kernel::Point_3>(file_edges, "\n"));
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Covariance
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      PointPMap point_pmap,
                      NormalPMap normal_pmap,
                      double R,
                      double r,
                      const Covariance &) {
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;

    vcm_estimate_normals(first, beyond,
                         point_pmap, normal_pmap,
                         R, r,
                         Kernel(),
                         Covariance());
}
/// @endcond

} // namespace CGAL

#endif // CGAL_VCM_ESTIMATE_NORMALS_H
