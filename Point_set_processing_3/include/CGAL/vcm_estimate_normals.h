#ifndef CGAL_VCM_ESTIMATE_NORMALS_H
#define CGAL_VCM_ESTIMATE_NORMALS_H

#include <CGAL/property_map.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/voronoi_covariance_3.h>
#include <CGAL/vcm_utilities.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <iterator>
#include <vector>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// @cond SKIP_IN_MANUAL
/// Computes the VCM for each point in the property map.
/// The matrix is computed by intersecting the Voronoi cell
/// of a point and a sphere whose radius is `R` and discretized
/// by `N` planes.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
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
/// @endcond

/// @cond SKIP_IN_MANUAL
// Convolve using a radius.
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
    typedef typename K::Point_3 Point;
    typedef CGAL::Search_traits_3<K> Traits;
    typedef CGAL::Kd_tree<Traits> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

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
/// @endcond

/// @cond SKIP_IN_MANUAL
// Convolve using neighbors.
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
              unsigned int nb_neighbors_convolve,
              const K &)
{
    typedef typename K::Point_3 Point;
    typedef CGAL::Search_traits_3<K> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;

    ForwardIterator it;

    // Search tree
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
        Neighbor_search search(tree, get(point_pmap, *it), nb_neighbors_convolve);
        std::vector<Point> nn;

        for (typename Neighbor_search::iterator nit = search.begin();
             nit != search.end();
             ++nit) {
            nn.push_back(nit->first);
        }

        Covariance m;
        for (size_t k = 0; k < nn.size(); ++k)
            m += cov[indices [nn[k]]];
        ncov.push_back(m);
    }
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// Compute the VCM and make the convolution using a radius.
template < class ForwardIterator,
           class PointPMap,
           class Kernel,
           class Covariance
>
void
vcm_offset_and_convolve (ForwardIterator first,
                         ForwardIterator beyond,
                         PointPMap point_pmap,
                         std::vector<Covariance> &ccov,
                         double R,
                         double r,
                         const Kernel & k)
{
    // First, compute the VCM for each point
    std::vector<Covariance> cov;
    size_t N = 20;
    internal::vcm_offset (first, beyond,
                          point_pmap,
                          cov,
                          R,
                          N,
                          k);

    // Then, convolve it (only when r != 0)
    if (r == 0) {
        std::copy(cov.begin(), cov.end(), std::back_inserter(ccov));
    } else {
        internal::vcm_convolve(first, beyond,
                               point_pmap,
                               cov,
                               ccov,
                               r,
                               k);
    }
}
/// @endcond

} // namespace internal

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// @cond SKIP_IN_MANUAL
/// Estimates normal directions of the `[first, beyond)` range of points
/// using the Voronoi Covariance Measure.
/// The output normals are randomly oriented.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
///        It can be omitted if ForwardIterator value_type is convertible to `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = `Kernel::Vector_3`.
/// @tparam Kernel Geometric traits class.
/// @tparam Covariance Covariance matrix type..
/// @pre If `nb_neighbors_convolve` is equal to -1, then the convolution is made using a radius.
/// On the contrary, if `nb_neighbors_convolve` is different from -1, the convolution is made using
/// this number of neighbors.

// This variant requires all of the parameters.
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
                      const Kernel & k, ///< geometric traits.
                      const Covariance &, ///< covariance matrix type.
                      int nb_neighbors_convolve = -1 ///< number of neighbors used during the convolution.
)
{
    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    if (nb_neighbors_convolve == -1) {
        internal::vcm_offset_and_convolve(first, beyond,
                                          point_pmap,
                                          cov,
                                          R,
                                          r,
                                          k);
    } else {
        internal::vcm_offset(first, beyond,
                             point_pmap,
                             cov,
                             R,
                             20,
                             k);

        std::vector<Covariance> ccov;
        internal::vcm_convolve(first, beyond,
                               point_pmap,
                               cov,
                               ccov,
                               (unsigned int) nb_neighbors_convolve,
                               k);

        cov.clear();
        std::copy(ccov.begin(), ccov.end(), std::back_inserter(cov));
    }

    // And finally, compute the normals
    int i = 0;
    for (ForwardIterator it = first; it != beyond; ++it) {
        Eigen::Vector3f enormal;
        Eigen::Vector3f dir;
        internal::extract_greater_eigenvector(cov[i], enormal);

        typename Kernel::Vector_3 normal(enormal[0],
                                         enormal[1],
                                         enormal[2]);
        put(normal_pmap, *it, normal);
        i++;
    }
}
/// @endcond

/// \ingroup PkgPointSetProcessing
/// Estimates normal directions of the `[first, beyond)` range of points
/// using the Voronoi Covariance Measure with a radius for the convolution.
/// The output normals are randomly oriented.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
///        It can be omitted if ForwardIterator value_type is convertible to `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = `Kernel::Vector_3`.

// This variant deduces the kernel from the point property map
// and uses a radius for the convolution.
template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double R, ///< offset radius: radius of the sphere to intersect the Voronoi cell with.
                      double r ///< convolution radius: all points in a sphere with this radius will be convolved.
)
{
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::FT FT;
    typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

    vcm_estimate_normals(first, beyond,
                         point_pmap, normal_pmap,
                         R, r,
                         Kernel(),
                         Covariance());
}

/// \ingroup PkgPointSetProcessing
/// Estimates normal directions of the `[first, beyond)` range of points
/// using the Voronoi Covariance Measure with a number of neighbors for the convolution.
/// The output normals are randomly oriented.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
///        It can be omitted if ForwardIterator value_type is convertible to `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = `Kernel::Vector_3`.

// This variant deduces the kernel from the point property map
// and uses a number of neighbors for the convolution.
template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double R, ///< offset radius: radius of the sphere to intersect the Voronoi cell with.
                      unsigned int nb_neighbors_convolve ///< number of neighbor points used for the convolution.
)
{
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::FT FT;
    typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

    vcm_estimate_normals(first, beyond,
                         point_pmap, normal_pmap,
                         R, 0,
                         Kernel(),
                         Covariance(),
                         nb_neighbors_convolve);
}

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map
// and use a radius for the convolution.
template < typename ForwardIterator,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      NormalPMap normal_pmap,
                      double R,
                      double r) {
    vcm_estimate_normals(first, beyond,
                         make_identity_property_map(typename std::iterator_traits<ForwardIterator>::value_type()),
                         normal_pmap,
                         R, r);
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map
// and use a number of neighbors for the convolution.
template < typename ForwardIterator,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      NormalPMap normal_pmap,
                      double R,
                      unsigned int nb_neighbors_convolve) {
    vcm_estimate_normals(first, beyond,
                         make_identity_property_map(typename std::iterator_traits<ForwardIterator>::value_type()),
                         normal_pmap,
                         R, nb_neighbors_convolve);
}
/// @endcond

} // namespace CGAL

#endif // CGAL_VCM_ESTIMATE_NORMALS_H
