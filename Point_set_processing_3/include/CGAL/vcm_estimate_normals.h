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
// Offset
template < typename ForwardIterator,
           typename PointPMap,
           class K,
           class Covariance
>
void
vcm_offset (ForwardIterator first,
            ForwardIterator beyond,
            PointPMap point_pmap,
            std::vector<Covariance> &cov,
            double R, size_t N,
            const K &)
{
    typedef CGAL::Delaunay_triangulation_3<K> DT;

    typename CGAL::Voronoi_covariance_3::Sphere_discretization<K> sphere(R, N);
    DT dt;
    ForwardIterator it;
    for (it = first; it != beyond; it++) {
        dt.insert(get(point_pmap, *it));
    }

    cov.clear();

    for (it = first; it != beyond; it++) {
        typename DT::Vertex_handle vh = dt.nearest_vertex(get(point_pmap, *it));
        Covariance c = CGAL::Voronoi_covariance_3::voronoi_covariance_3(dt, vh, sphere);
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

    Tree tree;
    for (it = first; it != beyond; it++) {
        tree.insert(get(point_pmap, *it));
    }

    std::map<Point, size_t> indices;
    size_t i = 0;
    for (it = first; it != beyond; it++) {
        indices[get(point_pmap, *it)] = i;
        i++;
    }

    ncov.clear();
    for (it = first; it != beyond; it++) {
        std::vector<Point> nn;
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere (get(point_pmap, *it), r));

        Covariance m;
        for (size_t k = 0; k < nn.size(); ++k)
            m += cov[indices [nn[k]]];
        ncov.push_back(m);
    }
}

namespace internal
{
    int
    max_vector3f (Eigen::Vector3f & v)
    {
        float max_eigenvalue = std::max(std::max(v[0], v[1]), v[2]);

        if (max_eigenvalue == v[0])
            return 0;

        if (max_eigenvalue == v[1])
            return 1;

        return 2;
    }

    template <class Covariance>
    bool
    extract_greater_eigenvector (Covariance &cov,
                                 Eigen::Vector3f &normal) {
        Eigen::Matrix3f m;

        // Construct covariance matrix
        m(0,0) = cov[0]; m(0,1) = cov[1]; m(0,2) = cov[2];
        m(1,1) = cov[3]; m(1,2) = cov[4];
        m(2,2) = cov[5];

        m(1, 0) = m(0,1); m(2, 0) = m(0,2); m(2,1) = m(1,2);

        // Diagonalizing the matrix
        Eigen::Vector3f eigenvalues;
        Eigen::Matrix3f eigenvectors;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(m);

        if (eigensolver.info() != Eigen::Success) {
            return false;
        }

        eigenvalues = eigensolver.eigenvalues();
        eigenvectors = eigensolver.eigenvectors();

        int index_max_eigenvalue = max_vector3f(eigenvalues);
        normal = eigenvectors.col(index_max_eigenvalue);

        return true;
    }
} // namespace internal

template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel,
           typename Covariance
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      PointPMap point_pmap,
                      NormalPMap normal_pmap,
                      double R,
                      double r,
                      const Kernel & /*kernel*/,
                      const Covariance &) {
    // First, compute the VCM for each point
    std::vector<Covariance> cov;
    size_t N = 20;
    vcm_offset (first,
                beyond,
                point_pmap,
                cov,
                R,
                N,
                Kernel());

    // Then, convolve it
    std::vector<Covariance> ccov;
    vcm_convolve(first,
                 beyond,
                 point_pmap,
                 cov,
                 ccov,
                 r,
                 Kernel());

    // And finally, compute the normals
    int  i = 0;
    for (ForwardIterator it = first; it != beyond; it++) {
        Eigen::Vector3f enormal;
        if (! internal::extract_greater_eigenvector(cov[i], enormal)) {
            std::cerr << "Error during extraction of normal" << std::endl;
            exit(1);
        }
        CGAL::Vector_3<Kernel> normal(enormal[0],
                                      enormal[1],
                                      enormal[2]);
        put(normal_pmap, *it, normal);
        i++;
    }
}

} // namespace CGAL

#endif // CGAL_VCM_ESTIMATE_NORMALS_H
