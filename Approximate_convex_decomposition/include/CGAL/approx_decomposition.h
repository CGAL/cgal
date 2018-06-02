#ifndef CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
#define CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/internal/Approximate_convex_decomposition/approx_decomposition.h>
#include <CGAL/internal/Approximate_convex_decomposition/concavity.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <map>

namespace CGAL
{


/*
 * @brief Function computing concavity value of a cluster.
 * Cluster is a subset of connected faces in a triangle mesh.
 *
 * @return concativity value (the largest squared distance from a point in a cluster to the projected point onto the convex hull of a mesh).
 */
template <class TriangleMesh, class FacePropertyMap,
//          class GeomTraits = CGAL::Simple_cartesian<double>
          class GeomTraits = CGAL::Exact_predicates_inexact_constructions_kernel
          >
double
concavity_value(const TriangleMesh& mesh,
                FacePropertyMap face_ids,
                std::size_t cluster_id,
                const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));

    internal::Concavity<TriangleMesh, GeomTraits> algorithm(mesh, traits);
    return algorithm.compute(face_ids, cluster_id);
}


/*
 * @brief Function computing the approximate convex decomposition of a triangular mesh.
 *
 * This function fills a property map which associates a cluster-id (in [0, 'number_of_clusters'-1] or -1) to each face.
 * If a face doesn't belong to any cluster, it's cluster-id is set to -1.
 * Clusters are subsets of connected faces in a triangular mesh which concavity values satisfy 'concavity_threshold'.
 *
 * @return number of clusters.
 */
template <class TriangleMesh, class FacePropertyMap,
          class PointPropertyMap = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
          class GeomTraits = typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
          >
std::size_t
convex_decomposition(const TriangleMesh& mesh,
                     FacePropertyMap face_ids,
                     double concavity_threshold = 100,
                     std::size_t min_number_of_clusters = 1,
                     PointPropertyMap ppmap = PointPropertyMap(),
                     const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));

    internal::Approx_decomposition<TriangleMesh, GeomTraits> algorithm(mesh, traits);
    return algorithm.decompose(face_ids, ppmap, concavity_threshold, min_number_of_clusters);
}


} //namespace CGAL

#endif //CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
