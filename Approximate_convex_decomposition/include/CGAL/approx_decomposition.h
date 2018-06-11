#ifndef CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
#define CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/internal/Approximate_convex_decomposition/approx_decomposition.h>
#include <CGAL/internal/Approximate_convex_decomposition/concavity.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL
{


/*
 * @brief Function computing concavity value of a cluster.
 * 
 * Cluster is a subset of connected faces in a triangle mesh.
 * Concavity value is the largest distance from a vertex in a cluster to the projected point onto the convex hull of a mesh.
 *
 * @return concativity value.
 */
template <class TriangleMesh, class FacePropertyMap,
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
 * @brief Function computing concavity value of a triangle mesh.
 * 
 * Concavity value is the largest distance from a vertex in a mesh to the projected point onto the convex hull of a mesh.
 *
 * @return concativity value.
 */
template <class TriangleMesh,
          class GeomTraits = CGAL::Exact_predicates_inexact_constructions_kernel
          >
double
concavity_value(const TriangleMesh& mesh,
                const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));

    internal::Concavity<TriangleMesh, GeomTraits> algorithm(mesh, traits);
    return algorithm.compute();
}


/*
 * @brief Function computing the approximate convex decomposition of a triangle mesh.
 *
 * This function fills a property map which associates a cluster-id (in the range [0, 'number_of_clusters'-1]) to each face.
 * Clusters are subsets of connected faces in a triangle mesh which concavity values are less than 'concavity_threshold'.
 *
 * @return number of clusters.
 */
template <class TriangleMesh, class FacePropertyMap,
          class GeomTraits = CGAL::Exact_predicates_inexact_constructions_kernel
          >
std::size_t
convex_decomposition(const TriangleMesh& mesh,
                     FacePropertyMap face_ids,
                     double concavity_threshold = 0.05,
                     std::size_t min_number_of_clusters = 1,
                     const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));
    CGAL_precondition(CGAL::num_faces(mesh) >= min_number_of_clusters);

    internal::Approx_decomposition<TriangleMesh, GeomTraits> algorithm(mesh, traits);
    return algorithm.decompose(face_ids, concavity_threshold, min_number_of_clusters);
}


} //namespace CGAL

#endif //CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
