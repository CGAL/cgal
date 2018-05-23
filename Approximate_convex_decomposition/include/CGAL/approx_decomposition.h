#ifndef CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
#define CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H

#define CGAL_DEFAULT_CONCAVITY_THRESHOLD 0.5

#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/internal/Approximate_convex_decomposition/approx_decomposition.h>
#include <CGAL/internal/Approximate_convex_decomposition/concavity.h>
#include <vector>
#include <map>

namespace CGAL
{

/*
 * @brief Function computing concativy of a cluster.
 * Cluster is a closed triangular mesh.
 *
 * @return concativity value.
 */
template <class ClusterMesh,
          class GeomTraits = typename Kernel_traits<typename ClusterMesh::Point_3>::Kernel
          >
double
cluster_concavity_value(const ClusterMesh& mesh,
                        const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));
    CGAL_precondition(CGAL::is_closed(mesh));

    internal::Concavity<GeomTraits> algorithm(traits);
    return algorithm.calc(mesh);
}

/*
 * @brief Function computing the decomposition of a triangular mesh.
 *
 * This function fills a property map which associates a cluster-id (in [0, 'number_of_clusters'-1] or -1) to each facet.
 * If a facet doesn't belong to any cluster, it's cluster-id is -1.
 * Cluster is a closed triangular mesh that satisfies 'concavity_threshold'.
 *
 * @return number of clusters.
 */
template <class TriangleMesh, class FacetPropertyMap,
          class PointPropertyMap = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
          class GeomTraits = typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
          >
std::size_t
convex_decomposition(const TriangleMesh& mesh,
                     FacetPropertyMap facet_ids,
                     double concavity_threshold = CGAL_DEFAULT_CONCAVITY_THRESHOLD,
                     PointPropertyMap ppmap = PointPropertyMap(),
                     const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));

    internal::Approx_decomposition<TriangleMesh, GeomTraits> algorithm(mesh, traits);
    return algorithm.decompose(facet_ids, ppmap, concavity_threshold);
}

/*
 * @brief Function constructing the clusters given a mesh and a cluster-id per facet (decomposition).
 *
 * @return none.
 */
 template <class TriangleMesh, class FacetPropertyMap, class ClusterMesh,
           class GeomTraits = typename Kernel_traits<typename TriangleMesh::Point_3>::Kernel
           >
 void
 convex_decomposition_clusters(const TriangleMesh& mesh,
                               std::vector<ClusterMesh>& clusters,
                               FacetPropertyMap facet_ids,
                               const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));

    internal::Approx_decomposition<TriangleMesh, GeomTraits> algorithm(mesh, traits);
    return algorithm.construct_clusters(clusters, facet_ids);
}

/*
 * @brief Function constructing the clusters given a triangular mesh.
 *
 * This function implicitly computes the decomposition of a triangular mesh and then fills array with clusters.
 * Cluster is a closed triangular mesh that satisfies 'concavity_threshold'.
 *
 * @return none
*/
template <class TriangleMesh, class ClusterMesh,
          class GeomTraits = typename Kernel_traits<typename TriangleMesh::Point_3>::Kernel
          >
void
convex_decomposition_clusters(const TriangleMesh& mesh,
                              std::vector<ClusterMesh>& clusters,
                              double concavity_threshold = CGAL_DEFAULT_CONCAVITY_THRESHOLD,
                              const GeomTraits& traits = GeomTraits())
{
    CGAL_precondition(CGAL::is_triangle_mesh(mesh));

    typedef typename TriangleMesh::Facet_const_handle Facet_const_handle;
    typedef std::map<Facet_const_handle, int> Facet_int_map;
    Facet_int_map facet_map;
    boost::associative_property_map<Facet_int_map> facet_property_map(facet_map);
    
    convex_decomposition(mesh, facet_property_map, concavity_threshold);
    convex_decomposition_clusters(mesh, clusters, facet_property_map);
}


} //namespace CGAL

#endif //CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
