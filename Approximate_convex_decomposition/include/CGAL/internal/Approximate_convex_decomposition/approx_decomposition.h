#ifndef CGAL_APPROX_DECOMPOSITION_H
#define CGAL_APPROX_DECOMPOSITION_H

namespace CGAL
{
namespace internal
{

template <class TriangleMesh,
          class GeomTraits
          >
class Approx_decomposition
{
public:
    Approx_decomposition(const TriangleMesh& mesh, const GeomTraits& traits)
    : m_mesh(mesh)
    , m_traits(traits)
    {}

    template <class FacetPropertyMap, class PointPropertyMap>
    std::size_t decompose(FacetPropertyMap facet_ids, PointPropertyMap point_ids, double concavity_threshold, std::size_t min_number_of_clusters)
    {
        //TODO: implement
        typedef typename TriangleMesh::Facet_const_iterator Iterator;
        for (Iterator it = m_mesh.facets_begin(); it != m_mesh.facets_end(); ++it)
        {   
            facet_ids[it] = -1;
        }

        return 0;
    }

private:
    const TriangleMesh& m_mesh;
    const GeomTraits& m_traits;
};        

}
}

#endif // CGAL_APPROX_DECOMPOSITION_H
