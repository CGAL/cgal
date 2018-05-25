#ifndef CGAL_CONCAVITY_H
#define CGAL_CONCAVITY_H

namespace CGAL
{
namespace internal
{

    template <class TriangleMesh, class GeomTraits>
    class Concavity
    {
    public:
        Concavity(const TriangleMesh& mesh, const GeomTraits& traits)
        : m_mesh(mesh)
        , m_traits(traits)
        {}

        template <class FacetPropertyMap>
        double compute(FacetPropertyMap facet_ids, std::size_t cluster_id)
        {
            //TODO: implement
            return 5;
        }

    private:
        const TriangleMesh& m_mesh;
        const GeomTraits& m_traits;
    };

}
}

#endif // CGAL_CONCAVITY_H
