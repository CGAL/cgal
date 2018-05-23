#ifndef CGAL_CONCAVITY_H
#define CGAL_CONCAVITY_H

namespace CGAL
{
namespace internal
{

    template <class GeomTraits>
    class Concavity
    {
    public:
        Concavity(const GeomTraits& traits)
        : m_traits(traits)
        {}

        template <class ClusterMesh>
        double calc(const ClusterMesh& mesh)
        {
            //TODO: implement
            return 5;
        }

    private:
        const GeomTraits& m_traits;
    };

}
}

#endif // CGAL_CONCAVITY_H
