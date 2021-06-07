#ifndef CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H
#define CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>

namespace CGAL {
namespace Barycentric_coordinates {

template<typename GeomTraits>
class Wachspress_coordinates_3{

    public:

        // Dihedral angle calculation
        using Dihedral_angle_3 = GeomTraits::Compute_approximate_dihedral_angle_3;

        // Number type
        typedef typename GeomTraits::FT FT;

        // Point type.
        typedef typename GeomTraits::Point_3 Point_3;

        // Mesh type
        typedef CGAL::Surface_mesh<Point_3> Mesh;
        
        // Constructor just initializing data
        Wachspress_coordinates_3(
            const Mesh& mesh,
            const GeomTraits traits = GeomTraits()) :
        m_mesh(mesh), 
        m_traits(traits){}

        // Compute weights for query point q
        template<OutputIterator>
        OutputIterator compute(
            const Point_3& query, OutIterator output, const bool normalize){
            
            return output;
        }
    
    
    private:

        // Fields
        const GeomTraits m_traits;
        const Mesh m_mesh;

        const Dihedral_angle_3 m_dihedral_angle_3;


};



}
}

#endif