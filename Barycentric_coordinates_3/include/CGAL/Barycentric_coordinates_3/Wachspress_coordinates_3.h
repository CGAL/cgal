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
        using Dihedral_angle_3 = typename GeomTraits::Compute_approximate_dihedral_angle_3;

        // Number type
        typedef typename GeomTraits::FT FT;

        // Point type.
        typedef typename GeomTraits::Point_3 Point_3;

        // Mesh type
        typedef CGAL::Surface_mesh<Point_3> Mesh;

        // Custom types
        typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
        typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
        typedef typename boost::property_map<Mesh,CGAL::vertex_point_t>::type Point_property_map; 
        typedef typename CGAL::Face_around_target_circulator<Mesh> Face_circulator;
        
        // Constructor just initializing data
        Wachspress_coordinates_3(
            const Mesh& mesh,
            const GeomTraits traits = GeomTraits()) :
        m_mesh(mesh), 
        m_traits(traits),
        m_dihedral_angle_3(m_traits.compute_approximate_dihedral_angle_3_object()){}

        //Test function to return dihedral angles of edges opposite to first vertex
        std::vector<FT> dihedral_first(){

            ppm = get(CGAL::vertex_point, m_mesh);
            vertex_descriptor vd = (vertices(m_mesh).begin())[0];

            Point_3 v0 = get(ppm, vd);

            m_face_circulator = Face_circulator(m_mesh.halfedge(vd), m_mesh);
            Face_circulator done(m_face_circulator);

            do{

                std::cout << *m_face_circulator++ << "\n"; 
            }while(m_face_circulator!=done);

            return std::vector<FT>();
        }

    
    
    private:

        // Fields
        const GeomTraits m_traits;
        const Mesh& m_mesh;
        const Dihedral_angle_3 m_dihedral_angle_3;

        // Custom variables
        Point_property_map ppm;
        Face_circulator m_face_circulator;


};


}
}

#endif