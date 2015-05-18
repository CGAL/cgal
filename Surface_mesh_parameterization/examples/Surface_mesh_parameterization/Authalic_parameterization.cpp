#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <cstdlib>
#include <fstream>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "  Discrete Authalic Parameterization" << std::endl;
    std::cerr << "  circle border" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 != 1)
    {
        std::cerr << "Usage: " << argv[0] << " input_file.off" << std::endl;
        return(EXIT_FAILURE);
    }

    // File name is:
    const char* input_filename  = argv[1];

    //***************************************
    // Read the mesh
    //***************************************

    // Read the mesh
    std::ifstream stream(input_filename);
    Polyhedron mesh;
    stream >> mesh;
    if(!stream || !mesh.is_valid() || mesh.empty())
    {
        std::cerr << "Error: cannot read OFF file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }

    //***************************************
    // Create Polyhedron adaptor
    // Note: no cutting => we support only
    // meshes that are topological disks
    //***************************************

    typedef std::map<vertex_descriptor,std::size_t> V_index_map;
    typedef std::map<halfedge_descriptor,std::size_t> H_index_map;
    typedef std::map<halfedge_descriptor,Point_2> H_uv_map;
    typedef boost::associative_property_map<V_index_map> V_index_pmap;
    typedef boost::associative_property_map<H_index_map> H_index_pmap;
    typedef boost::associative_property_map<H_uv_map> H_uv_pmap;
    V_index_map vim; V_index_pmap vipm(vim);
    H_index_map him; H_index_pmap hipm(him);
    H_uv_map huvm; H_uv_pmap huvpm(huvm);

    std::size_t i = 0;
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
      vim[vd] = i++;
    } 

    i = 0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
      him[hd] = i++;
    }

    typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron, V_index_pmap, H_index_pmap, H_uv_pmap>
                                            Parameterization_polyhedron_adaptor;
    Parameterization_polyhedron_adaptor mesh_adaptor(mesh, vipm,hipm,huvpm);

    //***************************************
    // Discrete Authalic Parameterization
    // (defaults are circular border and Eigen solver)
    //***************************************

    typedef CGAL::Discrete_authalic_parameterizer_3<Parameterization_polyhedron_adaptor>
                                                        Parameterizer;

    Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor, Parameterizer());
    switch(err) {
    case Parameterizer::OK: // Success
        break;
    case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
    case Parameterizer::ERROR_NON_TRIANGULAR_MESH:   
    case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:     
    case Parameterizer::ERROR_BORDER_TOO_SHORT:    
        std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
        break;
    default: // Error
        std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
        break;
    };

    //***************************************
    // Output
    //***************************************

    // Raw output: dump (u,v) pairs

    std::cout <<"OFF\n" << num_vertices(mesh) << " " << num_faces(mesh) << " 0\n";
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
      // (u,v) pair is stored in any halfedge
      std::cout << huvm[vd->halfedge()] << " 0" << std::endl;
    }
    BOOST_FOREACH(face_descriptor fd, faces(mesh)){
      std::cout << "3";
      BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(halfedge(fd,mesh),mesh)){
        std::cout << " " <<   mesh_adaptor.info(vd)->index();
      }
      std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
