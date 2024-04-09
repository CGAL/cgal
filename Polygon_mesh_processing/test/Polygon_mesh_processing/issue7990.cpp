#define CGAL_PMP_SMOOTHING_DEBUG 1
#define CGAL_PMP_SMOOTHING_OUTPUT_INTERMEDIARY_STEPS 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
//#include <CGAL/draw_surface_mesh.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <iostream>
#include <set>
#include <string>
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;
int main(int argc, char* argv[])
{
        //QCoreApplication::addLibraryPath("C:/Qt/5.15.2/msvc2019_64/plugins");
        const std::string filename = (argc > 1) ? argv[1] : "./data/issue7990.off";
        Mesh mesh;
        if(!PMP::IO::read_polygon_mesh(filename, mesh))
        {
                std::cerr << "Invalid input." << std::endl;
                return 1;
        }
        //draw(mesh);
        const unsigned int nb_iterations = (argc > 2) ? std::atoi(argv[2]) : 3;
        const double time = (argc > 3) ? std::atof(argv[3]) : 0.001;
        std::set<Mesh::Vertex_index> constrained_vertices;
        for(Mesh::Vertex_index v : vertices(mesh))
        {
              //  if(is_border(v, mesh))
              //          constrained_vertices.insert(v);
        }
        std::cout << "Constraining: " << constrained_vertices.size() << " border vertices" << std::endl;

        constrained_vertices.insert(Mesh::Vertex_index(7));
        constrained_vertices.insert(Mesh::Vertex_index(8));
        constrained_vertices.insert(Mesh::Vertex_index(9));
        constrained_vertices.insert(Mesh::Vertex_index(10));
        constrained_vertices.insert(Mesh::Vertex_index(25));
        constrained_vertices.insert(Mesh::Vertex_index(26));
        constrained_vertices.insert(Mesh::Vertex_index(27));
        constrained_vertices.insert(Mesh::Vertex_index(28));
        constrained_vertices.insert(Mesh::Vertex_index(13));
        constrained_vertices.insert(Mesh::Vertex_index(16));
        constrained_vertices.insert(Mesh::Vertex_index(19));
        constrained_vertices.insert(Mesh::Vertex_index(22));
        CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);


        std::cerr.rdbuf(std::cout.rdbuf());
        std::set<Mesh::Face_index> faces_to_smooth;
#if 0
        faces_to_smooth.insert(faces(mesh).begin(), faces(mesh).end()); // case 1: all faces -> success
#else
        faces_to_smooth.insert(Mesh::Face_index(12));
        faces_to_smooth.insert(Mesh::Face_index(13));
        faces_to_smooth.insert(Mesh::Face_index(14));
        faces_to_smooth.insert(Mesh::Face_index(15));
        faces_to_smooth.insert(Mesh::Face_index(16));
        faces_to_smooth.insert(Mesh::Face_index(17));

        faces_to_smooth.insert(Mesh::Face_index(22));
        faces_to_smooth.insert(Mesh::Face_index(23));
        faces_to_smooth.insert(Mesh::Face_index(24));
        faces_to_smooth.insert(Mesh::Face_index(25));
        faces_to_smooth.insert(Mesh::Face_index(26));
        faces_to_smooth.insert(Mesh::Face_index(27));

        faces_to_smooth.insert(Mesh::Face_index(32));
        faces_to_smooth.insert(Mesh::Face_index(33));
        faces_to_smooth.insert(Mesh::Face_index(34));
        faces_to_smooth.insert(Mesh::Face_index(35));
        faces_to_smooth.insert(Mesh::Face_index(36));
        faces_to_smooth.insert(Mesh::Face_index(37));
#endif
        std::cout << "Smoothing shape... (" << nb_iterations << " iterations)" << std::endl;
        PMP::smooth_shape(faces_to_smooth, mesh, time, CGAL::parameters::number_of_iterations(nb_iterations)
                                                                                                   .vertex_is_constrained_map(vcmap));
        CGAL::IO::write_polygon_mesh("./issue7990_smoothed.off", mesh, CGAL::parameters::stream_precision(17));
        std::cout << "Done!" << std::endl;
        //draw(mesh);


        return EXIT_SUCCESS;
}
