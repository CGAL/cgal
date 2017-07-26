#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef typename Mesh::vertex_index vertex_index;
typedef typename Mesh::face_index face_index;
typedef typename Mesh::halfedge_index halfedge_index;

std::set<face_index> select(Mesh mesh, double x0)
{
    std::set<face_index> f_selection;

    typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VertexPointMap;
    VertexPointMap vpmap = get(CGAL::vertex_point, mesh);


    for(face_index f : faces(mesh))
    {
        halfedge_index he = halfedge(f, mesh);
        for(vertex_index v : vertices_around_face(he, mesh))
        {
            if(get(vpmap, v).x() < x0)
                f_selection.insert(f);
            continue;
        }
    }

    return f_selection;
}


int main(int argc, char* argv[]){


    const char* filename = "data/sphere.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    std::cout<<"all faces: "<<faces(mesh).size()<<std::endl;
    std::set<face_index> selected_faces = select(mesh, 0);
    std::cout<<"selected faces: "<<selected_faces.size()<<std::endl;

    CGAL::Polygon_mesh_processing::isotropic_remeshing(selected_faces, 0.05, mesh);

    std::ofstream output;
    /*
    output.open("data/sphere_half-subsampled.off");
    output << mesh;
    output.close();
    */

    CGAL::Polygon_mesh_processing::curvature_flow(mesh);


    output.open("data/sphere-half-subsampled_smoothed.off");
    output << mesh;
    output.close();



    return 0;
}
