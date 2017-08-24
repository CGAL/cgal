#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>


#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif

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

    Mesh mesh;
    Mesh originalMesh;
    std::ifstream input;
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::ofstream output;
#endif

    // flat test
    input.open("data/polygon.off");
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }
    input.close();

    originalMesh = mesh;
    CGAL::Polygon_mesh_processing::curvature_flow_smoothing(mesh);

    double dist = CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance
        <TAG>(originalMesh, mesh, CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(1000));

    if(dist > 1e-3)
    {
        return EXIT_FAILURE;
    }

    // half-sphere test
    input.open("data/sphere.off");
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }
    input.close();


    // select half sphere
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout<<"all faces: "<<faces(mesh).size()<<std::endl;
#endif

    std::set<face_index> selected_faces = select(mesh, 0);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout<<"selected faces: "<<selected_faces.size()<<std::endl;
#endif
    CGAL::Polygon_mesh_processing::isotropic_remeshing(selected_faces, 0.1, mesh);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    output.open("data/half-sphere.off");
    output << mesh;
    output.close();
#endif

    originalMesh = mesh;
    CGAL::Polygon_mesh_processing::curvature_flow_smoothing(mesh);

    dist = CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance
        <TAG>(originalMesh, mesh, CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(1000));

    if(dist > 1e-1)
    {
        return EXIT_FAILURE;
    }

#ifdef CGAL_PMP_REMESHING_VERBOSE
    output.open("data/half-sphere_smoothed.off");
    output << mesh;
    output.close();
#endif

    return 0;
}
