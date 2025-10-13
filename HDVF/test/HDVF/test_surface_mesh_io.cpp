#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/HDVF/Surface_mesh_io.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;
typedef HDVF::Hdvf_traits_3<K> Traits;
typedef CGAL::Homological_discrete_vector_field::Surface_mesh_io<SurfaceMesh,Traits> SurfaceMeshIO;

int main() {
    SurfaceMesh mesh;
    std::string filename("data/three_triangles.off");

    CGAL::IO::read_polygon_mesh(filename, mesh);

    SurfaceMeshIO mesh_io(mesh);

    return 0;
}


