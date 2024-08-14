#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SurfaceMesh;
typedef CGAL::Nef_polyhedron_3<Kernel> NefPolyhedron;


int main()
{
    SurfaceMesh surfaceMesh;
    SurfaceMesh::Vertex_index v0 = surfaceMesh.add_vertex(Point_3(-1, 0, 0));
    SurfaceMesh::Vertex_index v1 = surfaceMesh.add_vertex(Point_3(1, 0, 0));
    SurfaceMesh::Vertex_index v2 = surfaceMesh.add_vertex(Point_3(0, 1, 0));

    SurfaceMesh::Vertex_index v3 = surfaceMesh.add_vertex(Point_3(-1, 0, 1));
    SurfaceMesh::Vertex_index v4 = surfaceMesh.add_vertex(Point_3(1, 0, 1));
    SurfaceMesh::Vertex_index v5 = surfaceMesh.add_vertex(Point_3(0, 1, 1));

    surfaceMesh.add_face(v0, v1, v2);
    surfaceMesh.add_face(v3, v4, v5);

    make_tetrahedron(Point_3(-1, 0, 10),
                     Point_3(1, 0, 10),
                     Point_3(0, 1, 10),
                     Point_3(-1, 0, 11),
                     surfaceMesh);

    std::cout << "Before conversion, number_of_faces: " << surfaceMesh.number_of_faces() << std::endl;

    NefPolyhedron nefPoly(surfaceMesh);
    std::cout << "NefPolyhedron, number_of_faces: " << nefPoly.number_of_facets() << std::endl;
    SurfaceMesh convertedSurfaceMesh;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(nefPoly, convertedSurfaceMesh, true);
    std::cout << "After conversion, number_of_faces: " << convertedSurfaceMesh.number_of_faces() << std::endl;
    std::ofstream("out.off") << convertedSurfaceMesh;
    assert(vertices(convertedSurfaceMesh).size()==10);
    assert(faces(convertedSurfaceMesh).size()==6);
    return EXIT_SUCCESS;
}
