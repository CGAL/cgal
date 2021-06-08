#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>

//Typedefs
using Kernel = CGAL::Simple_cartesian<double>;

using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Mesh =  CGAL::Surface_mesh<Point_3>;

int main(){

    // Tetrahedron Mesh
    Mesh ms;

    // Construct tetrahedron
    const Point_3 p0(0.0, 0.0, 0.0);
    const Point_3 p1(1.0, 0.0, 0.0);
    const Point_3 p2(0.0, 1.0, 0.0);
    const Point_3 p3(0.0, 0.0, 1.0);

    CGAL::make_tetrahedron(p0, p1, p2, p3, ms);

    CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Kernel> ws(ms, Kernel());
    ws.dihedral_first();

    return EXIT_SUCCESS;
}