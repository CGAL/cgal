#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>

//Typedefs
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

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

    // Instantiate some interior query points for which we compute coordinates.
    const std::vector<Point_3> queries = {  Point_3(0.25f , 0.25f, 0.25f), Point_3(0.3f, 0.2f, 0.3f),
                                            Point_3(0.1f, 0.1f, 0.1f)};

    CGAL::make_tetrahedron(p0, p1, p2, p3, ms);

    CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> ws(ms);

    for(auto q:queries){

        // Store results
        std::vector<FT> coordinates;
        ws(q, std::back_inserter(coordinates));

        std::cout << "Coordinates: " << "(" << q << ") ---->>> " << coordinates[0] << " " << 
        coordinates[1] << " " << coordinates[2] << " " << coordinates[3] << std::endl;
    }

    return EXIT_SUCCESS;
}