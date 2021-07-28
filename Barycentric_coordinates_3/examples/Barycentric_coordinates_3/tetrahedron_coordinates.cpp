#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

//Typedefs
using Kernel = CGAL::Simple_cartesian<double>;

using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;

int main(){

    // Construct tetrahedron
    const Point_3 p0(0.0, 0.0, 0.0);
    const Point_3 p1(1.0, 0.0, 0.0);
    const Point_3 p2(0.0, 1.0, 0.0);
    const Point_3 p3(0.0, 0.0, 1.0);

    // Instantiate some interior, boundary, and exterior query points for which we compute coordinates.
    const std::vector<Point_3> queries = {
        Point_3(0.25f , 0.25f, 0.25f), Point_3(0.3f, 0.2f, 0.3f),         // interior query points
        Point_3(0.1f, 0.1f, 0.1f), Point_3(0.2f, 0.5f, 0.3f),             // interior query points
        Point_3(0.0f , 0.0f, 0.5f), Point_3(0.4f, 0.4f, 0.0f),            // boundary query points
        Point_3(0.0f, 0.4f, 0.4f), Point_3(0.4f, 0.0f, 0.4f),             // boundary query points
        Point_3(0.5f, 0.5f, 0.5f), Point_3(2.0f, 0.0f, 0.0f),             // exterior query points
        Point_3(-1.0f, -1.0f, 1.0f), Point_3(0.5f, 0.5f, -2.0f)};           // exterior query point

    // Compute tetrahedra coordinates;
    std::vector<FT> coordinates;
    coordinates.reserve(queries.size()*4);

    // Output all tetrahedra coordinates.
    std::cout << std::endl << "tetrahedra coordinates (all queries): " << std::endl
    << std::endl;
    for (std::size_t i = 0; i < coordinates.size(); i += 3){
        std::cout <<
        coordinates[i + 0] << ", " <<
        coordinates[i + 1] << ", " <<
        coordinates[i + 2] << ", " <<
        coordinates[i + 3] << std::endl;
    }
    std::cout << std::endl;

    // Get a tuple of triangle coordinates for all query points
    for(std::size_t i = 0; i < queries.size(); i++){
        const auto tuple =
        CGAL::Barycentric_coordinates::tetrahedron_coordinates_in_tuple(p0, p1, p2, p3, queries[i]);

        std::cout << "tetrahedra coordinates (query " << i << "): " <<
            std::get<0>(tuple) << " " << std::get<1>(tuple) << " " <<
            std::get<2>(tuple) << " " << std::get<3>(tuple) << std::endl;
    }

  return EXIT_SUCCESS;

}
