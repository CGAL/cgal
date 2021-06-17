#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

//Typedefs
using Kernel = CGAL::Simple_cartesian<double>;

int main(){

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh =  typename CGAL::Surface_mesh<Point_3>;

  // Set cout precision
  std::cout.precision(20);

  // Regular tetrahedron
  Mesh tetrahedron;

  // Vertices
  const Point_3 p0(0.0, 0.0, 0.0);
  const Point_3 p1(1.0, 0.0, 0.0);
  const Point_3 p2(0.0, 1.0, 0.0);
  const Point_3 p3(0.0, 0.0, 1.0);

  // Construct tetrahedron
  CGAL::make_tetrahedron(p0, p1, p2, p3, tetrahedron);

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> ws(tetrahedron);

  std::vector<FT> tri_coordinates;
  std::vector<FT> wp_coordinates;

  // Tolerance
  const FT tol = FT(1.0e-13);

  // Sample points
  const FT step  = FT(1) / FT(500);
  const FT scale = FT(100);

  std::size_t count = 0;
  const FT limit = scale * step;

  for(FT x = step; x < limit; x += step){
    for(FT y = step; y < limit; y += step){
      for(FT z = step; z < FT(1) - x - y; z+= step){ // Excludes points inside faces

        const Point_3 query = Point_3(x, y, z);

        ws(query, std::back_inserter(wp_coordinates));

        CGAL::Barycentric_coordinates::tetrahedron_coordinates(p0, p1,
         p2, p3, query, std::back_inserter(tri_coordinates));

        assert(tri_coordinates[count + 0] >= FT(0) && tri_coordinates[count + 0] <= FT(1));
        assert(tri_coordinates[count + 1] >= FT(0) && tri_coordinates[count + 1] <= FT(1));
        assert(tri_coordinates[count + 2] >= FT(0) && tri_coordinates[count + 2] <= FT(1));
        assert(tri_coordinates[count + 3] >= FT(0) && tri_coordinates[count + 3] <= FT(1));

        assert(wp_coordinates[count + 0] >= FT(0) && wp_coordinates[count + 0] <= FT(1));
        assert(wp_coordinates[count + 1] >= FT(0) && wp_coordinates[count + 1] <= FT(1));
        assert(wp_coordinates[count + 2] >= FT(0) && wp_coordinates[count + 2] <= FT(1));
        assert(wp_coordinates[count + 3] >= FT(0) && wp_coordinates[count + 3] <= FT(1));

        assert(
          CGAL::abs(tri_coordinates[count + 0] - wp_coordinates[count + 0]) < tol &&
          CGAL::abs(tri_coordinates[count + 1] - wp_coordinates[count + 1]) < tol &&
          CGAL::abs(tri_coordinates[count + 2] - wp_coordinates[count + 2]) < tol &&
          CGAL::abs(tri_coordinates[count + 3] - wp_coordinates[count + 3]) < tol);

        count += 4;
      }
    }
  }

  std::cout << "test_wp_tetrahedron: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
