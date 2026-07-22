#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Timer = CGAL::Real_timer;
using Vertices = std::vector<Point_3>;
using SM = CGAL::Surface_mesh<Point_3>;

using WPC3 = CGAL::Barycentric_coordinates::Wachspress_coordinates_3<SM, Kernel>;
using MVC3 = CGAL::Barycentric_coordinates::Mean_value_coordinates_3<SM, Kernel>;
using DHC3 = CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<SM, Kernel>;

template<typename COORD>
double benchmark_overload(const FT scale){

  const std::size_t number_of_runs = 10;

  // Cube
  SM cube;
  Vertices cube_coords;

  cube_coords = {Point_3(1.0, 0.0, 0.0), Point_3(1.0, 1.0, 0.0),
                Point_3(0.0, 1.0, 0.0), Point_3(0.0, 0.0, 0.0),
                Point_3(0.0, 0.0, 1.0), Point_3(1.0, 0.0, 1.0),
                Point_3(1.0, 1.0, 1.0), Point_3(0.0, 1.0, 1.0)};

  CGAL::make_hexahedron(cube_coords[0], cube_coords[1], cube_coords[2], cube_coords[3],
                        cube_coords[4], cube_coords[5], cube_coords[6], cube_coords[7], cube);

  CGAL::Polygon_mesh_processing::triangulate_faces(cube);

  COORD bar_cube(cube);

  const FT step  = FT(1) / scale;
  const FT limit = step*scale;

  std::vector<FT> bar_coordinates_cube;
  bar_coordinates_cube.resize(8);

  Timer timer;
  double time = 0.0;

  // Sample interior points
  for(std::size_t i = 0; i < number_of_runs; i++){

    timer.start();
    for(FT x = step; x < limit - step; x += step){
      for(FT y = step; y < limit - step; y += step){
        for(FT z = step; z < limit - step; z += step){

          const Point_3 query(x, y, z);
          bar_cube(query, bar_coordinates_cube.begin());
        }
      }
    }
    timer.stop();
    time += timer.time();
    timer.reset();
  }

  const double mean_time =
    time / static_cast<double>(number_of_runs);

  return mean_time;
}

int main(){

  std::ofstream wp_file("wp_bench.txt");
  std::ofstream dh_file("dh_bench.txt");
  std::ofstream mv_file("mv_bench.txt");

  std::ofstream num_file("num_bench.txt");

  for(std::size_t i = 5; i <= 100; i += 5){

    wp_file << benchmark_overload<WPC3>(i) << std::endl;
    dh_file << benchmark_overload<DHC3>(i) << std::endl;
    mv_file << benchmark_overload<MVC3>(i) << std::endl;

    num_file << i << std::endl;

    std::cout << i << "\% complete" << std::endl;
  }

  wp_file.close();
  dh_file.close();
  mv_file.close();
  num_file.close();

  return EXIT_SUCCESS;
}
