#include "test_util.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/internal/implicit_shapes_helper.h>
#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/IO/Image_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Image_3.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <deque>
#include <iostream>

namespace IS = CGAL::Isosurfacing;

template <typename K>
void test_cube()
{
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Vector = typename K::Vector_3;

  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Value_function_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

  auto implicit_function = [](const Point& q) -> FT
  {
  // ---
    auto b1 = [](const Point& q) { return IS::Shapes::box<K>(Point(-1,-1,-1), Point(1,1,1), q); };
    auto s1 = [](const Point& q) { return IS::Shapes::sphere<K>(Point(0,0,0), 1, q); };
    auto cyl1 = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(0,0,1), 0.5, q); };
    auto cyl2 = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(0,1,0), 0.5, q); };
    auto cyl3 = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(1,0,0), 0.5, q); };

    auto b1_is1 = [&](const Point& q) { return IS::Shapes::shape_intersection<K>(b1, s1, q); };
    auto cyl1_cyl2 = [&](const Point& q) { return IS::Shapes::shape_union<K>(cyl1, cyl2, q); };
    auto cyl1_cyl2_cyl3 = [&](const Point& q) { return IS::Shapes::shape_union<K>(cyl1_cyl2, cyl3, q); };

    auto b1_ms1_mcyl1_mcyl2 = [&](const Point& q) { return IS::Shapes::shape_difference<K>(b1_is1, cyl1_cyl2_cyl3, q); };

    return b1_ms1_mcyl1_mcyl2(q);
  };

  const CGAL::Bbox_3 bbox = {-2., -2., -2., 2., 2., 2.};
  const Vector spacing { 0.05, 0.05, 0.05 };
  const FT isovalue = 0;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring (Implicit function) with isovalue = " << isovalue << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  Grid grid { bbox, spacing };
  Values values {implicit_function , grid };
  Domain domain { grid, values };

  Point_range points;
  Polygon_range polygons;
  IS::marching_cubes<CGAL::Sequential_tag>(domain, isovalue, points, polygons, CGAL::parameters::use_topologically_correct_marching_cubes(true));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << polygons.size() << std::endl;

  CGAL::IO::write_polygon_soup("MC.off", points, polygons, CGAL::parameters::stream_precision(17));
}

int main(int, char**)
{
  test_cube<CGAL::Simple_cartesian<double> >();
  //test_cube<CGAL::Epeck >();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
