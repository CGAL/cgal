#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/curvature.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/graph/graph_traits.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define ABS_ERROR 1e-6

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef K::FT FT;
typedef CGAL::Surface_mesh<K::Point_3> SMesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

struct Average_test_info
{
  FT mean_curvature_avg;
  FT gaussian_curvature_avg;
  FT tolerance = 0.9;

  Average_test_info(FT mean_curvature_avg,
                    FT gaussian_curvature_avg)
    : mean_curvature_avg(mean_curvature_avg),
      gaussian_curvature_avg(gaussian_curvature_avg)
  { }
};

bool passes_comparison(FT result, FT expected, FT tolerance)
{
  std::cout << "result: " << result << std::endl;
  std::cout << "expected: " << expected << std::endl;

  if(abs(expected) < ABS_ERROR && abs(result) < ABS_ERROR)
    return true; // expected 0, got 0
  else if (abs(expected) < ABS_ERROR)
    return false; // expected 0, got non-0

  return (std::min)(result, expected) / (std::max)(result, expected) > tolerance;
}

template <typename TriangleMesh>
void test_curvatures(std::string mesh_path,
                     Average_test_info test_info)
{
  std::cout << "test discrete curvatures of " << mesh_path << std::endl;
  std::cout << "mesh type: " << typeid(mesh_path).name() << std::endl;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  TriangleMesh tmesh;
  const std::filesystem::path filename = CGAL::data_file_path(mesh_path);

  if(!CGAL::IO::read_polygon_mesh(filename, tmesh) || faces(tmesh).size() == 0)
  {
    std::cerr << "Invalid input file." << std::endl;
    std::exit(1);
  }

  typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<FT>>::type
    mean_curvature_map = get(CGAL::dynamic_vertex_property_t<FT>(), tmesh),
    gaussian_curvature_map = get(CGAL::dynamic_vertex_property_t<FT>(), tmesh);

  PMP::discrete_mean_curvatures(tmesh, mean_curvature_map);
  PMP::discrete_Gaussian_curvatures(tmesh, gaussian_curvature_map);

  FT mean_curvature_avg = 0, gaussian_curvature_avg = 0;
  for(vertex_descriptor v : vertices(tmesh))
  {
    mean_curvature_avg += get(mean_curvature_map, v);
    gaussian_curvature_avg += get(gaussian_curvature_map, v);
  }

  mean_curvature_avg /= vertices(tmesh).size();
  gaussian_curvature_avg /= vertices(tmesh).size();

  std::cout << "checking mean curvature..." << std::endl;
  assert(passes_comparison(mean_curvature_avg, test_info.mean_curvature_avg, test_info.tolerance));

  std::cout << "checking Gaussian curvature..." << std::endl;
  assert(passes_comparison(gaussian_curvature_avg, test_info.gaussian_curvature_avg, test_info.tolerance));
}

template <typename PolygonMesh>
void test_angle_sums(const std::string mesh_path,
                     const std::vector<FT>& expected_values)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  PolygonMesh pmesh;
  const std::filesystem::path filename = CGAL::data_file_path(mesh_path);

  if(!CGAL::IO::read_polygon_mesh(filename, pmesh) || faces(pmesh).size() == 0)
  {
    std::cerr << "Invalid input file." << std::endl;
    std::exit(1);
  }

  std::size_t pos = 0;
  for(vertex_descriptor v : vertices(pmesh))
  {
    FT angle_sum = PMP::angle_sum(v, pmesh,
                                  CGAL::parameters::geom_traits(K())
                                                   .vertex_point_map(get(CGAL::vertex_point, pmesh)));
    assert(passes_comparison(angle_sum, expected_values[pos++], 0.9));
  }
}

// Builds a triangle fan made of `n` triangles around a central vertex (the first vertex)
// located at the origin. The i-th vertex of the ring is at angle 2*i*pi/n, at distance `radius(i)`
// from the z axis, and at height `height(theta)`.
template <typename PolygonMesh, typename Kernel, typename RadiusFunctor, typename HeightFunctor>
void make_fan(PolygonMesh& pmesh, const std::size_t n,
              const RadiusFunctor& radius, const HeightFunctor& height)
{
  std::vector<typename Kernel::Point_3> points;
  std::vector<std::array<std::size_t, 3> > faces;

  points.emplace_back(0, 0, 0);
  for(std::size_t i=0; i<n; ++i)
  {
    const double theta = 2 * CGAL_PI * double(i) / double(n);
    const double r = radius(i);
    points.emplace_back(r * std::cos(theta), r * std::sin(theta), height(theta));
    faces.push_back({0, i+1, (i+1)%n + 1});
  }

  PMP::polygon_soup_to_polygon_mesh(points, faces, pmesh);
  assert(is_valid_polygon_mesh(pmesh));
}

template <typename PolygonMesh, typename Kernel = K>
void test_saddle_fans()
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  auto unit_radius = [](std::size_t) { return 1.; };
  auto irregular_radius = [](std::size_t i) { return 1. + 0.37 * double((7*i) % 5); };

  const std::size_t n = 12;

  // flat
  {
    PolygonMesh pmesh;
    make_fan<PolygonMesh, Kernel>(pmesh, n, unit_radius, [](double) { return 0.; });
    assert(PMP::is_saddle_vertex(*vertices(pmesh).begin(), pmesh));
  }
  {
    PolygonMesh pmesh;
    make_fan<PolygonMesh, Kernel>(pmesh, n, irregular_radius, [](double) { return 0.; });
    assert(PMP::is_saddle_vertex(*vertices(pmesh).begin(), pmesh));
  }

  // convex (cone)
  {
    PolygonMesh pmesh;
    make_fan<PolygonMesh, Kernel>(pmesh, n, unit_radius, [](double) { return -0.5; });
    assert(!PMP::is_saddle_vertex(*vertices(pmesh).begin(), pmesh));
  }
  {
    PolygonMesh pmesh; // almost flat cone
    make_fan<PolygonMesh, Kernel>(pmesh, n, unit_radius, [](double) { return 1e-3; });
    assert(!PMP::is_saddle_vertex(*vertices(pmesh).begin(), pmesh));
  }

  // saddle and monkey saddle
  {
    PolygonMesh pmesh;
    make_fan<PolygonMesh, Kernel>(pmesh, n, unit_radius, [](double t) { return 0.5 * std::cos(2*t); });
    assert(PMP::is_saddle_vertex(*vertices(pmesh).begin(), pmesh));
  }
  {
    PolygonMesh pmesh;
    make_fan<PolygonMesh, Kernel>(pmesh, n, irregular_radius, [](double t) { return 0.5 * std::cos(3*t); });
    const vertex_descriptor v = *vertices(pmesh).begin();
    assert(PMP::is_saddle_vertex(v, pmesh));
    assert(PMP::is_saddle_vertex(v, pmesh, CGAL::parameters::geom_traits(Kernel())
                                                            .vertex_point_map(get(CGAL::vertex_point, pmesh))));
  }
}

template <typename PolygonMesh>
void test_saddle_vertices_of_closed_mesh(const std::string mesh_path)
{
  std::cout << "test saddle vertices of " << mesh_path << std::endl;

  PolygonMesh pmesh;
  const std::filesystem::path filename = CGAL::data_file_path(mesh_path);

  if(!CGAL::IO::read_polygon_mesh(filename, pmesh) || faces(pmesh).size() == 0)
  {
    std::cerr << "Invalid input file." << std::endl;
    std::exit(1);
  }

  // all vertices of these convex polyhedra are strictly convex
  for(auto v : vertices(pmesh))
    assert(!PMP::is_saddle_vertex(v, pmesh));
}

template <typename PolygonMesh>
void test_saddle_vertices_of_quad_grid()
{
  // 3x3 grid of planar quads: the 4 interior vertices are flat, the others are on the border
  std::vector<K::Point_3> points;
  std::vector<std::array<std::size_t, 4> > faces;
  for(std::size_t j=0; j<4; ++j)
    for(std::size_t i=0; i<4; ++i)
      points.emplace_back(double(i), double(j), 0);
  for(std::size_t j=0; j<3; ++j)
    for(std::size_t i=0; i<3; ++i)
      faces.push_back({4*j+i, 4*j+i+1, 4*(j+1)+i+1, 4*(j+1)+i});

  PolygonMesh pmesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, pmesh);

  std::size_t nb_saddles = 0;
  for(auto v : vertices(pmesh))
  {
    if(PMP::is_saddle_vertex(v, pmesh))
    {
      assert(!is_border(v, pmesh));
      ++nb_saddles;
    }
  }
  assert(nb_saddles == 4);
}

void test_saddle_vertices()
{
  std::cout << "test saddle vertices" << std::endl;

  test_saddle_fans<SMesh>();
  test_saddle_fans<Polyhedron>();
  test_saddle_fans<CGAL::Surface_mesh<EPECK::Point_3>, EPECK>();

  test_saddle_vertices_of_closed_mesh<SMesh>("meshes/regular_tetrahedron.off");
  test_saddle_vertices_of_closed_mesh<Polyhedron>("meshes/regular_tetrahedron.off");
  test_saddle_vertices_of_closed_mesh<SMesh>("meshes/cube_quad.off");
  test_saddle_vertices_of_closed_mesh<Polyhedron>("meshes/cube_quad.off");
  test_saddle_vertices_of_closed_mesh<SMesh>("meshes/sphere.off");

  test_saddle_vertices_of_quad_grid<SMesh>();
  test_saddle_vertices_of_quad_grid<Polyhedron>();
}

int main(int, char**)
{
  test_saddle_vertices();

  // testing on a simple sphere(r = 0.5), on both Polyhedron & SurfaceMesh:
  // Expected: Mean Curvature = 2, Gaussian Curvature = 4
  test_curvatures<Polyhedron>("meshes/sphere.off", Average_test_info(2, 4));
  test_curvatures<SMesh>("meshes/sphere.off", Average_test_info(2, 4));

  // testing on a simple sphere(r = 10), on both Polyhedron & SurfaceMesh:
  // Expected: Mean Curvature = 0.1, Gaussian Curvature = 0.01
  test_curvatures<Polyhedron>("meshes/sphere966.off", Average_test_info(0.1, 0.01));
  test_curvatures<SMesh>("meshes/sphere966.off", Average_test_info(0.1, 0.01));

  // testing on a simple half cylinder(r = 1), on both Polyhedron & SurfaceMesh:
  // Expected: Mean Curvature = 0.5, Gaussian Curvature = 0
  // To be tested once the discrete curvatures are well defined for boundary vertices
  // test_curvatures<Polyhedron>("meshes/cylinder.off", Average_test_info(0.5, 0));
  // test_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0));

  test_angle_sums<Polyhedron>("meshes/quad.off", std::vector<FT>(4, 90));
  test_angle_sums<SMesh>("meshes/quad.off", std::vector<FT>(4, 90));

  test_angle_sums<Polyhedron>("meshes/regular_tetrahedron.off", std::vector<FT>(4, 180));
  test_angle_sums<SMesh>("meshes/regular_tetrahedron.off", std::vector<FT>(4, 180));

  test_angle_sums<Polyhedron>("meshes/cube_quad.off", std::vector<FT>(8, 270));
  test_angle_sums<SMesh>("meshes/cube_quad.off", std::vector<FT>(8, 270));

  test_angle_sums<Polyhedron>("meshes/cube_poly.off", std::vector<FT>(8, 270));
  test_angle_sums<SMesh>("meshes/cube_poly.off", std::vector<FT>(8, 270));

  std::cout << "Done." << std::endl;
  return EXIT_SUCCESS;
}
