#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>

#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point;

typedef CGAL::Surface_mesh<Point>                             SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Mesh>
void read_mesh(const std::string filename, Mesh& mesh)
{
  std::ifstream input(filename);

  if (!input || !(input >> mesh))
  {
    std::cerr << "Error: can not read file.";
    std::exit(1);
  }
}

template <typename Mesh>
void test_smoothing(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::angle_and_area_smoothing(mesh);
  PMP::angle_and_area_smoothing(mesh, CGAL::parameters::number_of_iterations(10));
}

template <typename Mesh>
void test_angle_smoothing(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::angle_and_area_smoothing(mesh);
  PMP::angle_and_area_smoothing(mesh, CGAL::parameters::number_of_iterations(10)
                                                       .use_area_smoothing(false));
}

template <typename Mesh>
void test_area_smoothing(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::angle_and_area_smoothing(mesh);
  PMP::angle_and_area_smoothing(mesh, CGAL::parameters::number_of_iterations(10)
                                                       .use_angle_smoothing(false));
}

template <typename Mesh>
void test_angle_smoothing_without_projection(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::angle_and_area_smoothing(mesh, CGAL::parameters::do_project(false)
                                                       .use_area_smoothing(false));
}

template <typename Mesh>
void test_area_smoothing_without_projection(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::angle_and_area_smoothing(mesh, CGAL::parameters::do_project(false)
                                                       .use_angle_smoothing(false));
}

template<typename Mesh>
void test_tangential_relaxation(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);
  PMP::tangential_relaxation(
    vertices(mesh),
    mesh,
    CGAL::parameters::number_of_iterations(4)
                     .relax_constraints(false));
  PMP::tangential_relaxation(vertices(mesh), mesh);
  PMP::tangential_relaxation(mesh);

  // test allow_move_functor
  typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Point> >::type vpm =
    get(CGAL::dynamic_vertex_property_t<Point>(), mesh);
  for (auto v : vertices(mesh))
    put(vpm, v, get(CGAL::vertex_point, mesh, v));
  auto no_move = [](typename boost::graph_traits<Mesh>::vertex_descriptor, Point, Point)
  {
    return false;
  };
  PMP::tangential_relaxation(vertices(mesh), mesh, CGAL::parameters::allow_move_functor(no_move)
                                                        .vertex_point_map(vpm));
  for (auto v : vertices(mesh))
  {
    assert(get(vpm, v) == get(CGAL::vertex_point, mesh, v));
  }
}

template <typename Mesh>
void test_constrained_vertices(const std::string filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, mesh);

  std::set<vertex_descriptor> selected_vertices;
  std::map<vertex_descriptor, Point> initial_positions;
  for(vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      selected_vertices.insert(v);
      initial_positions[v] = get(vpmap, v);
    }
  }

  CGAL::Boolean_property_map<std::set<vertex_descriptor> > vcmap(selected_vertices);

  PMP::angle_and_area_smoothing(mesh, CGAL::parameters::vertex_is_constrained_map(vcmap));

  for(vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
      assert(initial_positions.at(v) == get(vpmap, v));
  }
}

int main(int /*argc*/, char** /*argv*/)
{
  const std::string filename_elephant = CGAL::data_file_path("meshes/elephant.off");
  const std::string filename_mannequin = CGAL::data_file_path("meshes/mannequin-devil.off");

  std::cout << "Test files: " << filename_elephant << " " << filename_mannequin << std::endl;

  // test with Surface_mesh
  test_smoothing<SurfaceMesh>(filename_elephant);
  test_angle_smoothing<SurfaceMesh>(filename_elephant);
  test_area_smoothing<SurfaceMesh>(filename_mannequin);
  test_angle_smoothing_without_projection<SurfaceMesh>(filename_elephant);
  test_area_smoothing_without_projection<SurfaceMesh>(filename_mannequin);
  test_constrained_vertices<SurfaceMesh>(filename_elephant);
  test_tangential_relaxation<SurfaceMesh>(filename_elephant);

  // test with Polyhedron
  test_smoothing<Polyhedron>(filename_elephant);
  test_angle_smoothing<Polyhedron>(filename_elephant);
  test_area_smoothing<Polyhedron>(filename_mannequin);
  test_angle_smoothing_without_projection<Polyhedron>(filename_elephant);
  test_area_smoothing_without_projection<Polyhedron>(filename_mannequin);
  test_constrained_vertices<Polyhedron>(filename_mannequin);
  test_tangential_relaxation<Polyhedron>(filename_elephant);

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
