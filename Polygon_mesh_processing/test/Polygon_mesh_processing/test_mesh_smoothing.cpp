#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>

#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point;

typedef CGAL::Surface_mesh<Point>                             SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Mesh>
void read_mesh(const char* filename,
               Mesh& mesh)
{
  std::ifstream input(filename);

  if (!input || !(input >> mesh))
  {
    std::cerr << "Error: can not read file.";
    std::exit(1);
  }
}

template <typename Mesh>
void test_smoothing(const char* filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::smooth(mesh);
  PMP::smooth(mesh, CGAL::parameters::number_of_iterations(10));
}

template <typename Mesh>
void test_angle_smoothing(const char* filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::smooth_angles(mesh);
  PMP::smooth_angles(mesh, CGAL::parameters::number_of_iterations(10));
}

template <typename Mesh>
void test_area_smoothing(const char* filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::smooth_areas(mesh);
  PMP::smooth_areas(mesh, CGAL::parameters::number_of_iterations(10));
}

template <typename Mesh>
void test_angle_smoothing_without_projection(const char* filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::smooth_angles(mesh, CGAL::parameters::do_project(false));
}

template <typename Mesh>
void test_area_smoothing_without_projection(const char* filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  PMP::smooth_areas(mesh, CGAL::parameters::do_project(false));
}

template <typename Mesh>
void test_constrained_vertices(const char* filename)
{
  Mesh mesh;
  read_mesh(filename, mesh);

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, mesh);

  double x_init, y_init, z_init;
  std::set<vertex_descriptor> selected_vertices;
  for(vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      selected_vertices.insert(v);
      x_init = get(vpmap, v).x();
      y_init = get(vpmap, v).y();
      z_init = get(vpmap, v).z();
    }
  }

  CGAL::Boolean_property_map<std::set<vertex_descriptor> > vcmap(selected_vertices);

  PMP::smooth_angles(mesh, CGAL::parameters::vertex_is_constrained_map(vcmap));
  PMP::smooth_areas(mesh, CGAL::parameters::vertex_is_constrained_map(vcmap));

  for(vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      assert(x_init == get(vpmap, v).x());
      assert(y_init == get(vpmap, v).y());
      assert(z_init == get(vpmap, v).z());
    }
  }
}

int main(int /*argc*/, char** /*argv*/)
{
  const char* filename_elephant = "data/elephant.off";
  const char* filename_mannequin = "data/mannequin-devil.off";

  // test with Surface_mesh
  test_smoothing<SurfaceMesh>(filename_elephant);
  test_angle_smoothing<SurfaceMesh>(filename_elephant);
  test_area_smoothing<SurfaceMesh>(filename_mannequin);
  test_angle_smoothing_without_projection<SurfaceMesh>(filename_elephant);
  test_area_smoothing_without_projection<SurfaceMesh>(filename_mannequin);
  test_constrained_vertices<SurfaceMesh>(filename_elephant);

  // test with Polyhedron
  test_smoothing<Polyhedron>(filename_elephant);
  test_angle_smoothing<Polyhedron>(filename_elephant);
  test_area_smoothing<Polyhedron>(filename_mannequin);
  test_angle_smoothing_without_projection<Polyhedron>(filename_elephant);
  test_area_smoothing_without_projection<Polyhedron>(filename_mannequin);
  test_constrained_vertices<Polyhedron>(filename_mannequin);

  return EXIT_SUCCESS;
}
