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

template <typename Mesh>
void test_angle_smoothing(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }
  input.close();

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh);

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      CGAL_assertion(p_c.x() == 0.7203429230262004);
      CGAL_assertion(p_c.y() == 0.5);
      CGAL_assertion(p_c.z() == 0);
      break;
    }
  }
}

template <typename Mesh>
void test_area_smoothing(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }
  input.close();

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_areas(mesh);

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      CGAL_assertion(p_c.x() == 0.6691415930575334);
      CGAL_assertion(p_c.y() == 0.5);
      CGAL_assertion(p_c.z() == 0);
      break;
    }
  }
}

template <typename Mesh>
void test_angle_smoothing_without_projection(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }
  input.close();

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh, CGAL::Polygon_mesh_processing::parameters::do_project(false));

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.59571844622769954);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 1.0652302640732678);
      break;
    }
  }
}

template <typename Mesh>
void test_area_smoothing_without_projection(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }
  input.close();

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_areas(mesh, CGAL::Polygon_mesh_processing::parameters::do_project(false));

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      CGAL_assertion(p_c.x() == 0.42183982448892759);
      CGAL_assertion(p_c.y() == 0.5);
      CGAL_assertion(p_c.z() == 0.87816017551107273);
      break;
    }
  }
}

template <typename Mesh>
void test_constrained_vertices(const char* filename)
{
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Error: can not read file.";
    return;
  }
  input.close();

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  double x_init, y_init, z_init;
  std::set<vertex_descriptor> selected_vertices;
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
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

  CGAL::Polygon_mesh_processing::smooth_angles(mesh,
        CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap));
  CGAL::Polygon_mesh_processing::smooth_areas(mesh,
        CGAL::Polygon_mesh_processing::parameters::vertex_is_constrained_map(vcmap));

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      CGAL_assertion(x_init == get(vpmap, v).x());
      CGAL_assertion(y_init == get(vpmap, v).y());
      CGAL_assertion(z_init == get(vpmap, v).z());
    }
  }
}

int main(int /*argc*/, char** /*argv*/)
{
  const char* filename_polygon = "data/simple_polygon.off";
  const char* filename_pyramid = "data/simple_pyramid.off";

  // test with Surface_mesh
  test_angle_smoothing<SurfaceMesh>(filename_polygon);
  test_area_smoothing<SurfaceMesh>(filename_polygon);
  test_constrained_vertices<SurfaceMesh>(filename_polygon);
  test_angle_smoothing_without_projection<SurfaceMesh>(filename_pyramid);
  test_area_smoothing_without_projection<SurfaceMesh>(filename_pyramid);

  // test with Polyhedron
  test_angle_smoothing<Polyhedron>(filename_polygon);
  test_area_smoothing<Polyhedron>(filename_polygon);
  test_constrained_vertices<Polyhedron>(filename_polygon);
  test_angle_smoothing_without_projection<Polyhedron>(filename_pyramid);
  test_area_smoothing_without_projection<Polyhedron>(filename_pyramid);

  return EXIT_SUCCESS;
}
