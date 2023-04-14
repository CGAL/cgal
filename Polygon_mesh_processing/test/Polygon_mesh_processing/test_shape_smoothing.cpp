#define CGAL_PMP_SMOOTHING_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

#include <CGAL/utility.h>
#include <CGAL/use.h>

#include <fstream>
#include <iostream>
#include <set>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef Kernel::FT                                                          FT;
typedef Kernel::Point_3                                                     Point;
typedef CGAL::Surface_mesh<Point>                                           SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>        Mesh_with_id;

namespace PMP = CGAL::Polygon_mesh_processing;

bool equal_doubles(double d1, double d2, double e)
{
  return (d1 > d2 - e) && (d1 < d2 + e);
}

template <typename Mesh>
void test_implicit_constrained_devil(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::cout << "-- test_implicit_constrained_devil --" << std::endl;
#endif

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, mesh);

  // max 'z' is 20 in the devil
  std::set<vertex_descriptor> selected_vertices;
  for(vertex_descriptor v : vertices(mesh))
  {
    const double z = get(vpmap, v).z();
    if(is_border(v, mesh) || z > 19.0)
      selected_vertices.insert(v);
  }

  std::cout << selected_vertices.size() << " constrained vertices" << std::endl;

  CGAL::Boolean_property_map<std::set<vertex_descriptor> > vcmap(selected_vertices);

  std::vector<Point> fixed_points(selected_vertices.size());
  std::size_t i = 0;
  for(vertex_descriptor v : selected_vertices)
    fixed_points[i++] = get(vpmap, v);

  const double time_step = 1.0;
  PMP::smooth_shape(mesh, time_step, CGAL::parameters::vertex_is_constrained_map(vcmap)
                                                      .number_of_iterations(2));

  i = 0;
  for(vertex_descriptor v : selected_vertices)
  {
    const Point p = get(vpmap, v);
    CGAL_USE(p);
    assert(equal_doubles(p.x(), fixed_points[i].x(), 1e-10));
    assert(equal_doubles(p.y(), fixed_points[i].y(), 1e-10));
    assert(equal_doubles(p.z(), fixed_points[i].z(), 1e-10));
    ++i;
  }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::ofstream out("output_implicit_constrained_devil.off");
  out << mesh;
  out.close();
#endif
}

template <typename Mesh>
void test_implicit_constrained_elephant(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::cout << "-- test_implicit_constrained_elephant --" << std::endl;
#endif

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, mesh);

  std::set<vertex_descriptor> selected_vertices;
  for(vertex_descriptor v : vertices(mesh))
  {
    const double y = get(vpmap, v).z();
    if(y > 0.14)
      selected_vertices.insert(v);
  }

  CGAL::Boolean_property_map<std::set<vertex_descriptor> > vcmap(selected_vertices);

  std::vector<Point> fixed_points(selected_vertices.size());
  std::size_t i = 0;
  for(vertex_descriptor v : selected_vertices)
    fixed_points[i++] = get(vpmap, v);

  const double time_step = 1.0;
  PMP::smooth_shape(mesh, time_step,
                    CGAL::parameters::vertex_is_constrained_map(vcmap)
                                     .number_of_iterations(1));

  i = 0;
  for(vertex_descriptor v : selected_vertices)
  {
    const Point p = get(vpmap, v);
    CGAL_USE(p);
    assert(equal_doubles(p.x(), fixed_points[i].x(), 1e-10));
    assert(equal_doubles(p.y(), fixed_points[i].y(), 1e-10));
    assert(equal_doubles(p.z(), fixed_points[i].z(), 1e-10));
    ++i;
  }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::ofstream out("output_implicit_constrained_elephant.off");
  out << mesh;
  out.close();
#endif
}

template <typename Mesh>
void test_implicit_unscaled_elephant(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::cout << "-- test_implicit_unscaled_elephant --" << std::endl;
#endif

  const FT ivol = PMP::volume(mesh);
  std::cout << "Input volume is " << ivol << std::endl;

  Mesh mesh_cpy(mesh);
  const double time_step = 0.001;
  PMP::smooth_shape(mesh_cpy, time_step, CGAL::parameters::number_of_iterations(5).do_scale(true));

  FT ovol = PMP::volume(mesh_cpy);
  std::cout << "With scaling, output volume is " << ovol << std::endl;
  assert(equal_doubles(ivol, ovol, 1e-10));

  PMP::smooth_shape(mesh, time_step, CGAL::parameters::number_of_iterations(5).do_scale(false));
  ovol = PMP::volume(mesh);
  std::cout << "Without scaling, output volume is " << ovol << std::endl;
}

template <typename Mesh>
void test_curvature_flow_time_step(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::cout << "-- test_curvature_flow_time_step --" << std::endl;
#endif

  const double time_step = 1e-15;
  PMP::smooth_shape(mesh, time_step);

#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::ofstream out("output_devil_time_step.off");
  out << mesh;
  out.close();
#endif
}

template <typename Mesh>
void test_curvature_flow(Mesh mesh)
{
#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::cout << "-- test_curvature_flow --" << std::endl;
#endif

  const double time_step = 1.0;
  PMP::smooth_shape(mesh, time_step);

#ifdef CGAL_PMP_SMOOTHING_DEBUG
  std::ofstream out("output_precision_elephant.off");
  out << mesh;
  out.close();
#endif
}

int main(int, char**)
{
  const std::string filename_devil = CGAL::data_file_path("meshes/mannequin-devil.off");
  const std::string filename_elephant = CGAL::data_file_path("meshes/elephant.off");

  std::ifstream input1(filename_devil);
  SurfaceMesh mesh_devil;
  if(!input1 || !(input1 >> mesh_devil))
  {
    std::cerr << "Error: cannot read file " << filename_devil << std::endl;
    return EXIT_FAILURE;
  }
  input1.close();

  std::ifstream input2(filename_elephant);
  SurfaceMesh mesh_elephant;
  if(!input2 || !(input2 >> mesh_elephant))
  {
    std::cerr << "Error: cannot read file " << filename_elephant << std::endl;
    return EXIT_FAILURE;
  }
  input2.close();

  test_curvature_flow_time_step<SurfaceMesh>(mesh_devil);
  test_curvature_flow<SurfaceMesh>(mesh_elephant);
  test_implicit_constrained_elephant<SurfaceMesh>(mesh_elephant);
  test_implicit_constrained_devil<SurfaceMesh>(mesh_devil);
  test_implicit_unscaled_elephant<SurfaceMesh>(mesh_elephant);

  input1.open(filename_devil);
  Mesh_with_id pl_mesh_devil;
  if(!input1 || !(input1 >> pl_mesh_devil)){
    std::cerr << "Error: cannot read file " << filename_devil << std::endl;
    return EXIT_FAILURE;
  }
  input1.close();

  // Polyhedron

  input2.open(filename_elephant);
  Mesh_with_id pl_mesh_elephant;
  if(!input2 || !(input2 >> pl_mesh_elephant))
  {
    std::cerr << "Error: cannot read file " << filename_elephant << std::endl;
    return EXIT_FAILURE;
  }
  input2.close();

  set_halfedgeds_items_id(pl_mesh_devil);
  set_halfedgeds_items_id(pl_mesh_elephant);

  test_curvature_flow_time_step<Mesh_with_id>(pl_mesh_devil);
  test_curvature_flow<Mesh_with_id>(pl_mesh_elephant);
  test_implicit_constrained_elephant<Mesh_with_id>(pl_mesh_elephant);
  test_implicit_constrained_devil<Mesh_with_id>(pl_mesh_devil);
  test_implicit_unscaled_elephant<Mesh_with_id>(pl_mesh_elephant);

  return EXIT_SUCCESS;
}
