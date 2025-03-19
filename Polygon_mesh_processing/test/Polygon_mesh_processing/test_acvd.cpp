#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/approximated_centroidal_Voronoi_diagram_remeshing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/subdivision_method_3.h>

#include <iostream>
#include <filesystem>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

namespace params = CGAL::parameters;

void run_test(std::string fname, std::size_t genus)
{
  Mesh mesh;
  CGAL::IO::read_polygon_mesh(fname, mesh);

  PMP::triangulate_faces(mesh);
  CGAL::Subdivision_method_3::Loop_subdivision(mesh, params::number_of_iterations(5));

  Mesh ref=mesh;
  CGAL::Bbox_3 bb_ref = PMP::bbox(ref);

  PMP::approximated_centroidal_Voronoi_diagram_remeshing(mesh, 1000);

  CGAL::Bbox_3 bb = PMP::bbox(ref);

  assert(-(vertices(mesh).size()-edges(mesh).size()+faces(mesh).size()-2)/2==genus);
  assert(vertices(mesh).size()==1000);

  double dx = bb.xmax()-bb.xmin(),
         dy = bb.ymax()-bb.ymin(),
         dz = bb.zmax()-bb.zmin();

  assert(std::abs(bb.xmax()-bb_ref.xmax()) < 0.01 * dx);
  assert(std::abs(bb.ymax()-bb_ref.ymax()) < 0.01 * dy);
  assert(std::abs(bb.zmax()-bb_ref.zmax()) < 0.01 * dz);

std::ofstream("/tmp/out_"+std::to_string(genus)+".off") << mesh;

  PMP::approximated_centroidal_Voronoi_diagram_remeshing(mesh, 4);
  assert(-(vertices(mesh).size()-edges(mesh).size()+faces(mesh).size()-2)/2==genus);

std::ofstream("/tmp/min_out_"+std::to_string(genus)+".off") << mesh;
}

int main()
{
  run_test(CGAL::data_file_path("meshes/tetrahedron.off"), 0);
  run_test(CGAL::data_file_path("meshes/torus_quad.off"), 1);
  run_test(CGAL::data_file_path("meshes/double-torus-example.off"), 2);

  return 0;
}

