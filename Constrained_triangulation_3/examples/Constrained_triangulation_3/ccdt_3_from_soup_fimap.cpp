#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/draw_constrained_triangulation_3.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <algorithm>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;
using face_descriptor = Surface_mesh::Face_index;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  auto filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cubes_one_patch_id_per_cc.off");

  CGAL::Surface_mesh<K::Point_3> mesh;
    if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  auto fpmap = mesh.add_property_map<face_descriptor, std::size_t>("f:cc_id").first;

  std::cout << "Read " << mesh.number_of_vertices() << " vertices and "
            << mesh.number_of_faces() << " polygons" << std::endl;

  PMP::connected_components(mesh, fpmap);

  std::cout << "Number of connected components: "
            << fpmap[*std::max_element(mesh.faces_begin(), mesh.faces_end(),
                                [&](face_descriptor f1, face_descriptor f2) {
                                  return fpmap[f1] < fpmap[f2];
                                })] + 1
            << std::endl;

  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);

  auto polygon_to_patch_id = [&](std::size_t i) {
    return fpmap[*std::next(mesh.faces_begin(), i)];
  };

  auto soup_fpmap = boost::make_function_property_map<std::size_t>(polygon_to_patch_id);
  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(
      points, polygons, CGAL::parameters::plc_face_id(soup_fpmap));

  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  std::ofstream ofs(argc > 2 ? argv[2] : "out.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, ccdt);

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}
