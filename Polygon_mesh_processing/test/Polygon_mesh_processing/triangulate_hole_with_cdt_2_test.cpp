#include <set>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>

#define CGAL_NO_CDT_2_WARNING

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

template<
class PolygonMesh,
class Halfedge_handle>
void detect_borders(
  PolygonMesh& pmesh,
  std::vector<Halfedge_handle>& borders) {

  typedef CGAL::Halfedge_around_face_circulator<PolygonMesh>
    Halfedge_around_facet_circulator;

  borders.clear();
  std::set<Halfedge_handle> border_map;
  for (Halfedge_handle h : halfedges(pmesh)) {
    if (
      face(h, pmesh) == boost::graph_traits<PolygonMesh>::null_face() &&
      border_map.find(h) == border_map.end()) {

      borders.push_back(h);
      Halfedge_around_facet_circulator hf_around_facet(h, pmesh);
      Halfedge_around_facet_circulator done(hf_around_facet);

      do {
        assert(border_map.insert(*hf_around_facet).second); // is insertion ok?
      } while (++hf_around_facet != done);
    }
  }
}

// This test is inspired by the issue: https://github.com/CGAL/cgal/issues/4464.
template<
typename PolygonMesh,
typename GeomTraits>
void test_triangulate_hole_with_cdt_2(
  const std::string kernel_name,
  int argc, char **argv,
  const std::string file_name,
  const std::size_t num_borders,
  const std::size_t num_patch_faces,
  const bool verbose) {

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor Face_handle;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_handle;

  // Reading the file.
  if (verbose) {
    std::cout << "test with the " << kernel_name << " kernel:" << std::endl;
  }
  PolygonMesh pmesh;
  std::string path = "data/" + file_name + ".off";
  std::ifstream in(path.c_str(), std::ios_base::in);
  CGAL::IO::set_ascii_mode(in);
  CGAL::IO::read_OFF(in, pmesh);
  in.close();
  if (verbose) {
    std::cout << "* finished reading the file" << std::endl;
  }

  // Detecting the hole borders.
  std::vector<Halfedge_handle> borders;
  detect_borders(pmesh, borders);
  if (verbose) {
    std::cout << "* number of detected borders: " <<
      borders.size() << std::endl;
  }
  assert(borders.size() == num_borders);

  // Triangulating the holes.
  std::vector<Face_handle> patch_faces;
  for (const Halfedge_handle& h : borders) {
    patch_faces.clear();
    CGAL::Polygon_mesh_processing::triangulate_hole(
      pmesh,
      h,
      std::back_inserter(patch_faces),
      CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
        get(CGAL::vertex_point, pmesh)).
        use_2d_constrained_delaunay_triangulation(true).
        geom_traits(GeomTraits()));

    if (verbose) {
      std::cout << "* number of faces in the constructed patch: " <<
        patch_faces.size() << std::endl;
    }
    assert(patch_faces.size() == num_patch_faces);
  }
  assert(pmesh.is_valid() && is_closed(pmesh));
  assert(CGAL::Polygon_mesh_processing::is_outward_oriented(pmesh,
    CGAL::parameters::all_default()));

  // Writing the file.
  if (verbose) {
    path = "";
    if (argc > 1) path = std::string(argv[1]);
    path += "4464_" + file_name + ".off";
    std::ofstream out(path.c_str(), std::ios_base::out);
    CGAL::IO::set_ascii_mode(out);
    CGAL::IO::write_OFF(out, pmesh);
    out.close();
    std::cout << "* finished writing the file" << std::endl;
  }
}

int main(int argc, char **argv) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel EI;
  typedef CGAL::Exact_predicates_exact_constructions_kernel EE;

  typedef CGAL::Surface_mesh<EI::Point_3> Surface_mesh_EI;
  typedef CGAL::Polyhedron_3<EE> Polyhedron_3_EE;

  // Checking on a data file with two planar, simple, and horizontal holes.
  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, "w_horizontal_hole", 2, 25, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, "w_horizontal_hole", 2, 25, false);
  std::cout <<
    "test_triangulate_hole_with_cdt_2: horizontal planar hole SUCCESS" << std::endl;

  // Checking on a data file with two planar, simple, and orthogonal holes.
  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, "w_orthogonal_hole", 2, 2, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, "w_orthogonal_hole", 2, 2, false);
  std::cout <<
    "test_triangulate_hole_with_cdt_2: orthogonal planar hole SUCCESS" << std::endl;

  // Checking on a data file with two planar, simple, and horizontal holes.
  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, "elephant_flat_hole", 1, 17, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, "elephant_flat_hole", 1, 17, false);
  std::cout <<
    "test_triangulate_hole_with_cdt_2: near flat hole SUCCESS" << std::endl;

  // Checking on a data file with two planar, simple, and horizontal holes.
  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, "elephant_concave_hole", 1, 24, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, "elephant_concave_hole", 1, 24, false);
  std::cout <<
    "test_triangulate_hole_with_cdt_2: concave hole SUCCESS" << std::endl;

  // Checking on a data file with one simple but not a planar hole.
  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, "elephant_curved_hole", 1, 19, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, "elephant_curved_hole", 1, 19, false);
  std::cout <<
    "test_triangulate_hole_with_cdt_2: curved hole SUCCESS" << std::endl;

  // Checking on a data file with one hole that is neither simple nor planar.
  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, "elephant_complex_hole", 1, 29, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, "elephant_complex_hole", 1, 29, false);
  std::cout <<
    "test_triangulate_hole_with_cdt_2: complex hole SUCCESS" << std::endl;

  return EXIT_SUCCESS;
}
