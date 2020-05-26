#include <set>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

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
  BOOST_FOREACH(Halfedge_handle h, halfedges(pmesh)) {
    if (
      face(h, pmesh) == boost::graph_traits<PolygonMesh>::null_face() &&
      border_map.find(h) == border_map.end()) {

      borders.push_back(h);
      Halfedge_around_facet_circulator hf_around_facet(h, pmesh);
      Halfedge_around_facet_circulator done(hf_around_facet);

      do {
        const bool insertion_ok = border_map.insert(*hf_around_facet).second;
        assert(insertion_ok);
      } while (++hf_around_facet != done);
    }
  }
}

// This test is inspired by the issue: https://github.com/CGAL/cgal/issues/4464.
template<
typename PolygonMesh,
typename GT>
void test_triangulate_hole_with_cdt_2(
  const std::string kernel_name,
  int argc, char **argv, const bool save) {

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor Face_handle;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_handle;

  // Reading the file.
  std::cout << "test with the " << kernel_name << " kernel:" << std::endl;
  PolygonMesh pmesh;
  std::string path = "data/w.off";
  std::ifstream in(path.c_str(), std::ios_base::in);
  CGAL::set_ascii_mode(in);
  CGAL::read_off(in, pmesh);
  in.close();
  std::cout << "* finished reading the file" << std::endl;

  // Detecting hole borders.
  std::vector<Halfedge_handle> borders;
  detect_borders(pmesh, borders);
  assert(borders.size() == 2);

  // Triangulating holes.
  std::vector<Face_handle> patch_faces;
  BOOST_FOREACH(Halfedge_handle h, borders) {
    patch_faces.clear();
    CGAL::Polygon_mesh_processing::triangulate_hole_with_cdt_2(
      pmesh,
      h,
      std::back_inserter(patch_faces),
      CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
        get(CGAL::vertex_point, pmesh)).
        geom_traits(GT()));

    assert(patch_faces.size() == 25);
    std::cout << "* number of faces in the constructed patch: " <<
      patch_faces.size() << std::endl;
  }
  assert(pmesh.is_valid() && is_closed(pmesh));

  // Writing the file.
  if (save) {
    path = "";
    if (argc > 1) path = std::string(argv[1]);
    path += "w-4464.off";
    std::ofstream out(path.c_str(), std::ios_base::out);
    CGAL::set_ascii_mode(out);
    CGAL::write_off(out, pmesh);
    out.close();
    std::cout << "* finished writing the file" << std::endl;
  }
}

int main(int argc, char **argv) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel EI;
  typedef CGAL::Exact_predicates_exact_constructions_kernel EE;

  typedef CGAL::Surface_mesh<EI::Point_3> Surface_mesh_EI;
  typedef CGAL::Polyhedron_3<EE> Polyhedron_3_EE;

  test_triangulate_hole_with_cdt_2<Surface_mesh_EI, EI>(
    "exact_inexact", argc, argv, false);
  test_triangulate_hole_with_cdt_2<Polyhedron_3_EE, EE>(
    "exact_exact", argc, argv, true);
}
