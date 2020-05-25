#include <vector>
#include <fstream>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel GT;
typedef CGAL::Polyhedron_3<GT> Polyhedron;

typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;

int main(int argc, char **argv) {

  Polyhedron polyhedron;
  std::string path = "/Users/monet/Documents/gf/bugs/issue-4464/data/data.off";
  std::ifstream in(path.c_str(), std::ios_base::in);
  CGAL::set_ascii_mode(in);
  CGAL::read_off(in, polyhedron);
  in.close();
  std::cout << "* finished reading the file" << std::endl;

  // Incrementally fill the holes.
  std::size_t num_borders = 0;
  for (Halfedge_iterator h = polyhedron.halfedges_begin();
  h != polyhedron.halfedges_end(); ++h) {
    if (h->is_border()) {

      ++num_borders;
      std::vector<Facet_handle> patch_faces;
      CGAL::Polygon_mesh_processing::triangulate_hole_with_cdt_2(
        polyhedron,
        h,
        std::back_inserter(patch_faces),
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
          get(CGAL::vertex_point, polyhedron)).
          geom_traits(GT()));

      assert(patch_faces.size() == 25);
      std::cout << "* number of faces in the constructed patch: " <<
        patch_faces.size() << std::endl;
    }
  }
  assert(num_borders == 2);

  path = "/Users/monet/Documents/fork/logs/result-4464.off";
  std::ofstream out(path.c_str(), std::ios_base::out);
  CGAL::set_ascii_mode(out);
  CGAL::write_off(out, polyhedron);
  out.close();
  std::cout << "* finished writing the file" << std::endl;
}
