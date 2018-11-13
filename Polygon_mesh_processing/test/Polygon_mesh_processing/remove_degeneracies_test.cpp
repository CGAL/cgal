#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair.h>

#include <iostream>
#include <fstream>

//note : when
//CGAL::get_default_random()::get_seed() = 1473902576
//the last test (on trihole.off) does not terminate
//

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;

typedef CGAL::Surface_mesh<K::Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor      face_descriptor;

void detect_degeneracies(const Surface_mesh& mesh)
{
  std::vector<face_descriptor> dfaces;

  PMP::degenerate_faces(mesh, std::back_inserter(dfaces));
  PMP::degenerate_faces(faces(mesh), mesh, std::back_inserter(dfaces));
  PMP::degenerate_faces(mesh, std::back_inserter(dfaces), CGAL::parameters::all_default());
  PMP::degenerate_faces(faces(mesh), mesh, std::back_inserter(dfaces), CGAL::parameters::all_default());
  assert(!dfaces.empty());

  std::set<edge_descriptor> dedges;
  PMP::degenerate_edges(mesh, std::inserter(dedges, dedges.end()));
  PMP::degenerate_edges(edges(mesh), mesh, std::inserter(dedges, dedges.begin()));
  PMP::degenerate_edges(mesh, std::inserter(dedges, dedges.end()), CGAL::parameters::all_default());
  PMP::degenerate_edges(edges(mesh), mesh, std::inserter(dedges, dedges.begin()), CGAL::parameters::all_default());
  assert(dedges.empty());
}

void fix_degeneracies(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  detect_degeneracies(mesh);

  CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);
  assert( CGAL::is_valid_polygon_mesh(mesh) );
}

int main()
{
  fix_degeneracies("data_degeneracies/degtri_2dt_1edge_split_twice.off");
  fix_degeneracies("data_degeneracies/degtri_four-2.off");
  fix_degeneracies("data_degeneracies/degtri_four.off");
  fix_degeneracies("data_degeneracies/degtri_on_border.off");
  fix_degeneracies("data_degeneracies/degtri_three.off");
  fix_degeneracies("data_degeneracies/degtri_single.off");
  fix_degeneracies("data_degeneracies/trihole.off");

  return 0;
}
