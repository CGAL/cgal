// #define CGAL_PMP_SMOOTHING_DEBUG
#define CGAL_PMP_COMPUTE_NORMAL_DEBUG
#define CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
#define CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <unordered_map>

namespace PMP = ::CGAL::Polygon_mesh_processing;
namespace CP = ::CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;

typedef CGAL::Surface_mesh<K::Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor      face_descriptor;

void fix_self_intersections(const char* mesh_filename,
                            const char* mesh_selection_filename = nullptr)
{
  std::cout << std::endl << "---------------" << std::endl;
  std::cout << "Test " << mesh_filename << std::endl;
  if(mesh_selection_filename)
    std::cout << "With selection " << mesh_selection_filename << std::endl;

  std::ifstream input(mesh_filename);
  Surface_mesh mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Error: " << mesh_filename << " is not a valid off file.\n";
    std::exit(1);
  }

  std::list<face_descriptor> selected_faces;

  if(mesh_selection_filename)
  {
    std::ifstream selection_input(mesh_selection_filename);
    std::string line;
    // skip the first line (faces are on the second line)
    if(!selection_input || !std::getline(selection_input, line) || !std::getline(selection_input, line))
    {
      std::cerr << "Error: could not read selection: " << mesh_selection_filename << std::endl;
      std::exit(1);
    }

    std::istringstream face_line(line);
    std::size_t face_id;
    while(face_line >> face_id)
      selected_faces.push_back(*(faces(mesh).begin() + face_id));
    std::cout << selected_faces.size() << " faces selected" << std::endl;

    PMP::experimental::remove_self_intersections(selected_faces, mesh);
  }
  else
  {
    PMP::experimental::remove_self_intersections(mesh);
  }

  std::ofstream out("post_repair.off");
  out.precision(17);
  out << mesh;
  out.close();

  assert(CGAL::is_valid_polygon_mesh(mesh));
  assert(!PMP::does_self_intersect(mesh));
}

void fix_local_self_intersections(const char* fname)
{
  std::cout << std::endl << "-----LOCAL------" << std::endl;
  std::cout << "Test " << fname << std::endl;

  std::ifstream input(fname);
  Surface_mesh mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Error: " << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  PMP::experimental::remove_self_intersections(mesh, CP::apply_per_connected_component(true));

  std::ofstream out("post_local_repair.off");
  out.precision(17);
  out << mesh;
  out.close();

  assert(CGAL::is_valid_polygon_mesh(mesh));
  assert(!PMP::does_self_intersect(mesh));
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

#if 0
  fix_local_self_intersections("data_repair/brain.off");
  fix_self_intersections("data_repair/brain.off")

  fix_self_intersections("data_repair/flute.off");
  fix_self_intersections("data_repair/dinosaur.off");
  fix_self_intersections("data_repair/hemispheres.off");

  // selection is adjacent to a self-intersection but does not contain any intersection
  fix_self_intersections("data_repair/brain.off", "data_repair/brain-complete.selection.txt");

  // selection covers nicely a self-intersection
  fix_self_intersections("data_repair/brain.off", "data_repair/brain-adjacent.selection.txt");

  // selection contains part of a self intersection
  fix_self_intersections("data_repair/brain.off", "data_repair/brain-partial.selection.txt");

  // selection contains disjoint parts of a self intersection
  fix_self_intersections("data_repair/brain.off", "data_repair/brain-disjoint.selection.txt");

  // Remove only self-intersections within a single hemisphere
  fix_self_intersections("data_repair/hemispheres.off");
  fix_self_intersections("data_repair/hemispheres.off", "data_repair/hemispheres-half.selection.txt");
#endif

  return 0;
}
