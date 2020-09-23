#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/IO/OFF_reader.h>

#include <vector>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>      Polyhedron;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/tet-shuffled.off";
  std::ifstream input(filename);

  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;

  if(!input || !CGAL::read_OFF(input, points, polygons) || points.empty())
  {
    std::cerr << "Cannot open file " << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

  Polyhedron mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  // Number the faces because 'orient_to_bound_a_volume' needs a face <--> index map
  int index = 0;
  for(Polyhedron::Face_iterator fb=mesh.facets_begin(), fe=mesh.facets_end(); fb!=fe; ++fb)
    fb->id() = index++;

  if(CGAL::is_closed(mesh))
    CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(mesh);

  std::ofstream out("tet-oriented1.off");
  out.precision(17);
  out << mesh;
  out.close();

  CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  std::ofstream out2("tet-oriented2.off");
  out2.precision(17);
  out2 << mesh;
  out2.close();

  return EXIT_SUCCESS;
}
