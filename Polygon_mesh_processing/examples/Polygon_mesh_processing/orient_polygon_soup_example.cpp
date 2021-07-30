#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <vector>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>      Polyhedron;

//Optional visitor for orientation to demonstrate usage
struct Visitor : public CGAL::Polygon_mesh_processing::Default_orientation_visitor
{
  void non_manifold_edge(const std::size_t& id1, const std::size_t& id2) final
  {
    std::cout<<"The edge "<<id1<<", "<<id2<<" is not manifold."<<std::endl;
  }
  void non_manifold_vertex(const std::size_t & id) final
  {
     std::cout<<"The vertex "<<id<<" is not manifold."<<std::endl;
  }
  void duplicated_vertex(const std::size_t& v1, const std::size_t& v2) final
  {
    std::cout<<v1<<" has been duplicated, the new id is "<<v2<<"."<<std::endl;
  }
  void point_id_in_polygon_updated(const std::size_t& p_id, const std::size_t& i1, const std::size_t& i2) final
  {
    std::cout<<"In the polygon "<<p_id<<"< the index "<<i1<<" has been replaced by "<<i2<<"."<<std::endl;
  }
  void polygon_orientation_reversed(const std::size_t& p_id) final
  {
    std::cout<<"The polygon "<<p_id<<" has been reversed."<<std::endl;
  }
  };

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/tet-shuffled.off";

  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;

  if(!CGAL::IO::read_polygon_soup(filename, points, polygons) || points.empty())
  {
    std::cerr << "Cannot open file " << std::endl;
    return EXIT_FAILURE;
  }

  Visitor visitor;
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons, CGAL::parameters::visitor(visitor));

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
