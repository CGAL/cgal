#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

template <typename K>
void test_orient(const char* file_name)
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  
  typedef typename K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Surface_mesh;
  
  //run test for a Polyhedron
  std::ifstream input(file_name);
  Polyhedron poly; // file should contain oriented polyhedron
  if ( !input || !(input >> poly) || poly.empty() )
  {
    std::cerr << "Error: cannot read Polyhedron : " << file_name << "\n";
    assert(false);
  }

  assert(CGAL::Polygon_mesh_processing::is_outward_oriented(poly));

  CGAL::Polygon_mesh_processing::reverse_face_orientations(poly);

  assert(!CGAL::Polygon_mesh_processing::is_outward_oriented(poly));

  std::cerr << "Polyhedron "<< file_name << " passed the test." << std::endl;
  poly.clear(); //empty memory
  input.close();

  //test again for a Surface_mesh
  Surface_mesh mesh;
  input.open(file_name);
  if (!(input >> mesh))
  {
    std::cerr << "Error: cannot read Surface_mesh : " << file_name << "\n";
    assert(false);
  }

  assert(CGAL::Polygon_mesh_processing::is_outward_oriented(mesh));

  CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);

  assert(!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh));

  std::cerr << "Surface_mesh " << file_name << " passed the test." << std::endl;
}

int main()
{

  test_orient<Epic>("data/elephant.off");
  test_orient<Epec>("data/elephant.off");

  std::cerr << "All done." << std::endl;
  return EXIT_SUCCESS;
}
