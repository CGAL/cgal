#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;

void test(const char* file_name)
{
  //run test for a Polyhedron
  std::ifstream input(file_name);
  Polyhedron poly; // file should contain oriented polyhedron
  if ( !input || !(input >> poly) || poly.empty() )
  {
    std::cerr << "Error: cannot read Polyhedron : " << file_name << "\n";
    CGAL_assertion(false);
  }

  bool before = CGAL::Polygon_mesh_processing::is_outward_oriented(poly);
  CGAL_assertion(before);

  poly.inside_out();

  bool after = CGAL::Polygon_mesh_processing::is_outward_oriented(poly);
  CGAL_assertion(!after);

  std::cerr << "Polyhedron "<< file_name << " passed the test." << std::endl;
  poly.clear(); //empty memory
  input.close();

  //test again for a Surface_mesh
  Surface_mesh mesh;
  input.open(file_name);
  if (!(input >> mesh))
  {
    std::cerr << "Error: cannot read Surface_mesh : " << file_name << "\n";
    CGAL_assertion(false);
  }
  before = CGAL::Polygon_mesh_processing::is_outward_oriented(mesh);
  CGAL_assertion(before);

  poly.inside_out();

  after = CGAL::Polygon_mesh_processing::is_outward_oriented(mesh);
  CGAL_assertion(!after);

  std::cerr << "Surface_mesh " << file_name << " passed the test." << std::endl;
}

int main()
{
  test("data/elephant.off");

  std::cerr << "All done." << std::endl;
}
