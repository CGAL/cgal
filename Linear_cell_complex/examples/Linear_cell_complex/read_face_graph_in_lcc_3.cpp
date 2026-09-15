#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3_to_lcc.h>
#include <CGAL/boost/graph/generators.h>


typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef CGAL::Polyhedron_3<LCC::Traits> Polyhedron;
typedef LCC::Point Point_3;
typedef CGAL::Surface_mesh<Point_3>  Surface_mesh;


int main()
{
  Polyhedron polyhedron;
  Surface_mesh surface_mesh;
  CGAL::make_tetrahedron(Point_3(0,0,0), Point_3(1,0,0),Point_3(0,1,0),Point_3(0,0,1), polyhedron);
  CGAL::make_tetrahedron(Point_3(0,0,0), Point_3(1,0,0),Point_3(0,1,0),Point_3(0,0,1), surface_mesh);


  LCC lccp, lccs;
  CGAL::polyhedron_3_to_lcc(lccp, polyhedron);
  CGAL::polyhedron_3_to_lcc(lccs, surface_mesh);

  std::cout << lccp << std::endl;
  std::cout << lccs << std::endl;
  return 0;
}