#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <CGAL/draw_linear_cell_complex.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> LCC_3;
using Point=LCC_3::Point;

//==============================================================================
int main()
{
  LCC_3 lcc;
  CGAL::Linear_cell_complex_incremental_builder_3<LCC_3> ib(lcc);
  
  ib.add_vertex(Point(0,0,0)); // vertex 0
  ib.add_vertex(Point(1,0,0)); // vertex 1
  ib.add_vertex(Point(1,1,0)); // vertex 2
  ib.add_vertex(Point(0,1,0)); // vertex 3

  ib.add_vertex(Point(0,1,1)); // vertex 4
  ib.add_vertex(Point(0,0,1)); // vertex 5
  ib.add_vertex(Point(1,0,1)); // vertex 6
  ib.add_vertex(Point(1,1,1)); // vertex 7

  // Create a cube
  ib.begin_surface();
  ib.add_facet({0,1,2,3}); // Create a new facet v1: given all of its indices
  ib.add_facet({1,0,5,6});
  ib.add_facet({2,1,6,7});
  ib.add_facet({3,2,7,4});

  ib.begin_facet(); // Create a new facet v2: begin facet
  ib.add_vertex_to_facet(0); // all incrementally its indices
  ib.add_vertex_to_facet(3);
  ib.add_vertex_to_facet(4);
  ib.add_vertex_to_facet(5);
  ib.end_facet(); // end facet
  
  ib.add_facet({5,4,7,6});
  
  ib.end_surface();

  ib.add_vertex(Point(-1, 0.5, 0.5)); // vertex 8
  
  // Create a pyramid, sharing one of its face with the cube
  ib.begin_surface();
  ib.add_facet({3,0,5,4});
  ib.add_facet({0,3,8});
  ib.add_facet({3,4,8});
  ib.add_facet({4,5,8});
  ib.add_facet({5,0,8});
  ib.end_surface();  

  // Draw the lcc and display its characteristics
  lcc.display_characteristics(std::cout)<<std::endl;
  CGAL::draw(lcc);
  
  return EXIT_SUCCESS;
}
