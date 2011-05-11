#include <CGAL/Triangulation_data_structure.h>
#include <iostream>
#include <vector>

int main()
{
  const int sdim = 27; // dimension of TDS with compile-time dimension
  const int ddim = 15; // dimension of TDS with dynamic dimension

  typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag>
          Dynamic_tds;
  typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<sdim> >
          Static_tds;
  typedef Static_tds::Face             Face;
  typedef Static_tds::Facet            Facet;
  typedef Static_tds::Vertex_handle    Vertex_handle;
  typedef Static_tds::Full_cell_handle Full_cell_handle;

  Static_tds  S(158);  // the argument is not taken into account.
  Dynamic_tds D(ddim); // the argument is taken into account.

  assert( sdim == S.ambient_dimension() );
  assert( ddim == D.ambient_dimension() );

  assert( -2 == S.current_dimension() );
  assert( S.is_valid() );

  std::vector<Vertex_handle> V(10);
  V[0] = S.insert_increase_dimension();
  assert( -1 == S.current_dimension() );

  for( int i = 1; i <= 5; ++i )
      V[i] = S.insert_increase_dimension(V[rand() % i]);

  assert( 4 == S.current_dimension() );
  assert( 6 == S.number_of_vertices() );
  assert( 6 == S.number_of_full_cells() );

  Full_cell_handle c = V[5]->full_cell();

  V[6] = S.insert_in_full_cell(c);
  
  assert( 7 == S.number_of_vertices() );
  assert( 10 == S.number_of_full_cells() );

  c = V[3]->full_cell();
  // ft will designate the Facet opposite to vertex 2 in c:
  Facet ft(c, 2);

  V[7] = S.insert_in_facet(ft);
  assert( 8 == S.number_of_vertices() );
  assert( 16 == S.number_of_full_cells() );

  c = V[3]->full_cell();
  // face will be the edge joining vertices 2 and 4 in full_cell c:
  Face face(c);
  face.set_index(0, 2);
  face.set_index(1, 4);

  V[8] = S.insert_in_face(face);
  assert( S.is_valid() );

  Full_cell_handle hole[2];
  hole[0] = V[8]->full_cell();
  hole[1] = hole[0]->neighbor(0);
  // ft will be a face on the boundary of hole[0] union hole[1]:
  ft = Facet(hole[0], 1);

  V[9] = S.insert_in_hole(hole, hole+2, ft);
  assert( S.is_valid() );

  std::cout << "Ok."<<std::endl;
  return 0;
}
