#include <CGAL/Triangulation_data_structure.h>

#include <vector>
#include <cassert>

int main()
{
  typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<7> > TDS;

  TDS S;
  assert( 7 == S.maximal_dimension() );
  assert( -2 == S.current_dimension() );
  assert( S.is_valid() );

  std::vector<TDS::Vertex_handle> V(10);
  V[0] = S.insert_increase_dimension(); //insert first vertex
  assert( -1 == S.current_dimension() );

  for( int i = 1; i <= 5; ++i )
      V[i] = S.insert_increase_dimension(V[0]);
  // the first 6 vertices have created a triangulation
  // of the 4-dimensional topological sphere
  // (the boundary of a five dimensional simplex).
  assert( 4 == S.current_dimension() );
  assert( 6 == S.number_of_vertices() );
  assert( 6 == S.number_of_full_cells() );

  TDS::Full_cell_handle c = V[5]->full_cell();
  V[6] = S.insert_in_full_cell(c);
  // full cell c is split in 5
  assert( 7 == S.number_of_vertices() );
  assert( 10 == S.number_of_full_cells() );

  c = V[3]->full_cell();
  TDS::Facet ft(c, 2); // the Facet opposite to vertex 2 in c
  V[7] = S.insert_in_facet(ft);
  // facet ft is split in 4 and the two incident cells are split accordingly
  assert( 8 == S.number_of_vertices() );
  assert( 16 == S.number_of_full_cells() );

  c = V[3]->full_cell();
  TDS::Face face(c);
  // meant to contain the edge joining vertices 2 and 4 of full_cell c
  face.set_index(0, 2); // namely vertex 2
  face.set_index(1, 4); // and vertex 4
  V[8] = S.insert_in_face(face);
  // face is split in 2, and all incident full cells also
  assert( S.is_valid() );

  TDS::Full_cell_handle hole[2];
  hole[0] = V[8]->full_cell();
  hole[1] = hole[0]->neighbor(0);
  // the hole is made of two adjacent full cells
  ft = TDS::Facet(hole[0], 1);  // a face on the boundary of hole[0]
  V[9] = S.insert_in_hole(hole, hole+2, ft);
  // the hole is triangulated by linking a new vertex to its boundary
  assert( S.is_valid() );
  return 0;
}
