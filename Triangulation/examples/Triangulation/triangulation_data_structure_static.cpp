#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/assertions.h>

#include <vector>

int main()
{
  typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<7> > TDS;

  TDS S;
  CGAL_assertion( 7 == S.maximal_dimension() );
  CGAL_assertion( -2 == S.current_dimension() );
  CGAL_assertion( S.is_valid() );

  std::vector<TDS::Vertex_handle> V(10);
  V[0] = S.insert_increase_dimension(); //insert first vertex
  CGAL_assertion( -1 == S.current_dimension() );

  for( int i = 1; i <= 5; ++i )
      V[i] = S.insert_increase_dimension(V[0]);
  // the first 6 vertices have created a triangulation
  // of the 4-dimensional topological sphere 
  // (the boundary of a five dimensional simplex).
  CGAL_assertion( 4 == S.current_dimension() );
  CGAL_assertion( 6 == S.number_of_vertices() );
  CGAL_assertion( 6 == S.number_of_full_cells() );

  TDS::Full_cell_handle c = V[5]->full_cell();
  V[6] = S.insert_in_full_cell(c);
  // full cell c is split in 5
  CGAL_assertion( 7 == S.number_of_vertices() );
  CGAL_assertion( 10 == S.number_of_full_cells() );

  c = V[3]->full_cell();
  TDS::Facet ft(c, 2); // the Facet opposite to vertex 2 in c
  V[7] = S.insert_in_facet(ft);
  // facet ft is split in 4 and the two incident cells are split accordingly
  CGAL_assertion( 8 == S.number_of_vertices() );
  CGAL_assertion( 16 == S.number_of_full_cells() );

  c = V[3]->full_cell();
  TDS::Face face(c);  
  // meant to contain the edge joining vertices 2 and 4 of full_cell c
  face.set_index(0, 2); // namely vertex 2
  face.set_index(1, 4); // and vertex 4
  V[8] = S.insert_in_face(face);
  // face is split in 2, and all incident full cells also
  CGAL_assertion( S.is_valid() );

  TDS::Full_cell_handle hole[2];
  hole[0] = V[8]->full_cell();
  hole[1] = hole[0]->neighbor(0);
  // the hole is made of two adjacent full cells
  ft = TDS::Facet(hole[0], 1);  // a face on the boundary of hole[0]
  V[9] = S.insert_in_hole(hole, hole+2, ft);
  // the hole is triangulated by linking a new vertex to its boundary
  CGAL_assertion( S.is_valid() );
  return 0;
}
