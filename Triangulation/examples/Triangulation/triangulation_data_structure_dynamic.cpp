#include <CGAL/Triangulation_data_structure.h>
#include <cassert>

#include <vector>

int main()
{
  const int ddim = 5; // dimension of TDS with dynamic dimension
  typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag>
          TDS;
  typedef TDS::Vertex_handle    Vertex_handle;
  TDS D(ddim); // the argument is taken into account.

  assert( ddim == D.maximal_dimension() );
  assert( -2 == D.current_dimension() );
  assert( D.is_valid() );
  std::vector<Vertex_handle> V(5);
  V[0] = D.insert_increase_dimension();
  V[1] = D.insert_increase_dimension(V[0]);
  V[2] = D.insert_increase_dimension(V[0]);
  V[3] = D.insert_increase_dimension(V[0]);
  V[4] = D.insert_in_full_cell(V[3]->full_cell());
  assert( 6 == D.number_of_full_cells() );
  assert( 2 == D.current_dimension() );
  assert( D.is_valid() );
  return 0;
}
