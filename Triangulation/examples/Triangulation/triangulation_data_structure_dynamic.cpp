#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/assertions.h>

#include <vector>

int main()
{
  const int ddim = 5; // dimension of TDS with dynamic dimension
  typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag>
          TDS;
  typedef TDS::Vertex_handle    Vertex_handle;
  TDS D(ddim); // the argument is taken into account.

  CGAL_assertion( ddim == D.maximal_dimension() );
  CGAL_assertion( -2 == D.current_dimension() );
  CGAL_assertion( D.is_valid() );
  std::vector<Vertex_handle> V(5);
  V[0] = D.insert_increase_dimension();
  V[1] = D.insert_increase_dimension(V[0]);
  V[2] = D.insert_increase_dimension(V[0]);
  V[3] = D.insert_increase_dimension(V[0]);
  V[4] = D.insert_in_full_cell(V[3]->full_cell());
  CGAL_assertion( 6 == D.number_of_full_cells() );
  CGAL_assertion( 2 == D.current_dimension() );
  CGAL_assertion( D.is_valid() );
  return 0;
}
