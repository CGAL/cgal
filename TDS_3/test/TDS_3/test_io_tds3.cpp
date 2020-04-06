#include <fstream>
#include <cassert>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Triangulation_data_structure_3<>      Tds;
typedef Tds::size_type                              size_type;
typedef Tds::Cell_handle                            Cell_handle;
typedef Tds::Vertex_handle                          Vertex_handle;

template <typename T1, typename T2>
struct Update_vertex
{
  typedef typename T1::Vertex                  V1;
  typedef typename T2::Vertex                  V2;

  V2 operator()(const V1&)
  {
    return V2();
  }

  void operator()(const V1&, V2&)
  {

  }


}; // end struct Update_vertex

template <typename T1, typename T2>
struct Update_cell
{
  typedef typename T1::Cell                  C1;
  typedef typename T2::Cell                  C2;

  C2 operator()(const C1&)
  {
    return C2();
  }

  void operator()(const C1&, C2&)
  {

  }


}; // end struct Update_vertex

int main()
{
  Tds T;
  std::vector<Vertex_handle> PV(7);
  PV[0] = T.insert_increase_dimension();
  // each of the following insertions of vertices increases the dimension
  for ( int i=1; i<5; i++ ) {
    PV[i] = T.insert_increase_dimension(PV[0]);
  }
  // we now have a simplex in dimension 4
  // cell incident to PV[0]
  Cell_handle c = PV[0]->cell();
  int ind=0;
  // PV[0] is the vertex of index ind in c
  // insertion of a new vertex in the facet opposite to PV[0]
  PV[5] = T.insert_in_facet(c, ind);
  // insertion of a new vertex in c
  PV[6] = T.insert_in_cell(c);
  std::ofstream out("tr");
  out << T;
  out.close();
  Tds T2;
  std::ifstream in("tr");
  T2.file_input<Tds,Update_vertex<Tds, Tds>, Update_cell<Tds, Tds> >(in);
  in.close();
  return 0;
}
