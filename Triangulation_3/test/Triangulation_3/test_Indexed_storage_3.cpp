#define CGAL_TDS_USE_INDEXED_STORAGE_3
#include <CGAL/Triangulation_3.h>
#include <cassert>

bool del = false;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_triangulation_3.h>
// Explicit instantiation of the whole class :
// template class CGAL::Triangulation_3<K>;

int main()
{
  typedef CGAL::Triangulation_data_structure_3<CGAL::VertexWithPoint<K>, CGAL::Cell<>> Tds;
  typedef CGAL::Triangulation_3<K,Tds>                               Cls3;

  _test_cls_triangulation_3( Cls3() );

  // Test operator== between triangulations of different Tds types.

  // typedef CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<CGAL::Triangulation_vertex_base_with_info_3<int, K> > > Cls3_2;

  // assert(Cls3() == Cls3_2());
  std::cout << "done" << std::endl;
  return 0;
}
