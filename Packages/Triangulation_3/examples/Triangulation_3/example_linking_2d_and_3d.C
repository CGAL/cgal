// Triangulation_3/example_linking_2d_and_3d.C

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <cassert>

template < typename T3, typename Vb = CGAL::Triangulation_ds_vertex_base_2<> >
class My_vertex_2;

template < typename T2, typename Vb = CGAL::Triangulation_ds_vertex_base_3<> >
class My_vertex_3;

typedef CGAL::Triangulation_data_structure_2<My_vertex_2<CGAL::Dummy_tds_3> >
        TDS_2;
typedef CGAL::Triangulation_data_structure_3<My_vertex_3<TDS_2> >
        TDS_3;


template < typename T3, typename Vb >
class My_vertex_2
  : public Vb
{
public:
  typedef typename Vb::Face_handle    Face_handle;

  template <typename TDS2>
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_2<TDS_3, Vb2>                        Other;
  };

  My_vertex_2() {}

  My_vertex_2(Face_handle f) : Vb(f) {}

  typename T3::Vertex_handle v3;
};

template < typename T2, typename Vb >
class My_vertex_3
  : public Vb
{
public:
  typedef typename Vb::Cell_handle    Cell_handle;

  template <typename TDS2>
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_3<T2, Vb2>                           Other;
  };

  My_vertex_3() {}

  My_vertex_3(Cell_handle c) : Vb(c) {}

  typename T2::Vertex_handle v2;
};


int main() {
  TDS_2 t2;
  TDS_3 t3;

  TDS_2::Vertex_handle v2 = t2.insert_dim_up();
  TDS_3::Vertex_handle v3 = t3.insert_increase_dimension();

  v2->v3 = v3;
  v3->v2 = v2;

  assert(t2.is_valid());
  assert(t3.is_valid());
  return 0;
}
