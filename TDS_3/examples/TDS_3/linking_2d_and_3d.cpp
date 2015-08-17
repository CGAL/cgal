#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <cassert>

// declare the 2D vertex base type, parametrized by some 3D TDS.
template < typename T3, typename Vb = CGAL::Triangulation_ds_vertex_base_2<> >
class My_vertex_2;

// declare the 3D vertex base type, parametrized by some 2D TDS.
template < typename T2, typename Vb = CGAL::Triangulation_ds_vertex_base_3<> >
class My_vertex_3;

// Then, we have to break the dependency cycle.

// we need to refer to a dummy 3D TDS.
typedef CGAL::Triangulation_ds_vertex_base_3<>::Triangulation_data_structure
        Dummy_tds_3;
// the 2D TDS, initially plugging a dummy 3D TDS in the vertex type
// (to break the dependency cycle).
typedef CGAL::Triangulation_data_structure_2<My_vertex_2<Dummy_tds_3> >  TDS_2;
// the 3D TDS, here we can plug the 2D TDS directly.
typedef CGAL::Triangulation_data_structure_3<My_vertex_3<TDS_2> >        TDS_3;


template < typename T3, typename Vb >
class My_vertex_2
  : public Vb
{
public:
  typedef typename Vb::Face_handle    Face_handle;

  template <typename TDS2>
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    // we also have to break the cycle here by hardcoding TDS_3 instead of T3.
    typedef My_vertex_2<TDS_3, Vb2>                        Other;
  };

  My_vertex_2() {}

  My_vertex_2(Face_handle f) : Vb(f) {}

  // we store a vertex handle of the 3D TDS.
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

  // we store a vertex handle of the 2D TDS.
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
