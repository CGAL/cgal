
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/boost/graph/graph_traits_Arrangement_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <map>
#include <boost/unordered_map.hpp>
#include <CGAL/boost/graph/helpers.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Arrangement_traits_2;
typedef CGAL::Arrangement_2<Arrangement_traits_2> Arrangement_2;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef CGAL::Triangulation_2<Kernel> Triangulation_2;
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_vertex_base_3<Kernel>,
    CGAL::Delaunay_triangulation_cell_base_3<Kernel>,
    CGAL::Parallel_tag>                          Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Triangulation_3;
#endif

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> Linear_cell_complex_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<3, 3> Linear_cell_complex_3_gmap;


template <typename P,typename Descriptor>
void
fct(const P& )
{
  std::map<Descriptor,int> M;
  Descriptor d;
  typename std::map<Descriptor,int>::const_iterator it = M.find(d);
  CGAL_USE(it);

  boost::unordered_map<Descriptor, int> U;
  U[d] = 12;
}

template<typename LCC>
void fct2()
{
  typedef typename LCC::Dart_handle dh;
  typedef typename LCC::Vertex_attribute_handle vh;

  { // For dart handle
  std::map<dh, int> M;
  dh e;
  typename std::map<dh, int>::const_iterator it = M.find(e);
  CGAL_USE(it);

  boost::unordered_map<dh, int> U;
  U[e] = 12;
  }

  { // For vertex attribute handle
  std::map<vh, int> M;
  vh e;
  typename std::map<vh, int>::const_iterator it = M.find(e);
  CGAL_USE(it);

  boost::unordered_map<vh, int> U;
  U[e] = 12;
  }
}

template <typename P>
void
fct3(const P& )
{
  typedef typename P::Vertex_handle vertex_descriptor;

  std::map<vertex_descriptor,int> M;
  vertex_descriptor vd;
  M.find(vd);

  boost::unordered_map<vertex_descriptor, int> U;
  U[vd] = 12;
}

template <class P>
void test_edge_hash_and_null(const P& p)
{
  typedef boost::graph_traits<P> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  for(halfedge_descriptor h : halfedges(p))
  {
    assert(
      hash_value( edge(h,p) ) ==
      hash_value( edge(opposite(h,p),p) ) );
  }
  face_descriptor f1=GT::null_face(), f2=GT::null_face();
  vertex_descriptor v1=GT::null_vertex(), v2=GT::null_vertex();
  halfedge_descriptor h1=GT::null_halfedge(), h2=GT::null_halfedge();

  assert(hash_value(f1)==hash_value(f2));
  assert(hash_value(v1)==hash_value(v2));
  assert(hash_value(h1)==hash_value(h2));
}

template <typename P>
void
fct4(const P& p)
{
  fct<P, typename boost::graph_traits<P>::vertex_descriptor>(p);
  fct<P, typename boost::graph_traits<P>::halfedge_descriptor>(p);
  fct<P, typename boost::graph_traits<P>::edge_descriptor>(p);
  fct<P, typename boost::graph_traits<P>::face_descriptor>(p);
  test_edge_hash_and_null(p);
}

int main()
{
  Arrangement_2 A;
  fct<Arrangement_2, boost::graph_traits<Arrangement_2>::vertex_descriptor>(A);
  fct<Arrangement_2, boost::graph_traits<Arrangement_2>::edge_descriptor>(A);

  Kernel::Point_3 p3;
  Polyhedron P;
  CGAL::make_triangle(p3,p3,p3,P);
  fct4(P);

  Surface_mesh S;
  CGAL::make_triangle(p3,p3,p3,S);
  fct4(S);

  Triangulation_2 T;
  T.insert(Kernel::Point_2(0,0));
  T.insert(Kernel::Point_2(0,1));
  T.insert(Kernel::Point_2(1,0));
  T.insert(Kernel::Point_2(1,1));
  fct4(T);

  fct2<Linear_cell_complex_3_cmap>();
  fct2<Linear_cell_complex_3_gmap>();

#ifdef CGAL_LINKED_WITH_TBB
  Triangulation_3 T3;
  fct3(T3);
#endif

  return 0;
}


