
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/boost/graph/graph_traits_Arrangement_2.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <map>
#include <boost/unordered_map.hpp>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Arrangement_traits_2;
typedef CGAL::Arrangement_2<Arrangement_traits_2> Arrangement_2;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef CGAL::Triangulation_2<Kernel> Triangulation_2;
typedef CGAL::Linear_cell_complex<3, 3> Linear_cell_complex_3;


template <typename P>
void
fct(const P& )
{
  typedef typename boost::graph_traits<P>::vertex_descriptor vertex_descriptor;

  std::map<vertex_descriptor,int> M;
  vertex_descriptor vd;
  M.find(vd);

  boost::unordered_map<vertex_descriptor, int> U;
  U[vd] = 12;
}

void fct2()
{
  typedef Linear_cell_complex_3::Dart_handle dh;
  typedef Linear_cell_complex_3::Vertex_attribute_handle vh;

  { // For dart handle
  std::map<dh, int> M;
  dh e;
  M.find(e);
  boost::unordered_map<dh, int> U;
  U[e] = 12;
  }

  { // For vertex attribute handle
  std::map<vh, int> M;
  vh e;
  M.find(e);
  boost::unordered_map<vh, int> U;
  U[e] = 12;
  }
}

int main()
{
  Arrangement_2 A;
  fct(A);

  Polyhedron P;
  fct(P);

  Surface_mesh S;
  fct(S);

  Triangulation_2 T;
  fct(T);

  fct2();
  
  return 0;
}


