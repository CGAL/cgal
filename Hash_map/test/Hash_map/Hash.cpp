
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <map>
#include <boost/unordered_map.hpp>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef CGAL::Triangulation_2<Kernel> Triangulation_2;



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


int main()
{
  Polyhedron P;
  fct(P);

  Surface_mesh S;
  fct(S);

  Triangulation_2 T;
  fct(T);

  return 0;
}


