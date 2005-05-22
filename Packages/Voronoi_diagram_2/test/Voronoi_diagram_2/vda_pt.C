#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include "vda_test.h"

#include <CGAL/Simple_cartesian.h>
#include "Delaunay_graph_concept.h"
#include "Voronoi_traits_concept.h"


typedef CGAL::Simple_cartesian<double>              K;
typedef CGAL::Delaunay_graph_concept<K>             DG;
typedef CGAL::Voronoi_traits_concept<DG>            VT;
typedef CGAL::Voronoi_diagram_adaptor_2<DG,VT>      VDA;


int main()
{
  typedef Project_site<DG::Vertex_handle,DG::Site_2>   Project_site;
  typedef Project_dual<VDA,DG::Site_2>                 Project_dual;
  typedef VDA_Tester<VDA,Project_site,Project_dual>    Tester;

  Project_site   project_site;
  Project_dual   project_dual;

  Tester test(project_site, project_dual);

  test.reset_timers();

  test("data/empty.cin");
  test("data/data1.pt.cin");

  test.print_times();

  return 0;
}
