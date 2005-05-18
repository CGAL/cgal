#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include "vda_test.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>
#include <CGAL/Apollonius_graph_Voronoi_traits_2.h>


typedef CGAL::Simple_cartesian<double>                        Rep;
typedef CGAL::Apollonius_graph_filtered_traits_2<Rep>         K;

#if 1 // definitions for hierarchy

#include <CGAL/Apollonius_graph_hierarchy_2.h>

typedef CGAL::Apollonius_graph_hierarchy_2<K>                 AG;
#else
typedef CGAL::Apollonius_graph_2<K>                           AG;
#endif
//typedef CGAL::Apollonius_graph_Voronoi_traits_2<AG>           VT;
typedef CGAL::Apollonius_graph_cached_Voronoi_traits_2<AG>    VT;
typedef CGAL::Voronoi_diagram_adaptor_2<AG,VT>                VDA;


int main()
{
  typedef Project_site<AG::Vertex_handle,AG::Site_2>    Project_site;
  typedef Project_ag_dual<VDA,AG::Site_2>               Project_ag_dual;
  typedef VDA_Tester<VDA,Project_site,Project_ag_dual>  Tester;

  Project_site      project_site;
  Project_ag_dual   project_ag_dual;

  Tester test(project_site, project_ag_dual);

  test("data/empty.cin");
  test("data/data1.ag.cin");
  test("data/data2.ag.cin");
  test("data/data3.ag.cin");
  test("data/data4.ag.cin");
  test("data/degenerate.ag.cin");

  return 0;
}
