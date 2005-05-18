#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include "vda_test.h"

#include <CGAL/MP_Float.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>
#include <CGAL/Segment_Voronoi_diagram_Voronoi_traits_2.h>

typedef CGAL::MP_Float  NT;
typedef CGAL::Ring_tag  MTag;

typedef CGAL::Simple_cartesian<NT>      K;
typedef CGAL::Simple_cartesian<double>  DK;
typedef
CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2<DK>
Gt;
//CGAL::Segment_Voronoi_diagram_traits_without_intersections_2<K,MTag>  Gt;

#if 1 // definitions for hierarchy

#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>

typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>                 SVD;
#else
typedef CGAL::Segment_Voronoi_diagram_2<Gt>                           SVD;
#endif
//typedef CGAL::Segment_Voronoi_diagram_Voronoi_traits_2<SVD>           VT;
typedef CGAL::Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD>    VT;
typedef CGAL::Voronoi_diagram_adaptor_2<SVD,VT>                       VDA;


int main()
{
  typedef Project_site<SVD::Vertex_handle,SVD::Site_2>   Project_site;
  typedef Project_primal<VDA,SVD::Point_2>               Project_primal;
  typedef VDA_Tester<VDA,Project_site,Project_primal>    Tester;

  Project_site    project_site;
  Project_primal  project_primal;

  Tester test(project_site, project_primal);

  test("data/empty.cin");
  test("data/complicated.svd.cin");
  test("data/non-degenerate.svd.cin");
  test("data/data0.svd.cin");
  test("data/data1.svd.cin");
  test("data/data2.svd.cin");
  test("data/data3.svd.cin");
  test("data/data4.svd.cin");
  test("data/data5.svd.cin");
  test("data/data6.svd.cin");
  test("data/data7.svd.cin");
  test("data/data8.svd.cin");
  test("data/data9.svd.cin");

  return 0;
}
