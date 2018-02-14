#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_plus_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_hierarchy_2.h>

#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Triangulation_2<Kernel> Triangulation;
typedef CGAL::Delaunay_triangulation_2<Kernel> DT2;
typedef CGAL::Regular_triangulation_2<Kernel> RT2;
typedef CGAL::Constrained_triangulation_2<Kernel> CT2;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel> CDT2;
typedef CGAL::Constrained_triangulation_plus_2<CDT2> CDTP2;
typedef CGAL::Triangulation_hierarchy_2<DT2> THDT2;
typedef CGAL::Triangulation_hierarchy_2<CDTP2> THCDTP2;

template<typename T>
void concept_check_triangulation() {
  boost::function_requires< boost::VertexListGraphConcept<T> >();
  boost::function_requires< boost::BidirectionalGraphConcept<T> >();
}

int main()
{
  concept_check_triangulation<Triangulation>();
  concept_check_triangulation<DT2>();
  concept_check_triangulation<RT2>();
  concept_check_triangulation<CT2>();
  concept_check_triangulation<CDT2>();
  concept_check_triangulation<CDTP2>();
  concept_check_triangulation<THDT2>();
  concept_check_triangulation<THCDTP2>();
  return 0;
}
