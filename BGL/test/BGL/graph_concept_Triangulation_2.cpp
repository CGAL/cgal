#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_plus_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_hierarchy_2.h>

#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double>                                Kernel;

typedef CGAL::Triangulation_2<Kernel>                                 Triangulation;
typedef CGAL::Delaunay_triangulation_2<Kernel>                        DT2;
typedef CGAL::Regular_triangulation_2<Kernel>                         RT2;
typedef CGAL::Constrained_triangulation_2<Kernel>                     CT2;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel>            CDT2;
typedef CGAL::Constrained_triangulation_plus_2<CDT2>                  CDTP2;

typedef CGAL::Triangulation_vertex_base_2<Kernel>                     Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                  TDS;
typedef CGAL::Exact_predicates_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;
typedef CGAL::Triangulation_hierarchy_2<CDT>                          THCDT2;
typedef CGAL::Constrained_triangulation_plus_2<THCDT2>                THCDTP2;

template<typename T>
void concept_check_triangulation()
{
  boost::function_requires< boost::GraphConcept<T> >();
  boost::function_requires< boost::IncidenceGraphConcept<T> >();
  boost::function_requires< boost::VertexListGraphConcept<T> >();
  boost::function_requires< boost::EdgeListGraphConcept<T> >();
  boost::function_requires< boost::BidirectionalGraphConcept<T> >();

  boost::function_requires< CGAL::HalfedgeGraphConcept<T> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<T> >();
  boost::function_requires< CGAL::FaceGraphConcept<T> >();
  boost::function_requires< CGAL::FaceListGraphConcept<T> >();

  // Triangulations are not mutable graphs
//  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<T> >();
//  boost::function_requires< CGAL::MutableFaceGraphConcept<T> >();

  // null
  boost::graph_traits<T>::null_vertex();
  boost::graph_traits<T>::null_halfedge();
  boost::graph_traits<T>::null_face();
}

int main()
{
  concept_check_triangulation<Triangulation>();
  concept_check_triangulation<DT2>();
  concept_check_triangulation<RT2>();
  concept_check_triangulation<CT2>();
  concept_check_triangulation<CDT2>();
  concept_check_triangulation<CDTP2>();
  concept_check_triangulation<THCDT2>();
  concept_check_triangulation<THCDTP2>();

  return 0;
}
