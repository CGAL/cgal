#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/boost/graph/graph_traits_Arrangement_2.h>
#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Cartesian<CGAL::Exact_rational>         Kernel;

typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

template<typename T>
void concept_check_arrangement()
{
  boost::function_requires< boost::GraphConcept<T> >();
  boost::function_requires< boost::IncidenceGraphConcept<T> >();
  boost::function_requires< boost::VertexListGraphConcept<T> >();
  boost::function_requires< boost::EdgeListGraphConcept<T> >();
  boost::function_requires< boost::BidirectionalGraphConcept<T> >();

//  boost::function_requires< CGAL::HalfedgeGraphConcept<T> >();
//  boost::function_requires< CGAL::HalfedgeListGraphConcept<T> >();
//  boost::function_requires< CGAL::FaceGraphConcept<T> >();
//  boost::function_requires< CGAL::FaceListGraphConcept<T> >();
//  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<T> >();
//  boost::function_requires< CGAL::MutableFaceGraphConcept<T> >();

  // null
  boost::graph_traits<T>::null_vertex();
//  boost::graph_traits<T>::null_halfedge();
//  boost::graph_traits<T>::null_face();
}

int main()
{
  concept_check_arrangement<Arrangement_2>();

  return 0;
}
