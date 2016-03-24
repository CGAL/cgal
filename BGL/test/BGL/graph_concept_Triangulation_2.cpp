#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Triangulation_2<Kernel> Triangulation;

template<typename T>
void concept_check_triangulation() {
  boost::function_requires< boost::VertexListGraphConcept<T> >();
  boost::function_requires< boost::BidirectionalGraphConcept<T> >();
}

int main()
{
  concept_check_triangulation<Triangulation>();
  return 0;
}
