
#include <boost/concept/assert.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/circulator.h>
#include <CGAL/Circulator/Circulator_concepts.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;
typedef SM::Halfedge_around_target_circulator Halfedge_around_source_circulator;
typedef SM::Halfedge_around_target_circulator Halfedge_around_target_circulator;
typedef SM::Halfedge_around_target_circulator Vertex_around_target_circulator;
typedef SM::Halfedge_around_target_circulator Face_around_target_circulator;
typedef SM::Halfedge_around_target_circulator Halfedge_around_face_circulator;
typedef SM::Halfedge_around_target_circulator Vertex_around_face_circulator;
typedef SM::Halfedge_around_target_circulator Face_around_face_circulator;


template <typename Circ>
void test()
{
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<Circ>));
  Circ circ;
  if(circ){}
  if(circ == nullptr){}
}

int main()
{
  test<Halfedge_around_source_circulator>();
  test<Halfedge_around_target_circulator>();
  test<Vertex_around_target_circulator>();
  test<Face_around_target_circulator>();

  test<Halfedge_around_face_circulator>();
  test<Vertex_around_face_circulator>();
  test<Face_around_face_circulator>();

  return 0;
}
