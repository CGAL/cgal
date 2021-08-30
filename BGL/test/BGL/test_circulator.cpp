#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/concept/assert.hpp>
#include <CGAL/Circulator/Circulator_concepts.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

typedef boost::graph_traits<Polyhedron> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::edge_descriptor edge_descriptor;
typedef GraphTraits::out_edge_iterator out_edge_iterator;
typedef GraphTraits::in_edge_iterator in_edge_iterator;

typedef CGAL::Halfedge_around_face_circulator<Polyhedron> halfedge_around_face_circulator;
typedef CGAL::Halfedge_around_target_circulator<Polyhedron> halfedge_around_target_circulator;
typedef CGAL::Vertex_around_target_circulator<Polyhedron> vertex_around_target_circulator;
typedef CGAL::Face_around_target_circulator<Polyhedron> face_around_target_circulator;

typedef CGAL::Halfedge_around_source_circulator<Polyhedron> halfedge_around_source_circulator;

typedef CGAL::Halfedge_around_target_iterator<Polyhedron> halfedge_around_target_iterator;
typedef CGAL::Halfedge_around_face_iterator<Polyhedron> halfedge_around_face_iterator;
typedef CGAL::Face_around_face_iterator<Polyhedron> face_around_face_iterator;
typedef CGAL::Vertex_around_target_iterator<Polyhedron> vertex_around_target_iterator;
int main(int argc, char* argv[])
{

  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_face_circulator>));
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_target_circulator>));
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<vertex_around_target_circulator>));
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<face_around_target_circulator>));
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_source_circulator>));

  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_source_circulator>));

   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<face_around_face_iterator>));
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<halfedge_around_face_iterator>));
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<halfedge_around_target_iterator>));
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<vertex_around_target_iterator>));

   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<in_edge_iterator>));
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<out_edge_iterator>));

  std::ifstream in((argc>1)?argv[1]:"data/cube.off");
  Polyhedron P;
  in >> P;

  halfedge_descriptor hd = *halfedges(P).first;
  {
    halfedge_around_face_circulator hafc(hd,P), done(hafc);

    do {
      std::cout << get(CGAL::vertex_point, P, target(*hafc,P)) << std::endl;
      ++hafc;
    }while(hafc != done);
  }

  {
    halfedge_around_target_circulator havc(hd,P), done(havc);
    vertex_descriptor vd = target(hd,P);
    do {
      halfedge_descriptor hd2 = *havc;
      assert(target(hd2,P) == vd);
      std::cout << get(CGAL::vertex_point, P, target(*havc,P)) << std::endl;
      ++havc;
    }while(havc != done);
  }
  {
    vertex_around_target_circulator havc(hd,P), done(havc);

    do {
      std::cout << get(CGAL::vertex_point, P, *havc) << std::endl;
      ++havc;
    }while(havc != done);
  }
  {
    face_around_target_circulator havc(hd,P), done(havc);

    do {
      //std::cout << get(CGAL::vertex_point, P, *havc) << std::endl;
      ++havc;
    }while(havc != done);
  }
  {
    halfedge_around_source_circulator havc(hd,P), done(havc);
    vertex_descriptor vd = source(hd,P);
    do {
      halfedge_descriptor hd2 = *havc;
      assert(source(hd2,P) == vd);
      std::cout << get(CGAL::vertex_point, P, target(*havc,P)) << std::endl;
      ++havc;
    }while(havc != done);
  }

  {
    halfedge_around_target_iterator vit, end;
    vertex_descriptor vd = target(hd,P);
    boost::tie(vit,end) = halfedges_around_target(hd,P);
    while(vit!= end) {
      halfedge_descriptor hd = *vit;
      assert(target(hd,P) == vd);
      std::cout << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
      ++vit;
    }
  }

   {
    halfedge_around_face_iterator vit, end;
    boost::tie(vit,end) = halfedges_around_face(hd,P);

    while(vit!= end) {
      halfedge_descriptor hd = *vit;
      std::cout << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
      ++vit;
    }
  }


  {
    out_edge_iterator ohi, end;
    for(boost::tie(ohi,end) = out_edges(target(hd,P),P); ohi != end; ++ohi){
      edge_descriptor ed = *ohi;
      halfedge_descriptor hd2 = halfedge(ed,P);
      std::cout << get(CGAL::vertex_point, P, target(hd2,P)) << std::endl;
    }
  }

  {
    for(edge_descriptor ed : out_edges(target(hd,P),P)){
      halfedge_descriptor hd2 = halfedge(ed,P);
      std::cout << get(CGAL::vertex_point, P, target(hd2,P)) << std::endl;
    }
  }
  return 0;
}
