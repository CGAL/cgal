#include <CGAL/boost/graph/graph_traits_Linear_cell_complex.h>
#include <CGAL/boost/graph/properties_Linear_cell_complex.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;

struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, Refs > Dart;
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attribute;
    typedef CGAL::Cell_attribute< Refs > Face_attribute;
    typedef CGAL::cpp11::tuple<Vertex_attribute, void, Face_attribute> Attributes;
  };
};

typedef CGAL::Linear_cell_complex<2, 3, MyTraits, Myitem> LCC;

typedef boost::graph_traits< LCC > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;
typedef LCC::Traits::Point Point_3;

template<typename Graph> 
void concept_check_polyhedron() {
  boost::function_requires< boost::GraphConcept<LCC> >();
  boost::function_requires< boost::VertexListGraphConcept<LCC> >();
  boost::function_requires< boost::EdgeListGraphConcept<LCC> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<LCC> >();
  boost::function_requires< boost::MutableGraphConcept<LCC> >();
  boost::function_requires< CGAL::FaceListGraphConcept<LCC> >();
  boost::function_requires< CGAL::HalfedgeGraphConcept<LCC> >();
  boost::function_requires< boost::IncidenceGraphConcept<LCC> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<LCC> >();
  boost::function_requires< boost::BidirectionalGraphConcept<LCC> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<LCC> >();

  // TODO Perhaps to remove ?
  boost::function_requires< CGAL::FaceGraphConcept<LCC> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<LCC> >();

  boost::function_requires< boost::concepts::PropertyGraph<
    LCC, halfedge_descriptor, boost::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    LCC, halfedge_descriptor, boost::halfedge_external_index_t> >();
  boost::function_requires< boost::PropertyGraphConcept<
    LCC, vertex_descriptor, boost::vertex_point_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    LCC, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    LCC, edge_descriptor, boost::edge_external_index_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    LCC, face_descriptor, boost::face_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    LCC, face_descriptor, boost::face_external_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
     LCC, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    LCC, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    LCC, vertex_descriptor, boost::vertex_external_index_t> >();
    
  // null
  boost::graph_traits<LCC>::null_vertex();
  boost::graph_traits<LCC>::null_face();
}

template<typename LCC>
void runtime_check_halfedgegraph()
{
  // u        v
  // +--------+
  // |\      /|
  // | \ f2 / |
  // |  \y /  |
  // | f3\/ f1|
  // |   /\   |
  // |  /  \  |
  // | / f4 \ |
  // |/      \|
  // +--------+
  // w        x
  LCC p;
  vertex_descriptor 
    u = add_vertex(Point_3(0,2,0), p),
    v = add_vertex(Point_3(2,2,0), p),
    w = add_vertex(Point_3(0,0,0), p),
    x = add_vertex(Point_3(2,0,0), p),
    y = add_vertex(Point_3(1,1,0), p);
  
  add_edge(v, u, p);
  add_edge(u, y, p);
  add_edge(y, v, p);
  add_face(p);

  add_edge(u, w, p);
  add_edge(w, y, p);
  add_face(p);

  add_edge(w, x, p);
  add_edge(x, y, p);
  add_face(p);

  add_edge(x, v, p);
  add_face(p);


  std::cout<<"num_edges(p)"<<num_edges(p)<<std::endl;
  std::cout<<"num_halfedges(p)"<<num_halfedges(p)<<std::endl;
  std::cout<<"num_faces(p)"<<num_faces(p)<<std::endl;
  

  assert(num_edges(p) ==  8);
  assert(num_halfedges(p) == 16);
  // WRONG ASSERT assert(num_faces(p) == 4);
}

int
main()
{
  LCC lcc;
  typedef typename boost::graph_traits<LCC>::halfedge_iterator iter;
  iter it;
    
    
    /*    CGAL::Halfedge_around_target_iterator<LCC> it;

    vd myvd;
    
    typedef typename boost::graph_traits<LCC>::in_edge_iterator iter;
    iter myiter;
  */ 
    // std::pair<iter, iter> p; // = in_edges(myvd, lcc);*/
    
  /*
  const LCC& lcc2 = lcc;
  
  iter_type a1(lcc);
  const iter_type a2(a1);
  iter_type a3(a2);
  //   a1=a2;


  CGAL::CMap_one_dart_per_cell_iterator<LCC, 1> toto1(lcc);
  const CGAL::CMap_one_dart_per_cell_iterator<LCC, 1> toto2(lcc);

  toto1=toto2;
  

  
  //iter_type it(lcc);
  
  std::pair<typename boost::graph_traits<LCC>::edge_iterator, typename boost::graph_traits<LCC>::edge_iterator> it=std::make_pair<iter_type, iter_type>
    (iter_type(lcc),
     iter_type(lcc)), it2=std::make_pair<iter_type, iter_type>
    (iter_type(lcc),
     iter_type(lcc));

  it=it2;
  
  iter_type t1(lcc);
  
  std::pair<iter_type, iter_type> p=edges(lcc);
  std::pair<iter_type, iter_type> p2=std::make_pair<iter_type, iter_type>(t1, t1);*/
  //p=std::make_pair(t1, t1); //p2;
 
 
  //iter_type(lcc), iter_type(lcc));
    //toto(lcc);

  //  std::pair<int, int> p;
  //  p = std::make_pair<int, int>(0,0);
  
  
  concept_check_polyhedron<LCC>();
  runtime_check_halfedgegraph<LCC>();
  std::cerr << "done\n";
  return 0;
}
