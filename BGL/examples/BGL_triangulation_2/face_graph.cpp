#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/graph/filtered_graph.hpp>
#include <boost/foreach.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Projection_traits_xy_3<Epic>  K;
typedef K::Point_2 Point;

typedef CGAL::Triangulation_2<K> Triangulation;


typedef boost::graph_traits<Triangulation>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Triangulation>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Triangulation>::halfedge_iterator halfedge_iterator;
typedef boost::graph_traits<Triangulation>::face_descriptor face_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Triangulation>::face_iterator face_iterator;
typedef boost::graph_traits<Triangulation>::edge_iterator edge_iterator;
typedef boost::graph_traits<Triangulation>::edge_descriptor edge_descriptor;

typedef std::map<vertex_descriptor,int> VertexIndexMap;
VertexIndexMap vertex_id_map;

typedef boost::associative_property_map<VertexIndexMap> VertexIdPropertyMap;
VertexIdPropertyMap vertex_index_pmap(vertex_id_map);

typedef std::map<halfedge_descriptor,int> HalfedgeIndexMap;
HalfedgeIndexMap halfedge_id_map;

typedef boost::associative_property_map<HalfedgeIndexMap> HalfedgeIdPropertyMap;
HalfedgeIdPropertyMap halfedge_index_pmap(halfedge_id_map);

int
main(int,char*[])
{
  Triangulation t;

  t.insert(Point(0.1,0,1));
  t.insert(Point(1,0,1));
  t.insert(Point(0.2,0.2, 2));
  t.insert(Point(0,1,2));
  t.insert(Point(0,2,3));

  vertex_iterator vit, ve;
  // Associate indices to the vertices
  int index = 0;
  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vit and ve
  for(boost::tie(vit,ve) = vertices(t); vit!=ve; ++vit ){
    vertex_descriptor  vd = *vit;
    if(! t.is_infinite(vd)){
      vertex_id_map[vd]= index++;
    }
  }

  std::cerr << index << " vertices" << std::endl;
  index = 0;
  face_iterator fit,fe;
  for(boost::tie(fit,fe) = faces(t); fit!= fe; ++fit){
    face_descriptor fd = *fit;
    halfedge_descriptor hd = halfedge(fd,t);
    halfedge_descriptor n = next(hd,t);
    
    halfedge_descriptor nn = next(n,t);
    if(next(nn,t) != hd){
      std::cerr << "the face is not a triangle" << std::endl;
    }
    
    ++index;
  }
  
  std::cerr << index << " faces" << std::endl;
  index = 0;

  edge_iterator eit,ee;
  for(boost::tie(eit,ee) = edges(t); eit!= ee; ++eit){
    edge_descriptor ed = *eit;
    vertex_descriptor vd = source(ed,t);
    CGAL_USE(vd);
    ++index;
  }

  std::cerr << index << " edges" << std::endl;
  index = 0;

  halfedge_iterator hit,he;
  for(boost::tie(hit,he) = halfedges(t); hit!= he; ++hit){
    halfedge_descriptor hd = *hit;
    vertex_descriptor vd = source(hd,t);
    CGAL_USE(vd);
    ++index;
  }
  std::cerr << index << " halfedges" << std::endl;

  std::cerr << num_vertices(t) << " " << num_edges(t) << " " << num_halfedges(t) << " " << num_faces(t) << std::endl;

  typedef boost::property_map<Triangulation, boost::vertex_point_t>::type Ppmap;
  Ppmap ppmap = get(boost::vertex_point, t);
 

  BOOST_FOREACH(vertex_descriptor vd, vertices_around_target(*vertices(t).first, t)){
    std::cout <<  ppmap[vd] << std::endl;
  }


  ppmap[*(++vertices(t).first)] = Point(78,1,2);
  std::cout << " changed point of vertex " << ppmap[*(++vertices(t).first)] << std::endl;

  return 0;
}
