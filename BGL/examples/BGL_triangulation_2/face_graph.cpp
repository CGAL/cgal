#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/graph/filtered_graph.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       Epic;
typedef CGAL::Projection_traits_xy_3<Epic>                        K;
typedef K::Point_2                                                Point;

typedef CGAL::Triangulation_2<K>                                  Triangulation;

typedef boost::graph_traits<Triangulation>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator       vertex_iterator;
typedef boost::graph_traits<Triangulation>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<Triangulation>::halfedge_iterator     halfedge_iterator;
typedef boost::graph_traits<Triangulation>::edge_iterator         edge_iterator;
typedef boost::graph_traits<Triangulation>::edge_descriptor       edge_descriptor;
typedef boost::graph_traits<Triangulation>::face_descriptor       face_descriptor;
typedef boost::graph_traits<Triangulation>::face_iterator         face_iterator;

typedef std::map<vertex_descriptor, int>                          VertexIndexMap;
typedef boost::associative_property_map<VertexIndexMap>           VertexIdPropertyMap;

typedef boost::property_map<Triangulation, boost::vertex_point_t>::type Ppmap;

int main(int /*argc*/, char** /*argc*/)
{
  Triangulation tr;

  tr.insert(Point(0.1,0,1));
  tr.insert(Point(1,0,1));
  tr.insert(Point(0.2,0.2, 2));
  tr.insert(Point(0,1,2));
  tr.insert(Point(0,2,3));

  // Associate indices to the vertices
  VertexIndexMap vertex_id_map;
  VertexIdPropertyMap vertex_index_pmap(vertex_id_map);
  int index = 0;

  for(vertex_descriptor vd : vertices(tr))
    vertex_id_map[vd] = index++;
  std::cerr << index << " vertices" << std::endl;

  index = 0;
  for(halfedge_descriptor hd : halfedges(tr))
  {
    vertex_descriptor vd = source(hd, tr);
    CGAL_USE(vd);
    ++index;
  }
  std::cerr << index << " halfedges" << std::endl;

  index = 0;
  for(edge_descriptor ed : edges(tr))
  {
    vertex_descriptor vd = source(ed, tr);
    CGAL_USE(vd);
    ++index;
  }
  std::cerr << index << " edges" << std::endl;

  index = 0;
  for(face_descriptor fd : faces(tr))
  {
    halfedge_descriptor hd = halfedge(fd, tr);
    CGAL_USE(hd);
    ++index;
  }
  std::cerr << index << " faces" << std::endl;

  std::cerr << num_vertices(tr) << " " << num_halfedges(tr) << " " << num_edges(tr) << " " << num_faces(tr) << std::endl;

  std::cout << "vertices incident to the first vertex:" << std::endl;
  Ppmap ppmap = get(boost::vertex_point, tr);
  for(vertex_descriptor vd : vertices_around_target(*vertices(tr).first, tr))
    std::cout <<  ppmap[vd] << std::endl;

  ppmap[*(++vertices(tr).first)] = Point(78, 1, 2);
  std::cout << "changed point of vertex " << ppmap[*(++vertices(tr).first)] << std::endl;

  return EXIT_SUCCESS;
}
