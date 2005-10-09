#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_HELPER_2_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_HELPER_2_H
#include <CGAL/KDS/basic.h>
#include <utility>
CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class TDS>
struct Triangulation_data_structure_helper_2 {
  typedef typename TDS::Edge Edge;
  typedef typename TDS::Vertex_handle Vertex_handle;
  typedef typename TDS::Face Face;
  typedef typename Face::Edge_label Edge_label;

  static Vertex_handle origin(const Edge &e){
    int o= e.first->ccw(e.second);
    return e.first->vertex(o);
  }

  static Vertex_handle destination(const Edge &e){
    int o= e.first->cw(e.second);
    return e.first->vertex(o);
  }

  static Edge_label get_directed_edge_label(const Edge &e){
    return e.first->get_edge_label(e.second);
  }

 static Edge_label get_undirected_edge_label(const Edge &e){
   if (get_directed_edge_label(e) != get_directed_edge_label(mirror_edge(e))){
     std::cerr << get_directed_edge_label(e) << " " 
	       <<  get_directed_edge_label(mirror_edge(e)) << std::endl;
     CGAL_precondition(get_directed_edge_label(e) == get_directed_edge_label(mirror_edge(e)));
   }
   return e.first->get_edge_label(e.second);
  }

  static Edge mirror_edge(const Edge &e){
    int i= e.first->mirror_index(e.second);
    return Edge(e.first->neighbor(e.second), i);
  }
  
  static void set_directed_edge_label(const Edge &e, const Edge_label& l){
    e.first->set_edge_label(e.second, l);
  }
  
  static void set_undirected_edge_label(const Edge &e, const Edge_label& l){
    set_directed_edge_label(e,l);
    set_directed_edge_label(mirror_edge(e),l);
  }

  static Vertex_handle third_vertex(const Edge &e){
    return e.first->vertex(e.second);
  }
  
  static Vertex_handle mirror_vertex(const Edge &e){
    return third_vertex(mirror_edge(e));
  }
};

CGAL_KDS_END_INTERNAL_NAMESPACE
#endif
