#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_3_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_3_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Delaunay_triangulation_visitor_base_3 {
  //typedef Tr Triangulation;
  Delaunay_triangulation_visitor_base_3(){}

  template <class Vertex_handle>
  void delete_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void new_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void change_vertex(Vertex_handle) {
  }

  template <class It>
  void new_cells(It, It) {
  }

  template <class It>
  void delete_cells(It, It) {
  }

  template <class Edge>
  void pre_edge_flip(Edge){

  }
  template <class Edge>
  void post_facet_flip(Edge){

  }

  template <class Facet>
  void pre_facet_flip(Facet){

  }
  
  template <class Facet>
  void post_edge_flip(Facet){
  }
};

CGAL_KDS_END_NAMESPACE

#endif
