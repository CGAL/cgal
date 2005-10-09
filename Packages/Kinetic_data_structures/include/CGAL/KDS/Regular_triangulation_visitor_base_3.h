#ifndef CGAL_KDS_REGULAR_TRIANGULATION_3_WATCHER_BASE_H
#define CGAL_KDS_REGULAR_TRIANGULATION_3_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Delaunay_triangulation_visitor_base_3.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Regular_triangulation_visitor_base_3: 
  public Delaunay_triangulation_visitor_base_3 {
  //typedef Tr Triangulation;
  Regular_triangulation_visitor_base_3(){}

  template <class Key, class Cell>
  void pre_move(Key, Cell){}

  template <class Key, class Cell>
  void post_move(Key, Cell){}
  
  template <class Key, class Cell>
  void pre_push(Key, Cell){}

  template <class Vertex_handle>
  void post_push(Vertex_handle){}

  template <class Vertex_handle>
  void pre_pop(Vertex_handle){}

  template <class Key, class Cell>
  void post_pop(Key, Cell){}
};

CGAL_KDS_END_NAMESPACE

#endif
