#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Delaunay_triangulation_visitor_base_2 {
  Delaunay_triangulation_visitor_base_2(){}
  template <class Ok>
  void delete_vertex(Ok) {
  }

  template <class Ok>
  void new_vertex(Ok) {
  }

  template <class Ok>
  void change_vertex(Ok) {
  }

  template <class It>
  void new_faces(It, It) {
  }

  template <class It>
  void delete_faces(It, It) {
  }

  template <class E>
  void pre_flip(E){

  }
  template <class E>
  void post_flip(E){

  }
};

CGAL_KDS_END_NAMESPACE

#endif
