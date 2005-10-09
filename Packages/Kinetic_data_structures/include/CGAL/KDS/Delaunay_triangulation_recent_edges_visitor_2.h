#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_2_RE_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_2_RE_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>
#include <set>

CGAL_KDS_BEGIN_NAMESPACE

template <class Triangulation>
struct Delaunay_triangulation_recent_edges_visitor_2 {
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Vertex_handle VH;
  Delaunay_triangulation_recent_edges_visitor_2(){}

  void delete_vertex(VH) {
    recent_.clear();
  }
  void new_vertex(VH) {
    recent_.clear();
  }

  void change_vertex(VH vh) {
    recent_.clear();
    typename Triangulation::Edge_circulator ec= vh->incident_edges(), ef=ec;
    if (ec != NULL){
      do {
	recent_.insert(*ec);
	++ec;
      } while (ec != ef);
    }
  }

  template <class It>
  void new_faces(It, It) {
  }

  template <class It>
  void delete_faces(It, It) {
  }

  void pre_flip(Edge){
    recent_.clear();
  }
  void post_flip(Edge e){
    recent_.insert(e);
  }

  typedef typename std::set<Edge>::const_iterator iterator;
  iterator begin()  const {
    return recent_.begin();
  }
  iterator end()  const {
    return recent_.end();
  }

  bool contains(Edge e) const {
    return recent_.find(e) != recent_.end();
  }

  std::set<Edge> recent_;
};

CGAL_KDS_END_NAMESPACE

#endif
