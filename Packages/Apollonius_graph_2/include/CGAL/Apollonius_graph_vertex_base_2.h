#ifndef CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H

#include <list>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE


template <class Gt, bool Store_trivial = true>
class Apollonius_graph_vertex_base_2
  :  private Triangulation_vertex_base_2<Gt>
{
public:
  typedef Gt                                    Geom_traits;
  typedef typename Gt::Weighted_point           Weighted_point;
  typedef Apollonius_graph_vertex_base_2<Gt>    Vertex_base;

  typedef std::list<Weighted_point>       Weighted_point_list;
  typedef Weighted_point_list::iterator   Weighted_point_list_iterator;

private:
  typedef Triangulation_vertex_base_2<Gt>       Vbase;

private:
  Weighted_point_list   weighted_point_list;

public:
  Apollonius_graph_vertex_base_2() : Vbase() {}

  Apollonius_graph_vertex_base_2(const Weighted_point& p, void* f = NULL) 
    : Vbase(p, f) {}

  ~Apollonius_graph_vertex_base_2()
  {
    clear_weighted_point_list();
  }

  inline
  Weighted_point point() const { return Vbase::point(); }

  inline void* face() const { return Vbase::face(); }

  inline void set_point(const Weighted_point& p) {
    Vbase::set_point(p);
  }

  inline void set_face(void* f) { Vbase::set_face(f); }

  inline bool is_valid(bool verbose, int level) const {
    return Vbase::is_valid(verbose, level);
  }

  inline
  void add_weighted_point(const Weighted_point & p)
  {
    if ( Store_trivial ) {
      weighted_point_list.push_back(p);
    }
  }

  inline
  unsigned int number_of_weighted_points() const
  {
    return weighted_point_list.size();
  }

  inline
  Weighted_point_list_iterator  weighted_points_begin()
  {
    return weighted_point_list.begin();
  }

  inline
  Weighted_point_list_iterator  weighted_points_end()
  {
    return weighted_point_list.end();
  }

  inline
  void clear_weighted_point_list()
  {
    weighted_point_list.clear();
  }

};

CGAL_END_NAMESPACE 

#endif // CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
