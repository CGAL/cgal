#ifndef CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H

#include <list>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE


template <class Gt, bool StoreHidden = true>
class Apollonius_graph_vertex_base_2
  : private Triangulation_vertex_base_2<Gt>
{
public:
  // TYPES
  //------
  typedef Gt                               Geom_traits;
  typedef typename Gt::Weighted_point      Weighted_point;
private:
  // local types
  typedef std::list<Weighted_point>        Container;
  typedef Triangulation_vertex_base_2<Gt>  Vbase;
public:
  // TYPES (continued)
  //------------------
  typedef Container                        Hidden_weighted_point_container;
  typedef Container::iterator              Hidden_weighted_point_iterator;

public:
  // CREATION
  //---------
  Apollonius_graph_vertex_base_2() : Vbase() {}

  Apollonius_graph_vertex_base_2(const Weighted_point& p, void* f = NULL) 
    : Vbase(p, f) {}

  ~Apollonius_graph_vertex_base_2()
  {
    clear_hidden_weighted_point_container();
  }


  // ACCESS METHODS
  //---------------
  inline Weighted_point point() const { return Vbase::point(); }

  inline void* face() const { return Vbase::face(); }

  inline unsigned int number_of_hidden_weighted_points() const {
    return weighted_point_list.size();
  }

  inline
  Hidden_weighted_point_iterator hidden_weighted_points_begin() { 
    return weighted_point_list.begin();
  }

  inline
  Hidden_weighted_point_iterator hidden_weighted_points_end() {
    return weighted_point_list.end();
  }

public:
  // SETTING AND UNSETTING
  //----------------------
  inline void set_point(const Weighted_point& p) {
    Vbase::set_point(p);
  }

  inline void set_face(void* f) { Vbase::set_face(f); }

  inline
  void add_hidden_weighted_point(const Weighted_point & p)
  {
    if ( StoreHidden ) {
      weighted_point_list.push_back(p);
    }
  }

  inline
  void clear_hidden_weighted_point_container()
  {
    weighted_point_list.clear();
  }

public:
  // VALIDITY CHECK
  inline bool is_valid(bool verbose, int level) const {
    return Vbase::is_valid(verbose, level);
  }

private:
  // class variables
  Container weighted_point_list;
};

CGAL_END_NAMESPACE 

#endif // CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
