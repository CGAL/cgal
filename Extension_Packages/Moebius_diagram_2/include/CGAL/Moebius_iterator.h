#ifndef CHR_MOEBIUS_ITERATOR_H
#define CHR_MOEBIUS_ITERATOR_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

CGAL_BEGIN_NAMESPACE

template <class Tds_>
class Moebius_circulator_3
{
 public:
  typedef Tds_ Tds;
  typedef Tds::Facet Facet;
  typedef Tds::Edge Edge;
  typedef Tds::Cell_handle Cell_handle;
  typedef Tds::Vertex_handle Vertex_handle;
  typedef Tds::Cell_circulator Cell_circulator;

  Moebius_circulator_3 ()
    : _s (), _t (),
    _prev (), _curr (), _next (),
    _fprev (), _fcurr (), _fnext (),
    _vprev (), _vcurr (), _vnext ()
    {}

  Moebius_circulator_3 (Cell_handle c, int s, int t)
    : _s(c->vertex(s)), _t(c->vertex(t)),
      _prev (c), _curr (next (_prev, _s, _t)), _next (next (_curr, _s, _t)),
      _fprev (facet (_prev, _s, _t)), _fcurr (facet (_curr, _s, _t)), 
      _fnext (facet (_next, _s, _t)),
      _vprev (vertex (_prev, _s, _t)), _vcurr (vertex (_curr, _s, _t)), 
      _vnext (vertex (_next, _s, _t)),
      
    {
      CGAL_triangulation_precondition( c != NULL &&
				       s >= 0 && s < 4 &&
				       t >= 0 && t < 4 );
    }

  Moebius_circulator_3 (const Edge & e)
    : _s (e.first->vertex(e.second)), _t (e.first->vertex(e.third)),
      _prev (e.first), _curr (next (_prev, _s, _t)), 
      _next (next (_curr, _s, _t)),
      _fprev (facet (_prev, _s, _t)), _fcurr (facet (_curr, _s, _t)), 
      _fnext (facet (_next, _s, _t)),
      _vprev (vertex (_prev, _s, _t)), _vcurr (vertex (_curr, _s, _t)), 
      _vnext (vertex (_next, _s, _t)),
    {
      CGAL_triangulation_precondition( e.first != NULL &&
				       e.second >=0 && e.second < 4 &&
				       e.third  >=0 && e.third  < 4);
    }
  
  //  Cell_handle prev_cell () const { return _prev; }
  Cell_handle cell () const { return _curr; }
  Cell_handle next_cell () const { return _next; }

  const Facet &prev_facet () const { return _fprev; }
  const Facet &facet () const { return _fcurr; }
  const Facet &next_facet () const { return _fnext; }

  Vertex_handle prev_vertex () const { return _vprev; }
  Vertex_handle vertex () const { return _vcurr; }
  Vertex_handle next_vertex () const { return _vnext; }

  Moebius_circulator_3 operator++ () {
    CGAL_triangulation_precondition( _prev != NULL );
    _prev = _curr; _curr = _next; _next = next (_curr, _s, _t);
    _fprev = _fcurr; _fcurr = _fnext; _fnext = facet (_next, _s, _t);
    _vprev = _vcurr; _vcurr = _vnext; _vnext = vertex (_next, _s, _t);
    return *this;
  }

 private:
  Vertex_handle _s;
  Vertex_handle _t;
  Cell_handle _prev, _curr, _next;
  Facet _fprev, _fcurr, _fnext;
  Vertex_handle _vprev, _vcurr, _vnext;

 protected:
  static Cell_handle next (const Cell_handle &c, 
			   const Vertex_handle &s, 
			   const Vertex_handle &t) {
    return c->neighbor (next_around_edge (c->index (s), c->index (t)));
  }

  static const Facet &facet (const Cell_handle &c, 
			     const Vertex_handle &s, 
			     const Vertex_handle &t) {
    return Facet (c, next_around_edge (c->index(s), pos->index(t)));
  }
    
  static Vertex_handle vertex (const Cell_handle &c, 
			       const Vertex_handle &s, 
			       const Vertex_handle &t) {
    return c->vertex (next_around_edge (c->index(t), pos->index(s)));
  }

  static int next_around_edge (const int i, const int j)
    {
      return Triangulation_utils_3::next_around_edge (i, j);
    }
};

CGAL_END_NAMESPACE

#endif// CHR_MOEBIUS_ITERATOR_H
