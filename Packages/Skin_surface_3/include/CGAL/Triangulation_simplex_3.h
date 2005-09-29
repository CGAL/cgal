#ifndef TRIANGULATION_SIMPLEX_3_H
#define TRIANGULATION_SIMPLEX_3_H

CGAL_BEGIN_NAMESPACE

template < class Triangulation >
class Triangulation_simplex_3 {
  typedef Triangulation             T;
  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Edge          Edge;
  typedef typename T::Facet         Facet;
  typedef typename T::Cell_handle Cell_handle;
public:
  Triangulation_simplex_3() : ref(-1), ch() {}
	
  Triangulation_simplex_3(Vertex_handle v) {
    ch = v->cell();
    ref = (ch->index(v) << 2); /* dim == 0 */
    CGAL_assertion (ch != NULL);
  }
  Triangulation_simplex_3(Edge e) {
    ch = e.first;
    ref = (((e.third<< 2) + e.second) << 2) + 1; /* dim */
    CGAL_assertion (ch != NULL);
  }
  Triangulation_simplex_3(Facet f) {
    ch = f.first;
    ref = (f.second << 2) + 2; /* dim */
    CGAL_assertion (ch != NULL);
  }
  Triangulation_simplex_3(Cell_handle c) {
    ch = c;
    ref = 3; /* dim */
    CGAL_assertion (ch != NULL);
  }

  operator Vertex_handle ()
  {
    assert(dimension() == 0);
    return ch->vertex(index(0));
  }

  operator Edge ()
  {
    assert(dimension() == 1);
    return Edge(ch,index(0),index(1));
  }

  operator Facet ()
  {
    assert(dimension() == 2);
    return Facet(ch,index(0));
  }

  operator Cell_handle ()
  {
    assert(dimension() == 3);
    return ch;
  }
  // returns the dimension of the simplex
  int dimension () const {
    return (ref & 3);
  }
  template < class T >
  friend bool operator==(Triangulation_simplex_3<T> s0, Triangulation_simplex_3<T> s1);
  template < class T >
  friend bool operator< (Triangulation_simplex_3<T> s0, Triangulation_simplex_3<T> s1);
	
private:
  int ref;        // storage iijjdd (index i, index j, dimension of simplex)
  Cell_handle ch; // Corresponding cell handle

  inline int index(int i) const {
    CGAL_assertion (i==0 || ((i==1) && (dimension()==1)));
    return (ref >> (2*(i+1))) & 3;
  }
};

///////////////////////////////
// Simplex functions
///////////////////////////////
template < class T >
bool
operator==(Triangulation_simplex_3<T> s0, Triangulation_simplex_3<T> s1) {
  typedef Triangulation_simplex_3<T>             Sim;
  if (s0.dimension() != s1.dimension()) return false;
	
  typename Sim::Cell_handle neighbor;
	
  switch (s0.dimension()) {
  case (0): // Vertex
    return (s0.ch->vertex(s0.index(0)) == s1.ch->vertex(s1.index(0)));
  case (1): // Edge
    return ((s0.ch->vertex(s0.index(0)) == s1.ch->vertex(s1.index(0)) &&
	     s0.ch->vertex(s0.index(1)) == s1.ch->vertex(s1.index(1))) ||
	    (s0.ch->vertex(s0.index(1)) == s1.ch->vertex(s1.index(0)) &&
	     s0.ch->vertex(s0.index(0)) == s1.ch->vertex(s1.index(1))));
  case (2):
    if (s0.ch == s1.ch && s0.index(0) == s1.index(0)) {
      return true;
    }
			
    neighbor = s0.ch->neighbor(s0.index(0));
    if (neighbor == s1.ch &&
	neighbor->index(s0.ch) == s1.index(0)) {
      return true;
    }
    return false;
  case (3):
    return (&(*s0.ch) == &(*s1.ch));
  }
  CGAL_assertion(false);
  return false;
}

template < class T >
bool
operator<(Triangulation_simplex_3<T> s0, Triangulation_simplex_3<T> s1) {
  typedef Triangulation_simplex_3<T>             Sim;

  if (s0 == s1) return false;
  if (s0.dimension() < s1.dimension()) return true;
  if (s0.dimension() > s1.dimension()) return false;
	
  // Dimensions are equal, compare the memory addresses of the simplices
  typename Sim::Cell_handle ch1, ch2;
  typename Sim::Vertex_handle vh1, vh2, vh3, vh4;
  switch (s0.dimension()) {
  case (0): // Vertex
    // Vertextices are not equal
    return (&(*s0.ch->vertex(s0.index(0))) <
	    &(*s1.ch->vertex(s1.index(0))));
  case (1): // Edge
    vh1 = s0.ch->vertex(s0.index(0));
    vh2 = s0.ch->vertex(s0.index(1));
    vh3 = s1.ch->vertex(s1.index(0));
    vh4 = s1.ch->vertex(s1.index(1));
			
    if (std::min(&(*vh1), &(*vh2)) < std::min(&(*vh3), &(*vh4)))
      return true;
			
    if (std::min(&(*vh1), &(*vh2)) > std::min(&(*vh3), &(*vh4)))
      return false;
			
    if (std::max(&(*vh1), &(*vh2)) < std::max(&(*vh3), &(*vh4)))
      return true;
			
    return false;
  case (2): // Facet
    ch1 = s0.ch->neighbor(s0.index(0));
    ch2 = s1.ch->neighbor(s1.index(0));
			
    if (std::min(&(*s0.ch), &(*ch1)) < std::min(&(*s1.ch), &(*ch2)))
      return true;
			
    if (std::min(&(*s0.ch), &(*ch1)) > std::min(&(*s1.ch), &(*ch2)))
      return false;
			
    if (std::max(&(*s0.ch), &(*ch1)) < std::max(&(*s1.ch), &(*ch2)))
      return true;
			
    return false;
  case (3): // Cell
    return (&(*s0.ch) < &(*s1.ch));
  }
  CGAL_assertion(0);
  return false;
}

CGAL_END_NAMESPACE

#endif // TRIANGULATION_SIMPLEX_3_H
