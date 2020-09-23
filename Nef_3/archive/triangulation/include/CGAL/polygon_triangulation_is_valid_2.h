#ifndef POLYGON_TRIANGULATION_IS_VALID_2_H
#define POLYGON_TRIANGULATION_IS_VALID_2_H

#include <map>
#include <list>

#undef _DEBUG
#define _DEBUG 17
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <class Traits>
struct Less_xy_segment_2
{
  template <class Segment_2>
  bool operator()( const Segment_2& s, const Segment_2& r) {
    typedef typename Traits::Compare_xy_2 Compare_xy_2;
    Compare_xy_2 compare_xy_2;
    return ((compare_xy_2( s.source(), r.source()) == EQUAL) ?
            (compare_xy_2( s.target(), r.target()) == SMALLER) :
            (compare_xy_2( s.source(), r.source()) == SMALLER));
  }
};

template <class Traits>
class Segment_lexicographically_xy_ordered_2
{
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Compare_xy_2 Compare_xy_2;
 public:
  //Segment_lexicographically_xy_ordered_2() : p(), q() {}
  Segment_lexicographically_xy_ordered_2( const Point_2& pi,
                                          const Point_2& qi) {
    Compare_xy_2 compare_xy_2;
    if( compare_xy_2( pi, qi) == SMALLER) {
      p = pi;
      q = qi;
    }
    else {
      p = qi;
      q = pi;
    }
  }
  Point_2 source() const { return p; }
  Point_2 target() const { return q; }
 private:
  Point_2 p, q;
};

template <class BoundaryIterator, class TriangleIterator, class Traits>
bool is_polygon_triangulation_valid_2( BoundaryIterator boundary_begin,
                                       BoundaryIterator boundary_end,
                                       TriangleIterator triangle_begin,
                                       TriangleIterator triangle_end,
                                       const Traits& traits) {
  typedef typename BoundaryIterator::value_type Circulator;
  typedef Segment_lexicographically_xy_ordered_2<Traits> Segment_2;
  typedef std::map< Segment_2, int, Less_xy_segment_2<Traits> > Counter;
  typedef typename Counter::iterator Counter_iterator;

  Counter counter;
  // count the number of faces incident to every edge
  for( TriangleIterator ti = triangle_begin; ti != triangle_end; ++ti) {
    counter[Segment_2( ti->vertex(0), ti->vertex(1))]++;
    counter[Segment_2( ti->vertex(1), ti->vertex(2))]++;
    counter[Segment_2( ti->vertex(2), ti->vertex(0))]++;
  }
  // check that boundary edges have only one incident face
  // and increment their incident faces counter by one
  BoundaryIterator bi = boundary_begin;
  Circulator c(*bi), cprev(c), cend(c);
  cprev--;
  CGAL_For_all( c, cend) { // outer boundary
    Segment_2 s( *cprev, *c);
    if( counter[s] != 1)
      return false;
    counter[s]++;
    cprev = c;
  }
  for( ++bi; bi != boundary_end; ++bi) {
    Circulator c(*bi), cprev(c), cend(c);
    cprev--;
    CGAL_For_all( c, cend) { // inner boundary
      Segment_2 s( *c, *cprev);
      if( counter[s] != 1)
        return false;
      counter[s]++;
      cprev = c;
    }
  }
  // verify that every edge has then two incident faces
  for( Counter_iterator ci = counter.begin(); ci != counter.end(); ++ci) {
    if( ci->second != 2)
      return false;
  }
  return true;
}

}

#endif // POLYGON_TRIANGULATION_IS_VALID_2_H
