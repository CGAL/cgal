// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_2_algorithms.C
// source        :
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#include <CGAL/Polygon_2_algorithms.h>

#include <cstdlib>
#include <algorithm>
#include <set>
#include <vector>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Simplicity_test_2
//-----------------------------------------------------------------------//
// The simplicity test is implemented as a class.
// The algorithm used is a sweep line algorithm. The sweep line is a
// horizontal line that sweeps from top (big y) to bottom.
// In the sweep status the edges are ordered from left (small) to right
// (big).

template <class ForwardIterator, class Traits>
class Simplicity_test_2 {
  protected:
    std::vector<ForwardIterator> d_index;
    // the attribute d_index is just a mapping between the integers and the
    // sequence of points

    int d_eventpoint;
    // the index of the current event point
    // the current sweepline is the horizontal line through this point

    const Traits& d_traits;
    // the traits class for polygons

  public:
    typedef typename Traits::Point_2 Point_2;

    Simplicity_test_2(const Traits& tr): d_traits(tr) {}
    ~Simplicity_test_2() {}

    const Traits& traits() const { return d_traits; }

    const Point_2& Vertex(int i) const { return *d_index[i]; }
    int NumberOfVertices() const { return d_index.size(); }

    const Point_2& EventPoint() const { return Vertex(d_eventpoint); }
    // return the current event point

    bool Test(ForwardIterator first, ForwardIterator last);
    // tests if the polygon with points in the range [first,last) is simple

    bool EdgesDoIntersect(int e1, int e2) const;
    // tests if the edges e1 and e2 have an intersection
    // N.B. the common vertex of two consecutive edges is not counted
    // as an intersection!

    bool VertexCompare(int i, int j) const;
    // compares the (lexicographical) order of vertex(i) and vertex(j)

    class VertexComp {
      protected:
        const Simplicity_test_2<ForwardIterator, Traits>* s;
      public:
        VertexComp() {}
        VertexComp(
          const Simplicity_test_2<ForwardIterator, Traits>* s0): s(s0)
        {}
        bool operator() (int i, int j) const { return s->VertexCompare(i,j); }
    };

    bool has_on_left_side(const Point_2& p,
                          const Point_2& q,
                          const Point_2& r ) const
    // returns true if point p is left of the point w, where w is the leftmost
    // intersection point of the horizontal line through p and the line
    // segment qr
    // N.B. if p lies on the segment qr, the result is indeterminate.
    {
      Comparison_result qr = d_traits.compare_y(q,r);
      if (qr == EQUAL)
        return (d_traits.compare_x(p,q) == SMALLER);
      else
        return ( d_traits.is_negative(d_traits.cross_product_2(p-q,r-q)) ==
                 (qr == SMALLER)                                       );
    }

    bool has_y_overlap(const Point_2& p,
                       const Point_2& q,
                       const Point_2& r ) const
    // returns true if the horizontal line through p intersects the segments qr
    {
      Comparison_result pq = d_traits.compare_y(p,q);
      Comparison_result pr = d_traits.compare_y(p,r);
      return (pq != pr) || (pq == EQUAL);
    }

    bool EdgeCompare(int e1, int e2) const;
    // computes the order of two edges e1 and e2 on the current sweepline

    bool edge_compare_consecutive(int e1, int e2) const;
    // computes the order of two edges e1 and e2 that share a vertex

    bool edge_compare_non_consecutive(int e1, int e2) const;
    // computes the order of two edges e1 and e2 that do not share a vertex

    bool consecutive_edges(int e1, int e2) const;
    // true if the edges e1 and e2 are not equal but share a vertex

    class EdgeComp {
      protected:
        const Simplicity_test_2<ForwardIterator, Traits>* s;
      public:
        EdgeComp() {}
        EdgeComp(
          const Simplicity_test_2<ForwardIterator, Traits>* s0): s(s0)
        {}
        bool operator() (int i, int j) const { return s->EdgeCompare(i,j); }
    };

    class EventQueue {
      //-----------------------------------------------------------------//
      // g++ 2.7.2 seems to have problems with the following typedef
      //
      // typedef set<int,VertexComp>::const_iterator const_iterator;
      //-----------------------------------------------------------------//

      protected:
        std::set<int,VertexComp> queue;
      public:
        EventQueue(Simplicity_test_2<ForwardIterator, Traits>* s)
          : queue(VertexComp(s)) {}
        bool insert(int i) { return queue.insert(i).second; }
        bool empty() const      { return queue.empty(); }
        int pop() {
          int Result = *(queue.begin());
          queue.erase(queue.begin());
          return Result;
        }
#ifdef CGAL_POLYGON_DEBUG
        void Show() const {
          cout << "    event queue: ";

          typename std::set<int,VertexComp>::const_iterator i;
          for (i = queue.begin(); i != queue.end(); ++i)
            cout << *i << " ";
          cout << endl;
        }
#endif
    };

    class SweepStatus {
      //-----------------------------------------------------------------//
      // g++ 2.7.2 seems to have problems with the following typedef
      //
      // typedef std::set<int,EdgeComp>::const_iterator const_iterator;
      //-----------------------------------------------------------------//

      protected:
        std::set<int,EdgeComp> status;
        // if i is an element of status, it means that 

        std::vector<typename std::set<int,EdgeComp>::iterator> index;
        // the iterators of the edges are stored to enable fast deletion

        const Simplicity_test_2<ForwardIterator, Traits>* s;
        // store a pointer to the Simplicity_test_2 class, to enable
        // access to the vertices

      public:
        SweepStatus(
          const Simplicity_test_2<ForwardIterator, Traits>* s0, int n)
          : status(EdgeComp(s0)), s(s0)
        {
	  index.insert(index.end(),n,status.end());
        }

        bool is_valid()
        // A necessary condition for the sweep status to be valid is that
        //
        // 1) every edge in the status intersects the current sweepline
        // 2) the edges are ordered along the current sweepline
        {
          int n = s->NumberOfVertices();

          typename std::set<int,EdgeComp>::const_iterator i;
          for (i = status.begin(); i != status.end(); ++i) {
            int v1 = *i;
            int v2 = (v1<n-1) ? v1+1 : v1+1-n; 
            // edge(v1) = (vertex(v1), vertex(v2))

            Comparison_result c1 =
              s->traits().compare_y(s->Vertex(v1), s->EventPoint());

            Comparison_result c2 =
              s->traits().compare_y(s->Vertex(v2), s->EventPoint());

            if (c1 == SMALLER && c2 == SMALLER) return false;
            if (c1 == LARGER && c2 == LARGER) return false;
          }

          return true;
        }

#ifdef CGAL_POLYGON_DEBUG
        void Show() {
          cout << "    sweep status: ";
          typename std::set<int,EdgeComp>::const_iterator i;
          for (i = status.begin(); i != status.end(); ++i)
            cout << *i << " ";
          cout << endl;
        }
#endif // CGAL_POLYGON_DEBUG

        void insert(int e)
        {
#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl << "    inserting edge " << e << " into sweep status" << endl;
}
#endif // CGAL_POLYGON_DEBUG

          index[e] = status.insert(e).first;
        }

        void erase(int e)
        {
#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl << "    removing  edge " << e << " from sweep status" << endl;
}
#endif // CGAL_POLYGON_DEBUG

          status.erase(index[e]);
        }

        int replace(int e1, int e2)
        {
#ifdef CGAL_POLYGON_DEBUG
{
    cout << endl << "    replacing edge " << e1 << " by edge "<< e2
      << " in sweep status" << endl;
}
#endif // CGAL_POLYGON_DEBUG
	    typename std::set<int,EdgeComp>::iterator cur = index[e1];
	    status.erase(cur++);
	    index[e2] = status.insert(cur, e2);
            return e2;
        }

        int left(int e) const
        { typename std::set<int,EdgeComp>::const_iterator i = index[e];
          return (i == status.begin()) ? -1 : *(--i);
        }

        int right(int e) const
        { typename std::set<int,EdgeComp>::const_iterator i = index[e]; ++i;
          return (i == status.end()) ? -1 : *i;
        }
    };
};

template <class ForwardIterator, class Traits>
inline
bool Simplicity_test_2<ForwardIterator, Traits>::VertexCompare(
  int i, int j) const
{
  return !d_traits.lexicographically_yx_smaller_or_equal(Vertex(i), Vertex(j));
}

template <class ForwardIterator, class Traits>
inline
bool Simplicity_test_2<ForwardIterator, Traits>::EdgeCompare(
  int e1, int e2) const
{
  // Edges must always be compared in the same order! This is to avoid problems
  // with overlapping edges:
  //
  //                             + v0
  //                            /
  //          -----------------+-v2----------------- sweepline
  //                          /
  //                         /
  //                        /
  //                       /
  //                      + v1
  //
  // In this case the order of the edges e0 = (v0,v1) and e1 = (v1,v2) on the
  // sweepline is indeterminate. However, it will only be detected that these
  // edges overlap after they have been inserted in the sweep state. To make
  // sure that this insertion is done correctly, the choice for the order
  // between e0 and e1 needs to be made consistently.

  bool Result;

  if (consecutive_edges(e1,e2)) {
    if (e1 < e2)
      Result = edge_compare_consecutive(e1,e2);
    else
      Result = !edge_compare_consecutive(e2,e1);
  }
  else {
    if (e1 < e2)
      Result = edge_compare_non_consecutive(e1,e2);
    else
      Result = !edge_compare_non_consecutive(e2,e1);
  }
  return Result;
}

template <class ForwardIterator, class Traits>
bool Simplicity_test_2<ForwardIterator, Traits>::edge_compare_consecutive(
  int e1, int e2) const
// This function is used to compare two edges that share a vertex:
//
//                             +                                |
//                            / \                               |
//                        e1 /   \ e2                           |
//                          /     \                             |
//                         /       \                            |
//           -------------+---------+-------------  sweepline   |
//                       /           \                          |
//                      +             +                         |
//
// Preconditions: 1) the shared vertex is lexicographically smaller or
//                   lexicographically bigger than both endpoints of the two
//                   edges (this condition is always satisfied in the
//                   sweepline algorithm)
//
//                2) both edges intersect the current sweepline

{
  int n = NumberOfVertices();
  int f1 = (e1<n-1) ? e1+1 : e1+1-n;  // edge(e1) = (vertex(e1), vertex(f1))
  int f2 = (e2<n-1) ? e2+1 : e2+1-n;  // edge(e2) = (vertex(e2), vertex(f2))

  if (f1 == e2)
    if (d_traits.compare_y(Vertex(e2), Vertex(f2)) != EQUAL)
       return has_on_left_side(Vertex(e1), Vertex(e2), Vertex(f2));
    else if (d_traits.compare_x(Vertex(e2), Vertex(f2)) == SMALLER)
    // Precondition 1) implies that e1 is on or above line (e2, f2).
    // If above the line, then say segment e1 is smaller.  If on the line,
    // vertex e1 is to the right of e2; by convention make the edge
    // with the smaller other endpoint the smaller one.
       return (d_traits.compare_y(Vertex(e1), Vertex(f2)) != EQUAL) ||
              (d_traits.compare_x(Vertex(e1), Vertex(f2)) == SMALLER);
    else // vertex e1 is on or below line (e2, f2); edge e1 is smaller only
         // if vertex e1 is on the line and to the left of f2
       return (d_traits.compare_y(Vertex(e1), Vertex(f2)) == EQUAL) &&
              (d_traits.compare_x(Vertex(e1), Vertex(f2)) == SMALLER);
  else // f2 and e1 are the same
    if (d_traits.compare_y(Vertex(e2), Vertex(f2)) != EQUAL)
       return has_on_left_side(Vertex(f1), Vertex(f2), Vertex(e2));
    else if (d_traits.compare_x(Vertex(e2), Vertex(f2)) == LARGER)
       return (d_traits.compare_y(Vertex(f1), Vertex(e2)) != EQUAL) ||
              (d_traits.compare_x(Vertex(f1), Vertex(e2)) == SMALLER);
    else
       return (d_traits.compare_y(Vertex(f1), Vertex(e2)) == EQUAL) &&
              (d_traits.compare_x(Vertex(f1), Vertex(e2)) == SMALLER);
}

template <class ForwardIterator, class Traits>
bool
Simplicity_test_2<ForwardIterator, Traits>::edge_compare_non_consecutive(
  int e1, int e2) const
{
  int n = NumberOfVertices();
  int f1 = (e1<n-1) ? e1+1 : e1+1-n;  // edge(e1) = (vertex(e1), vertex(f1))
  int f2 = (e2<n-1) ? e2+1 : e2+1-n;  // edge(e2) = (vertex(e2), vertex(f2))

  if (has_y_overlap(Vertex(e1), Vertex(e2), Vertex(f2)))
    return has_on_left_side(Vertex(e1), Vertex(e2), Vertex(f2));

  if (has_y_overlap(Vertex(f1), Vertex(e2), Vertex(f2)))
    return has_on_left_side(Vertex(f1), Vertex(e2), Vertex(f2));

  // if the vertices from edge e1 do not have y-overlap with edge e2, then
  // the vertices from edge e2 must have y_overlap with edge e1
  return !has_on_left_side(Vertex(e2), Vertex(e1), Vertex(f1));
}

template <class ForwardIterator, class Traits>
bool
Simplicity_test_2<ForwardIterator, Traits>::Test(ForwardIterator first,
                                                      ForwardIterator last)
{
  int n = 0;

  EventQueue events(this);
  while (first != last) {
    d_index.push_back(first);
    if (!events.insert(n++)) // if two vertices coincide...
      return false;
    ++first;
  }

  if (d_index.size() < 3)
    return true;

#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl;
  cout << "--- Simplicity test ----------------------------" << endl;
  cout << endl;
  cout << "Vertices:" << endl;
  typedef typename std::vector<ForwardIterator>::size_type Size_type;
  Size_type i;
  for (i=0; i<d_index.size(); i++)
    cout << i << " " << Vertex(i) << endl;
  cout << endl;
  events.Show();
}
#endif // CGAL_POLYGON_DEBUG

  SweepStatus status(this,n);

#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl;
  status.Show();
}
#endif // CGAL_POLYGON_DEBUG

  while (!events.empty()) {

#ifdef CGAL_POLYGON_DEBUG
  CGAL_polygon_assertion(status.is_valid());
#endif // CGAL_POLYGON_DEBUG

    int i = events.pop();
    d_eventpoint = i;

#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl;
  cout << "--- Event point: " << d_eventpoint << " on sweep line {y = "
       << Vertex(d_eventpoint).y() << "}" << endl;
  cout << endl;
  status.Show();
}
#endif // CGAL_POLYGON_DEBUG

    int prev = (i>0) ? i-1 : i-1+n;
    int next = (i<n-1) ? i+1 : i+1-n;

    bool prev_less_than_i = VertexCompare(i,prev);
    bool next_less_than_i = VertexCompare(i,next);
    if (prev_less_than_i != next_less_than_i) {
      int e = prev_less_than_i ? status.replace(i,prev) :
        status.replace(prev,i);
      CGAL_polygon_assertion(status.is_valid());

#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl;
  status.Show();
}
#endif // CGAL_POLYGON_DEBUG

      // check for intersections of newly inserted edge e with neighbors
      int left = status.left(e);
      if ((left >= 0) && (EdgesDoIntersect(left,e))) return false;

      int right = status.right(e);
      if ((right >= 0) && (EdgesDoIntersect(e,right))) return false;
    }
    else if (prev_less_than_i) {
      int e1 = prev;
      int e2 = i;

      status.insert(e1);
      status.insert(e2);
      CGAL_polygon_assertion(status.is_valid());

#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl;
  status.Show();
}
#endif // CGAL_POLYGON_DEBUG

      // check for intersections of newly inserted edges e1 and e2 with
      // neighbors
      int left, right;
      left = status.left(e1);
      if ((left >= 0) && (EdgesDoIntersect(left,e1))) return false;

      right = status.right(e1);
      if ((right >= 0) && (EdgesDoIntersect(e1,right))) return false;

      left = status.left(e2);
      if ((left >= 0) && (left != e1) && (EdgesDoIntersect(left,e2)))
        return false;

      right = status.right(e2);
      if ((right >= 0) && (right != e1) &&(EdgesDoIntersect(e2,right)))
        return false;
    }
    else {
       // check for intersections between edges that become new neighbors 
       // in the sweep status due to the deletion 
       int left, right;
       if (status.left(prev) == i)
       {
          left = status.left(i);
          right = status.right(prev);
       }
       else
       {
          left = status.left(prev);
          right = status.right(i);
       }
       status.erase(prev);
       status.erase(i);
       if (left >=0 && right >=0 && EdgesDoIntersect(left, right))
          return false;
      CGAL_polygon_assertion(status.is_valid());

#ifdef CGAL_POLYGON_DEBUG
{
  cout << endl;
  status.Show();
}
#endif // CGAL_POLYGON_DEBUG
    }
  }

  return true;
}

template <class ForwardIterator, class Traits>
bool Simplicity_test_2<ForwardIterator, Traits>::EdgesDoIntersect(
  int e1, int e2) const
{
#ifdef CGAL_POLYGON_DEBUG
{
  cout << "      intersecting edges " << e1 << " and " << e2 << endl;
}
#endif

  int n = NumberOfVertices();
  int f1 = (e1<n-1) ? e1+1 : e1+1-n;  // edge(e1) = (vertex(e1), vertex(f1))
  int f2 = (e2<n-1) ? e2+1 : e2+1-n;  // edge(e2) = (vertex(e2), vertex(f2))

  bool Result;
  if (consecutive_edges(e1,e2))
    Result = d_traits.have_equal_direction(Vertex(f1) - Vertex(e1),
                                           Vertex(e2) - Vertex(f2) );
  else
    Result = d_traits.do_intersect(Vertex(e1),
                                   Vertex(f1),
                                   Vertex(e2),
                                   Vertex(f2));

  // N.B. traits() instead of d_traits gives an error with g++ 2.7.2

  return Result;
}

template <class ForwardIterator, class Traits>
inline
bool Simplicity_test_2<ForwardIterator, Traits>::consecutive_edges(
  int e1, int e2) const
{
  int n = NumberOfVertices();
  return ( abs(e2-e1) == 1 || abs(e2-e1) == n-1 );
}

//-----------------------------------------------------------------------//
//                          is_simple_2
//-----------------------------------------------------------------------//
// uses Traits::compare_x
//      Traits::compare_y
//      Traits::cross_product_2
//      Traits::do_intersect
//      Traits::have_equal_direction
//      Traits::is_negative
//      Traits::lexicographically_yx_smaller_or_equal

template <class ForwardIterator, class Traits>
bool is_simple_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits)
{
  Simplicity_test_2<ForwardIterator, Traits> test(traits);
  return test.Test(first, last);
}

//-----------------------------------------------------------------------//
//                          left_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy

template <class ForwardIterator, class Traits>
ForwardIterator left_vertex_2(ForwardIterator first,
                                   ForwardIterator last,
                                   const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_xy Less_xy;
  return std::min_element(first, last, Less_xy());
}

//-----------------------------------------------------------------------//
//                          right_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy

template <class ForwardIterator, class Traits>
ForwardIterator right_vertex_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_xy Less_xy;
  return std::max_element(first, last, Less_xy());
}

//-----------------------------------------------------------------------//
//                          top_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_yx

template <class ForwardIterator, class Traits>
ForwardIterator top_vertex_2(ForwardIterator first,
                                  ForwardIterator last,
                                  const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_yx Less_yx;
  return std::max_element(first, last, Less_yx());
}

//-----------------------------------------------------------------------//
//                          bottom_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_yx

template <class ForwardIterator, class Traits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
                                     ForwardIterator last,
                                     const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_yx Less_yx;
  return std::min_element(first, last, Less_yx());
}

//-----------------------------------------------------------------------//
//                          bbox_2
//-----------------------------------------------------------------------//

template <class InputIterator>
Bbox_2 bbox_2(InputIterator first, InputIterator last)
{
  CGAL_polygon_precondition(first != last);
  Bbox_2 result = (*first).bbox();

  while (++first != last)
    result = result + (*first).bbox();

  return result;
}

//-----------------------------------------------------------------------//
//                          area_2
//-----------------------------------------------------------------------//
// uses Traits::determinant_2

// template <class ForwardIterator, class FT, class Traits>
// void area_2(ForwardIterator first,
//                 ForwardIterator last,
//                 FT& result,
//                 const Traits& traits)
//{

//-----------------------------------------------------------------------//
//                          is_convex_2
//-----------------------------------------------------------------------//
// uses Traits::lexicographically_xy_smaller
//      Traits::orientation

template <class ForwardIterator, class Traits>
bool is_convex_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits)
{
  ForwardIterator previous = first;
  if (previous == last) return true;

  ForwardIterator current = previous; ++current;
  if (current == last) return true;

  ForwardIterator next = current; ++next;
  if (next == last) return true;

  // initialization
  bool HasClockwiseTriples = false;
  bool HasCounterClockwiseTriples = false;
  bool Order = traits.lexicographically_xy_smaller(*previous, *current);
  int NumOrderChanges = 0;

  do {
    switch (traits.orientation(*previous, *current, *next)) {
      case CLOCKWISE:
        HasClockwiseTriples = true;
        break;
      case COUNTERCLOCKWISE:
        HasCounterClockwiseTriples = true;
        break;
      default:
	;
    }

    bool NewOrder = traits.lexicographically_xy_smaller(*current, *next);
    if (Order != NewOrder) NumOrderChanges++;

    if (NumOrderChanges > 2) {
#ifdef CGAL_POLYGON_DEBUG
cout << "too many order changes: not convex!" << endl;
#endif
      return false;
    }

    if (HasClockwiseTriples && HasCounterClockwiseTriples) {
#ifdef CGAL_POLYGON_DEBUG
cout << "polygon not locally convex!" << endl;
#endif
      return false;
    }

    previous = current;
    current = next;
    ++next;
    if (next == last) next = first;
    Order = NewOrder;
  }
  while (previous != first);

  return true;
}

//-----------------------------------------------------------------------//
//                          oriented_side_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy
//      Traits::compare_x
//      Traits::compare_y
//      Traits::determinant_2
//      Traits::orientation
//      Traits::sign

template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point& point,
                                        const Traits& traits)
{
  Oriented_side result;

  Orientation o = orientation_2(first, last, traits);
  CGAL_polygon_assertion(o != COLLINEAR);

  Bounded_side b = bounded_side_2(first, last, point, traits);
  switch (b) {
    case ON_BOUNDARY:
      result = ON_ORIENTED_BOUNDARY;
      break;

    case ON_BOUNDED_SIDE:
      result = (o == CLOCKWISE) ?
               ON_NEGATIVE_SIDE :
               ON_POSITIVE_SIDE;
      break;

    case ON_UNBOUNDED_SIDE:
      result = (o == CLOCKWISE) ?
               ON_POSITIVE_SIDE :
               ON_NEGATIVE_SIDE;
      break;
  }

  return result;
}

//-----------------------------------------------------------------------//
//                          bounded_side_2
//-----------------------------------------------------------------------//
// uses Traits::compare_x
//      Traits::compare_y
//      Traits::determinant_2
//      Traits::sign
//
// returns ON_BOUNDED_SIDE, ON_BOUNDARY or ON_UNBOUNDED_SIDE

template <class ForwardIterator, class Point, class Traits>
Bounded_side bounded_side_2(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point& point,
                                      const Traits& traits)
{
  ForwardIterator current = first;
  if (current == last) return ON_UNBOUNDED_SIDE;

  ForwardIterator next = current; ++next;
  if (next == last) return ON_UNBOUNDED_SIDE;

  bool IsInside = false;
  Comparison_result CompareCurrent = traits.compare_y(*current, point);

  do // check if the segment (current,next) intersects
     // the ray { (t,y) | t >= point.x() }
  {
    Comparison_result CompareNext = traits.compare_y(*next, point);

    switch (CompareCurrent) {
      case SMALLER:
        switch (CompareNext) {
          case SMALLER:
            break;
          case EQUAL:
            switch (traits.compare_x(point, *next)) {
              case SMALLER: IsInside = !IsInside; break;
              case EQUAL:   return ON_BOUNDARY;
              case LARGER:  break;
            }
            break;
          case LARGER:
            if (point.x() < std::min((*current).x(), (*next).x())) {
              IsInside = !IsInside;
            }
            else if (point.x() <= std::max((*current).x(),(*next).x())) {
              switch (traits.sign(traits.determinant_2(point,
                                                       *current,
                                                       *next)))
              {
                case 0: return ON_BOUNDARY;
                case 1: IsInside = !IsInside; break;
              }
            }
            break;
        }
        break;
      case EQUAL:
        switch (CompareNext) {
          case SMALLER:
            switch (traits.compare_x(point, *current)) {
              case SMALLER: IsInside = !IsInside; break;
              case EQUAL:   return ON_BOUNDARY;
              case LARGER:  break;
            }
            break;
          case EQUAL:
            if ( (std::min((*current).x(), (*next).x()) <= point.x()) &&
                 (point.x() <= std::max((*current).x(), (*next).x()))    ) {
              return ON_BOUNDARY;
            }
            break;
          case LARGER:
            if (point.x() == (*current).x()) {
              return ON_BOUNDARY;
            }
            break;
        }
        break;
      case LARGER:
        switch (CompareNext) {
          case SMALLER:
            if (point.x() < std::min((*current).x(), (*next).x())) {
              IsInside = !IsInside;
            }
            else if (point.x() <= std::max((*current).x(),(*next).x())) {
              switch (traits.sign(traits.determinant_2(point,
                                                       *current,
                                                       *next)))
              {
                case -1: IsInside = !IsInside; break;
                case  0: return ON_BOUNDARY;
              }
            }
            break;
          case EQUAL:
            if (point.x() == (*next).x()) {
              return ON_BOUNDARY;
            }
            break;
          case LARGER:
            break;
        }
        break;
    }

    current = next;
    CompareCurrent = CompareNext;
    ++next;
    if (next == last) next = first;   
  }
  while (current != first);

  return IsInside ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
}

//-----------------------------------------------------------------------//
//                          orientation_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy
//      Traits::orientation

template <class ForwardIterator, class Traits>
Orientation orientation_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits& traits)
{
  CGAL_polygon_precondition(is_simple_2(first, last, traits));

  ForwardIterator i = left_vertex_2(first, last, traits);

  ForwardIterator prev = (i == first) ? last : i;
  --prev;

  ForwardIterator next = i;
  ++next;
  if (next == last)
    next = first;

  // if the range [first,last) contains less than three points, then some
  // of the points (prev,i,next) will coincide

  // return the orientation of the triple (prev,i,next)
  return traits.orientation(*prev, *i, *next);
}

CGAL_END_NAMESPACE

