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

#ifndef CGAL_POLYGON_2_ALGORITHMS_H
#include <CGAL/Polygon_2_algorithms.h>
#endif // CGAL_POLYGON_2_ALGORITHMS_H

#ifndef CGAL_PROTECT_CSTDLIB
#include <cstdlib>
#define CGAL_PROTECT_CSTDLIB
#endif
#ifndef CGAL_PROTECT_ALGORITHM
#include <algorithm>
#define CGAL_PROTECT_ALGORITHM
#endif
#ifndef CGAL_PROTECT_SET
#include <set>
#define CGAL_PROTECT_SET
#endif
#ifndef CGAL_PROTECT_VECTOR
#include <vector>
#define CGAL_PROTECT_VECTOR
#endif

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          CGAL_Simplicity_test_2
//-----------------------------------------------------------------------//
// The simplicity test is implemented as a class.

template <class ForwardIterator, class Traits>
class CGAL_Simplicity_test_2 {
  private:
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

    CGAL_Simplicity_test_2(const Traits& tr): d_traits(tr) {}
    ~CGAL_Simplicity_test_2() {}

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
      private:
        const CGAL_Simplicity_test_2<ForwardIterator, Traits>* s;
      public:
        VertexComp() {}
        VertexComp(
          const CGAL_Simplicity_test_2<ForwardIterator, Traits>* s0): s(s0)
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
      CGAL_Comparison_result qr = d_traits.compare_y(q,r);
      if (qr == CGAL_EQUAL)
        return (d_traits.compare_x(p,q) == CGAL_SMALLER);
      else
        return ( d_traits.is_negative(d_traits.cross_product_2(p-q,r-q)) ==
                 (qr == CGAL_SMALLER)                                       );
    }

    bool has_y_overlap(const Point_2& p,
                       const Point_2& q,
                       const Point_2& r ) const
    // returns true if the horizontal line through p intersects the segments qr
    {
      CGAL_Comparison_result pq = d_traits.compare_y(p,q);
      CGAL_Comparison_result pr = d_traits.compare_y(p,r);
      return (pq != pr) || (pq == CGAL_EQUAL);
    }

    bool EdgeCompare(int e1, int e2) const;
    // computes the order of two edges e1 and e2 on the current sweepline

    bool EdgeCompareShared(int e1, int e2) const;
    // computes the order of two edges e1 and e2 that share a vertex

    bool EdgeCompareNonShared(int e1, int e2) const;
    // computes the order of two edges e1 and e2 that do not share a vertex

    bool EdgesShareVertex(int e1, int e2) const;
    // true if the edges e1 and e2 share a vertex

    class EdgeComp {
      private:
        const CGAL_Simplicity_test_2<ForwardIterator, Traits>* s;
      public:
        EdgeComp() {}
        EdgeComp(
          const CGAL_Simplicity_test_2<ForwardIterator, Traits>* s0): s(s0)
        {}
        bool operator() (int i, int j) const { return s->EdgeCompare(i,j); }
    };

    class EventQueue {
      //-----------------------------------------------------------------//
      // g++ 2.7.2 seems to have problems with the following typedef
      //
      // typedef set<int,VertexComp>::const_iterator const_iterator;
      //-----------------------------------------------------------------//

      private:
        std::set<int,VertexComp> queue;
      public:
        EventQueue(CGAL_Simplicity_test_2<ForwardIterator, Traits>* s)
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

      private:
        std::set<int,EdgeComp> status;
        // if i is an element of status, it means that 

        std::vector<typename std::set<int,EdgeComp>::const_iterator> index;
        // the iterators of the edges are stored to enable fast deletion

        const CGAL_Simplicity_test_2<ForwardIterator, Traits>* s;
        // store a pointer to the CGAL_Simplicity_test_2 class, to enable
        // access to the vertices

      public:
        SweepStatus(
          const CGAL_Simplicity_test_2<ForwardIterator, Traits>* s0, int n)
          : status(EdgeComp(s0)), s(s0)
        {
          index.reserve(n);
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

            CGAL_Comparison_result c1 =
              s->traits().compare_y(s->Vertex(v1), s->EventPoint());

            CGAL_Comparison_result c2 =
              s->traits().compare_y(s->Vertex(v2), s->EventPoint());

            if (c1 == CGAL_SMALLER && c2 == CGAL_SMALLER) return false;
            if (c1 == CGAL_LARGER && c2 == CGAL_LARGER) return false;
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
          // Ideally we would like to directly replace edge e1 with edge e2 by
          // putting *(index[e1]) = e2. However, this is not supported by STL
          // sets.
          erase(e1);
          insert(e2);
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
bool CGAL_Simplicity_test_2<ForwardIterator, Traits>::VertexCompare(
  int i, int j) const
{
  return !d_traits.lexicographically_yx_smaller_or_equal(Vertex(i), Vertex(j));
}

template <class ForwardIterator, class Traits>
inline
bool CGAL_Simplicity_test_2<ForwardIterator, Traits>::EdgeCompare(
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

  if (EdgesShareVertex(e1,e2)) {
    if (e1 < e2)
      Result = EdgeCompareShared(e1,e2);
    else
      Result = !EdgeCompareShared(e2,e1);
  }
  else {
    if (e1 < e2)
      Result = EdgeCompareNonShared(e1,e2);
    else
      Result = !EdgeCompareNonShared(e2,e1);
  }

#ifdef CGAL_POLYGON_DEBUG
{
  char c = Result ? '<' : '>';
  cout << "      edge " << e1 << " " << c << " edge " << e2 << endl;
}
#endif // CGAL_POLYGON_DEBUG

  return Result;
}

template <class ForwardIterator, class Traits>
bool CGAL_Simplicity_test_2<ForwardIterator, Traits>::EdgeCompareShared(
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
    return has_on_left_side(Vertex(e1), Vertex(e2), Vertex(f2));
  else
    return has_on_left_side(Vertex(f1), Vertex(f2), Vertex(e2));
}

template <class ForwardIterator, class Traits>
bool
CGAL_Simplicity_test_2<ForwardIterator, Traits>::EdgeCompareNonShared(
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
CGAL_Simplicity_test_2<ForwardIterator, Traits>::Test(ForwardIterator first,
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
  std::vector<ForwardIterator>::size_type i;
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
      status.erase(prev);
      status.erase(i);
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
bool CGAL_Simplicity_test_2<ForwardIterator, Traits>::EdgesDoIntersect(
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
  if (EdgesShareVertex(e1,e2))
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
bool CGAL_Simplicity_test_2<ForwardIterator, Traits>::EdgesShareVertex(
  int e1, int e2) const
{
  int n = NumberOfVertices();
  return ( abs(e2-e1) == 1 || abs(e2-e1) == n-1 );
}

//-----------------------------------------------------------------------//
//                          CGAL_is_simple_2
//-----------------------------------------------------------------------//
// uses Traits::compare_x
//      Traits::compare_y
//      Traits::cross_product_2
//      Traits::do_intersect
//      Traits::have_equal_direction
//      Traits::is_negative
//      Traits::lexicographically_yx_smaller_or_equal

template <class ForwardIterator, class Traits>
bool CGAL_is_simple_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits)
{
  CGAL_Simplicity_test_2<ForwardIterator, Traits> test(traits);
  return test.Test(first, last);
}

//-----------------------------------------------------------------------//
//                          CGAL_left_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy

template <class ForwardIterator, class Traits>
ForwardIterator CGAL_left_vertex_2(ForwardIterator first,
                                   ForwardIterator last,
                                   const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_xy Less_xy;
  return min_element(first, last, Less_xy());
}

//-----------------------------------------------------------------------//
//                          CGAL_right_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy

template <class ForwardIterator, class Traits>
ForwardIterator CGAL_right_vertex_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_xy Less_xy;
  return max_element(first, last, Less_xy());
}

//-----------------------------------------------------------------------//
//                          CGAL_top_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_yx

template <class ForwardIterator, class Traits>
ForwardIterator CGAL_top_vertex_2(ForwardIterator first,
                                  ForwardIterator last,
                                  const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_yx Less_yx;
  return max_element(first, last, Less_yx());
}

//-----------------------------------------------------------------------//
//                          CGAL_bottom_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_yx

template <class ForwardIterator, class Traits>
ForwardIterator CGAL_bottom_vertex_2(ForwardIterator first,
                                     ForwardIterator last,
                                     const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_yx Less_yx;
  return min_element(first, last, Less_yx());
}

//-----------------------------------------------------------------------//
//                          CGAL_bbox_2
//-----------------------------------------------------------------------//

template <class InputIterator>
CGAL_Bbox_2 CGAL_bbox_2(InputIterator first, InputIterator last)
{
  CGAL_polygon_precondition(first != last);
  CGAL_Bbox_2 result = (*first).bbox();

  while (++first != last)
    result = result + (*first).bbox();

  return result;
}

//-----------------------------------------------------------------------//
//                          CGAL_area_2
//-----------------------------------------------------------------------//
// uses Traits::determinant_2

template <class ForwardIterator, class FT, class Traits>
void CGAL_area_2(ForwardIterator first,
                 ForwardIterator last,
                 FT& result,
                 const Traits& traits)
{
  result = FT(0);

  // check if the polygon is empty
  if (first == last) return;

  ForwardIterator second = first; ++second;

  // check if the polygon has only one point
  if (second == last) return;

  ForwardIterator third = second;

  while (++third != last) {
    result = result + traits.determinant_2(*first, *second, *third);
    second = third;
  }

  result = result / FT(2);
}

//-----------------------------------------------------------------------//
//                          CGAL_is_convex_2
//-----------------------------------------------------------------------//
// uses Traits::lexicographically_xy_smaller
//      Traits::orientation

template <class ForwardIterator, class Traits>
bool CGAL_is_convex_2(ForwardIterator first,
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
      case CGAL_CLOCKWISE:
        HasClockwiseTriples = true;
        break;
      case CGAL_COUNTERCLOCKWISE:
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
//                          CGAL_oriented_side_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy
//      Traits::compare_x
//      Traits::compare_y
//      Traits::determinant_2
//      Traits::orientation
//      Traits::sign

template <class ForwardIterator, class Point, class Traits>
CGAL_Oriented_side CGAL_oriented_side_2(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point& point,
                                        const Traits& traits)
{
  CGAL_Oriented_side result;

  CGAL_Orientation o = CGAL_orientation_2(first, last, traits);
  CGAL_polygon_assertion(o != CGAL_COLLINEAR);

  CGAL_Bounded_side b = CGAL_bounded_side_2(first, last, point, traits);
  switch (b) {
    case CGAL_ON_BOUNDARY:
      result = CGAL_ON_ORIENTED_BOUNDARY;
      break;

    case CGAL_ON_BOUNDED_SIDE:
      result = (o == CGAL_CLOCKWISE) ?
               CGAL_ON_NEGATIVE_SIDE :
               CGAL_ON_POSITIVE_SIDE;
      break;

    case CGAL_ON_UNBOUNDED_SIDE:
      result = (o == CGAL_CLOCKWISE) ?
               CGAL_ON_POSITIVE_SIDE :
               CGAL_ON_NEGATIVE_SIDE;
      break;
  }

  return result;
}

//-----------------------------------------------------------------------//
//                          CGAL_bounded_side_2
//-----------------------------------------------------------------------//
// uses Traits::compare_x
//      Traits::compare_y
//      Traits::determinant_2
//      Traits::sign
//
// returns CGAL_ON_BOUNDED_SIDE, CGAL_ON_BOUNDARY or CGAL_ON_UNBOUNDED_SIDE

template <class ForwardIterator, class Point, class Traits>
CGAL_Bounded_side CGAL_bounded_side_2(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point& point,
                                      const Traits& traits)
{
  ForwardIterator current = first;
  if (current == last) return CGAL_ON_UNBOUNDED_SIDE;

  ForwardIterator next = current; ++next;
  if (next == last) return CGAL_ON_UNBOUNDED_SIDE;

  bool IsInside = false;
  CGAL_Comparison_result CompareCurrent = traits.compare_y(*current, point);

  do // check if the segment (current,next) intersects
     // the ray { (t,y) | t >= point.x() }
  {
    CGAL_Comparison_result CompareNext = traits.compare_y(*next, point);

    switch (CompareCurrent) {
      case CGAL_SMALLER:
        switch (CompareNext) {
          case CGAL_SMALLER:
            break;
          case CGAL_EQUAL:
            switch (traits.compare_x(point, *next)) {
              case CGAL_SMALLER: IsInside = !IsInside; break;
              case CGAL_EQUAL:   return CGAL_ON_BOUNDARY;
              case CGAL_LARGER:  break;
            }
            break;
          case CGAL_LARGER:
            if (point.x() < min((*current).x(), (*next).x())) {
              IsInside = !IsInside;
            }
            else if (point.x() <= max((*current).x(),(*next).x())) {
              switch (traits.sign(traits.determinant_2(point,
                                                       *current,
                                                       *next)))
              {
                case 0: return CGAL_ON_BOUNDARY;
                case 1: IsInside = !IsInside; break;
              }
            }
            break;
        }
        break;
      case CGAL_EQUAL:
        switch (CompareNext) {
          case CGAL_SMALLER:
            switch (traits.compare_x(point, *current)) {
              case CGAL_SMALLER: IsInside = !IsInside; break;
              case CGAL_EQUAL:   return CGAL_ON_BOUNDARY;
              case CGAL_LARGER:  break;
            }
            break;
          case CGAL_EQUAL:
            if ( (min((*current).x(), (*next).x()) <= point.x()) &&
                 (point.x() <= max((*current).x(), (*next).x()))    ) {
              return CGAL_ON_BOUNDARY;
            }
            break;
          case CGAL_LARGER:
            if (point.x() == (*current).x()) {
              return CGAL_ON_BOUNDARY;
            }
            break;
        }
        break;
      case CGAL_LARGER:
        switch (CompareNext) {
          case CGAL_SMALLER:
            if (point.x() < min((*current).x(), (*next).x())) {
              IsInside = !IsInside;
            }
            else if (point.x() <= max((*current).x(),(*next).x())) {
              switch (traits.sign(traits.determinant_2(point,
                                                       *current,
                                                       *next)))
              {
                case -1: IsInside = !IsInside; break;
                case  0: return CGAL_ON_BOUNDARY;
              }
            }
            break;
          case CGAL_EQUAL:
            if (point.x() == (*next).x()) {
              return CGAL_ON_BOUNDARY;
            }
            break;
          case CGAL_LARGER:
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

  return IsInside ? CGAL_ON_BOUNDED_SIDE : CGAL_ON_UNBOUNDED_SIDE;
}

//-----------------------------------------------------------------------//
//                          CGAL_orientation_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy
//      Traits::orientation

template <class ForwardIterator, class Traits>
CGAL_Orientation CGAL_orientation_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits& traits)
{
  CGAL_polygon_precondition(CGAL_is_simple_2(first, last, traits));

  ForwardIterator i = CGAL_left_vertex_2(first, last, traits);

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

