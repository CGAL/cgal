// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Simplicity_sweep_2.h
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// chapter       : Geometric Object Generators
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// implementation: Simplicity Sweep for generating random simple polygons
// ============================================================================
#ifndef CGAL_SIMPLICITY_SWEEP_2_H
#define CGAL_SIMPLICITY_SWEEP_2_H

#include <algorithm>
#include <set>
#include <vector>
#include <assert.h>

namespace CGAL {
template <class ForwardIterator, class Traits> class Simplicity_test_2;
}

#include <CGAL/polygon_assertions.h>
#include <CGAL/Simplicity_test.h>
#include <CGAL/Random.h>

namespace CGAL {

//-----------------------------------------------------------------------//
//                          Simplicity_sweep_2
//-----------------------------------------------------------------------//

template <class ForwardIterator, class Traits>
class Simplicity_sweep_2 : public Simplicity_test_2<ForwardIterator, Traits>
{
    typedef Simplicity_test_2<ForwardIterator, Traits>  Simplicity_test;
  public:
    Simplicity_sweep_2(const Traits& tr): Simplicity_test(tr) {}

    bool sweep(ForwardIterator first, ForwardIterator last, 
               std::pair<ForwardIterator, ForwardIterator>& reverse_range);
    // tests if the polygon with points in the range [first,last) is simple
    // and, if not, returns a pair of iterators that correspond to endpoints 
    // of a randomly chosen pair of intersecting edges.  The vertices spanned
    // by this range of iterators should be reversed to eliminate the
    // intersection.
};



// returns true if a change was made during the sweep (so the polygon
// was not simple) and false otherwise (indicating a simple polygon)
template <class ForwardIterator, class Traits>
bool
Simplicity_sweep_2<ForwardIterator, Traits>::sweep(
                    ForwardIterator first,
                    ForwardIterator last,
                    std::pair<ForwardIterator, ForwardIterator>& reverse_range)
{
  int n = 0;

  std::vector< std::pair<int, int> > intersections;

  typename Simplicity_test::EventQueue events(this);
  
  d_index.clear();
  while (first != last) {
    d_index.push_back(first);
    events.insert(n++); 
    ++first;
  }


  if (n < 3)
    return false;

  typename Simplicity_test::SweepStatus status(this,n);

  while (!events.empty()) 
  {
    int i = events.pop();
    d_eventpoint = i;

    int prev = (i>0) ? i-1 : i-1+n;
    int next = (i<n-1) ? i+1 : i+1-n;

    bool prev_less_than_i = VertexCompare(i,prev);
    bool next_less_than_i = VertexCompare(i,next);
    if (prev_less_than_i != next_less_than_i) 
    {
      int e = prev_less_than_i ? status.replace(i,prev) :
        status.replace(prev,i);
      CGAL_polygon_assertion(status.is_valid());

      // check for intersections of newly inserted edge e with neighbors
      int left = status.left(e);
      if ((left >= 0) && (EdgesDoIntersect(left,e))) 
      {  
         assert(left != e);
         intersections.push_back(std::pair<int, int>(left,e));
      }

      int right = status.right(e);
      if ((right >= 0) && (EdgesDoIntersect(e,right))) 
      {  
         assert(e != right);
         intersections.push_back(std::pair<int, int>(e,right));
      }
    }
    else if (prev_less_than_i) 
    {
      int e1 = prev;
      int e2 = i;

      status.insert(e1);
      status.insert(e2);
      CGAL_polygon_assertion(status.is_valid());

      // check for intersections of newly inserted edges e1 and e2 with
      // neighbors
      int left, right;
      left = status.left(e1);
      if ((left >= 0) && (EdgesDoIntersect(left,e1))) 
      {
         assert(left != e1);
         intersections.push_back(std::pair<int, int>(left,e1));
      }

      right = status.right(e1);
      if ((right >= 0) && (EdgesDoIntersect(e1,right))) 
      {
         assert(e1 != right);
         intersections.push_back(std::pair<int, int>(e1,right));
      }

      left = status.left(e2);
      if ((left >= 0) && (left != e1) && (EdgesDoIntersect(left,e2)))
      {
         assert(left != e2);
         intersections.push_back(std::pair<int, int>(left,e2));
      }

      right = status.right(e2);
      if ((right >= 0) && (right != e1) && (EdgesDoIntersect(e2,right)))
      {
         assert(e2 != right);
         intersections.push_back(std::pair<int, int>(e2,right));
      }
    }
    else 
    {
      if (status.left(prev) == i)
      {
         int left = status.left(i);
         int right = status.right(prev);
         status.erase(prev);
         status.erase(i);
         if (left >=0 && right >=0 && EdgesDoIntersect(left, right))
         {
            assert(left != right);
            intersections.push_back(std::pair<int, int>(left, right));
         }
      }
      else if (status.left(i) == prev)
      {
         int left = status.left(prev);
         int right = status.right(i);
         status.erase(prev);
         status.erase(i);
         if (left >=0 && right >=0 && EdgesDoIntersect(left, right))
         {
            assert(left != right);
            intersections.push_back(std::pair<int, int>(left, right));
         }
      }
      else
      {
         int left1 = status.left(i);
         int right1 = status.right(i);
         status.erase(i);
         if (left1 >=0 && right1 >=0 && EdgesDoIntersect(left1, right1))
         {
            assert(left1 != right1);
            intersections.push_back(std::pair<int, int>(left1, right1));
         }
         int left2 = status.left(prev);
         int right2 = status.right(prev);
         status.erase(prev);
         if (left2 >=0 && right2 >=0 && EdgesDoIntersect(left2, right2))
         {
            assert(left2 != right2);
            intersections.push_back(std::pair<int, int>(left2, right2));
         }
      }

      CGAL_polygon_assertion(status.is_valid());
    }
  }
  if (intersections.empty())
  {
      return false;
  }
  else
  {
     int i = CGAL::Random().get_int(0,intersections.size());
     if (intersections[i].first < intersections[i].second)
     {
        reverse_range.first = d_index[intersections[i].first];
        reverse_range.first++; // must be possible since second is larger
        reverse_range.second = d_index[intersections[i].second];
     }
     else
     {
        reverse_range.first = d_index[intersections[i].second];
        reverse_range.first++; // must be possible since first is larger
        reverse_range.second = d_index[intersections[i].first];
     }
     return true;
  }
}

}

#endif // CGAL_SIMPLICITY_SWEEP_2_H
