// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : pierce_rectangles_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Piercing Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_PIERCE_RECTANGLES_2_H)
#define CGAL_PIERCE_RECTANGLES_2_H 1

#include <CGAL/optimisation_assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/function_objects.h>
#include <CGAL/Transform_iterator.h>
#include <algo.h>
#include <vector.h>
#ifdef CGAL_PCENTER_WINDOW_TRACE
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>
#endif // CGAL_PCENTER_WINDOW_TRACE

//!!! to function_objects.h
template < class T1, class T2 >
struct CGAL_Has_on_unbounded_side
: public binary_function< T1, T2, bool >
{
  bool
  operator()( const T1& a, const T2& b) const
  { return a.has_on_unbounded_side( b); }
};
template < class T1, class T2 >
struct CGAL_Has_on_bounded_side
: public binary_function< T1, T2, bool >
{
  bool
  operator()( const T1& a, const T2& b) const
  { return a.has_on_bounded_side( b); }
};
template < class T1, class T2 >
struct CGAL_Has_on_boundary
: public binary_function< T1, T2, bool >
{
  bool
  operator()( const T1& a, const T2& b) const
  { return a.has_on_boundary( b); }
};

//!!! STL-extensions
template < class T >
struct CGAL_Wastebasket : public output_iterator
{
  typedef CGAL_Wastebasket< T > iterator;

  iterator
  operator=( const T&)
  { return *this; }

  iterator
  operator*()
  { return *this; }

  iterator
  operator++()
  { return *this; }

  iterator
  operator++( int)
  { return *this; }
};
template < class ForwardIterator, class OutputIterator, class Predicate >
OutputIterator
CGAL_remove_copy_if_adjacent( ForwardIterator first,
                              ForwardIterator last,
                              OutputIterator result,
                              Predicate pred)
// copy [first, last) to result such that for all x in result
// pred( *x, *++x) is not true by dropping the second of any
// two adjacent items for which pred( *x, *++x) is true
{
  if ( first == last)
    return result;
  ForwardIterator next = first;
  while ( ++next != last)
    if ( !pred( *first, *next)) {
      *result++ = *first;
      first = next;
    }
  *result++ = *first;
  return result;
}

template < class ForwardIterator, class BinaryPredicate >
ForwardIterator
CGAL_remove_if_adjacent( ForwardIterator first,
                         ForwardIterator last,
                         BinaryPredicate pred)
// transform [first, last) such that for all x
// pred( *x, *++x) is not true by dropping the second of
// any two adjacent items for which pred( *x, *++x) is true
{
  first = adjacent_find( first, last, pred);
  if ( first == last)
    return last;
  ForwardIterator save_first = first++;
  // documented in CD2:
  CGAL_optimisation_assertion( first != last);
  ForwardIterator next =
    find_if(
      ++first,
      last,
      compose1( logical_not< bool >(),
                bind1st( pred, *save_first)));
  return CGAL_remove_copy_if_adjacent( next, last, ++save_first, pred);
}
template < class ForwardIterator >
pair< ForwardIterator, ForwardIterator >
CGAL_min_max_element( ForwardIterator first,
                      ForwardIterator last)
{
  typedef pair< ForwardIterator, ForwardIterator > FP;
  if ( first == last)
    return FP( first, first);
  FP result( first, first);
  while ( ++first != last) {
    if ( *first < *(result.first))
      result.first = first;
    if ( *first > *(result.second))
      result.second = first;
  }
  return result;
}

template < class ForwardIterator, class CompareMin, class CompareMax >
pair< ForwardIterator, ForwardIterator >
CGAL_min_max_element( ForwardIterator first,
                      ForwardIterator last,
                      CompareMin comp_min,
                      CompareMax comp_max)
{
  typedef pair< ForwardIterator, ForwardIterator > FP;
  if ( first == last)
    return FP( first, first);
  FP result( first, first);
  while ( ++first != last) {
    if ( comp_min( *first, *(result.first)))
      result.first = first;
    if ( comp_max( *first, *(result.second)))
      result.second = first;
  }
  return result;
}

template < class _Traits, class _RandomAccessIC >
class CGAL__Loc_domain {
public:
  // ---------------------------------------------
  // types:

  typedef _Traits                            Traits;
  typedef _RandomAccessIC                    RandomAccessIC;
  typedef typename _Traits::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename _Traits::Point_2          Point_2;
  typedef typename _Traits::FT               FT;
  typedef typename _Traits::Xmin             Xmin;
  typedef typename _Traits::Xmax             Xmax;
  typedef typename _Traits::Ymin             Ymin;
  typedef typename _Traits::Ymax             Ymax;

  // ---------------------------------------------
  // creation:

  CGAL__Loc_domain( RandomAccessIC b, RandomAccessIC e)
  {
    CGAL_optimisation_precondition( b != e);
  
    RandomAccessIC i( b);
    r[0] = r[1] = r[2] = r[3] = i;
    while ( ++i != e) {
      check_and_update_border( 0, i);
      check_and_update_border( 1, i);
      check_and_update_border( 2, i);
      check_and_update_border( 3, i);
    }
  
    CGAL_optimisation_postcondition(
      !(*r[0]).has_on_bounded_side( vertex( 0)) &&
      !(*r[1]).has_on_bounded_side( vertex( 0)));
    CGAL_optimisation_postcondition(
      !(*r[1]).has_on_bounded_side( vertex( 1)) &&
      !(*r[2]).has_on_bounded_side( vertex( 1)));
    CGAL_optimisation_postcondition(
      !(*r[2]).has_on_bounded_side( vertex( 2)) &&
      !(*r[3]).has_on_bounded_side( vertex( 2)));
    CGAL_optimisation_postcondition(
      !(*r[3]).has_on_bounded_side( vertex( 3)) &&
      !(*r[0]).has_on_bounded_side( vertex( 3)));
  }

  // ---------------------------------------------
  // access operations:

  RandomAccessIC
  operator[]( int i) const
  // return defining rectangle (0 <-> left, 1 <-> bottom)
  {
    CGAL_optimisation_precondition( i >= 0 && i < 4);
    return r[i];
  }

  FT
  xmin() const
  // returns x-coordinate of left border
  { return Xmax()(*r[0]); }

  FT
  ymin() const
  // returns y-coordinate of lower border
  { return Ymax()(*r[1]); }

  FT
  xmax() const
  // returns x-coordinate of right border
  { return Xmin()(*r[2]); }

  FT
  ymax() const
  // returns y-coordinate of upper border
  { return Ymin()(*r[3]); }

  Point_2
  vertex( int i) const
  // PRE: 0 <= i < 4
  // POST: returns one of the four vertices:
  //  0 <=> lower left
  //  1 <=> lower right
  //  2 <=> upper right
  //  3 <=> upper left
  {
    CGAL_optimisation_precondition( i >= 0 && i < 4);
    typedef typename _Traits::Build_point Build;
    switch( i) {
    case 0:
      return Build()( xmin(), ymin());
    case 1:
      return Build()( xmax(), ymin());
    case 2:
      return Build()( xmax(), ymax());
    // case 3:
    //  return Point_2( xmin(), ymax());
    }
    return Build()( xmin(), ymax());
  }

  Point_2
  min() const
  // return lexicographically smallest vertex
  // (in analogy to CGAL_Iso_rectangle)
  { return vertex( 0); }

  Point_2
  max() const
  // return lexicographically largest vertex
  // (in analogy to CGAL_Iso_rectangle)
  { return vertex( 2); }

  // ---------------------------------------------
  // update operations:
  
  void
  init( int i, RandomAccessIC f)
  // PRE: 0 <= i < 4
  // POST: set defining rectangle r[i] to f
  {
    CGAL_optimisation_precondition( i >= 0 && i < 4);
    r[i] = f;
  }
  
  void
  check_and_update_border( int i, RandomAccessIC a)
  // PRE: 0 <= i < 4
  // POST: check whether rectangle *a invalidates
  //   border i. if so, update defining rectangle r[i]
  //   to a.
  {
    CGAL_optimisation_precondition( i >= 0 && i < 4);
    switch( i) {
    case 0:
      if ( Xmax()(*a) < Xmax()(*r[0]))
        r[0] = a;
      break;
    case 1:
      if ( Ymax()(*a) < Ymax()(*r[1]))
        r[1] = a;
      break;
    case 2:
      if ( Xmin()(*a) > Xmin()(*r[2]))
        r[2] = a;
      break;
    case 3:
      if ( Ymin()(*a) > Ymin()(*r[3]))
        r[3] = a;
      break;
    }
  }

private:
  // pointer to defining rectangles:
  RandomAccessIC r[4];
}; // class CGAL__Loc_domain
template < class _Traits, class _RandomAccessIC >
class CGAL__Rectangle_partition {
public:
  typedef _Traits                            Traits;
  typedef _RandomAccessIC                    RandomAccessIC;
  typedef typename _Traits::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename _Traits::FT               FT;
  typedef typename _Traits::Xmin             Xmin;
  typedef typename _Traits::Xmax             Xmax;
  typedef typename _Traits::Ymin             Ymin;
  typedef typename _Traits::Ymax             Ymax;
  typedef typename _Traits::Build_point      Build_point;
  typedef typename _Traits::Build_rectangle  Build_rectangle;

  // constants to identify the eleven partition sets:
  // BL and TL are sorted according to right sides
  //   in increasing order,
  // BR, TR and BT are sorted according to left sides
  //   in decreasing order,
  // LR is sorted according to top sides
  //   in increasing order and
  // TR, TL, L, R, NO are not sorted
  enum set_id { BL, TL, BT, BR, TR, LR, B, T, L, R, NO};

  typedef vector< Iso_rectangle_2 >             Container;
  typedef typename Container::iterator          iterator;
  typedef typename Container::reverse_iterator  reverse_iterator;
  typedef typename Container::const_iterator    const_iterator;

  // constructor: partition the range [f,l) of rectangles
  // according to d
  // precondition:
  //   d is the location domain associated with [f, l)
  CGAL__Rectangle_partition(
    RandomAccessIC f,
    RandomAccessIC l,
    const CGAL__Loc_domain< _Traits, RandomAccessIC >& d);

  // ---------------------------------------------------
  // access functions to the partition sets:
  
  const_iterator
  begin( set_id i) const
  {
    CGAL_optimisation_precondition( i >= BL && i <= NO);
    const_iterator p( s[i].begin());
    if ( i <= LR)
      return ++p;
    else
      return p;
  }
  
  const_iterator
  end( set_id i) const
  {
    CGAL_optimisation_precondition( i >= BL && i <= NO);
    const_iterator p( s[i].end());
    if ( i <= LR)
      return --p;
    else
      return p;
  }
  
  bool
  is_empty( set_id i) const
  {
    CGAL_optimisation_precondition( i >= BL && i <= NO);
    CGAL_optimisation_assertion( begin( i) <= end( i));
    return begin( i) == end( i);
  }
  const_iterator
  first_right_of( set_id i, FT v, const_iterator p) const
  // PRE: i == BL or i == BT, p in [ begin(i), end(i) ] and
  //   i == BT ==>
  //     *p lies to the left of (not including) the line x = v.
  // i == BL:
  //   returns the maximal q from [ p, end(i) ]
  //   such that all rectangles from [ p, q )
  //   lie to the left of (not including) the line x = v.
  // i == BT:
  //   returns the minimal q from [ begin(i), p ]
  //   such that all rectangles from [ q, end(i) )
  //   lie to the left of (not including) the line x = v.
  // TIME: O( #elements in set i between p and q)
  {
    CGAL_optimisation_precondition(
      (i == BL || i == BT && xmax( *p) < v) &&
      p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition_code( const_iterator pp( p);)
    // s[BL] is sorted acc. to right sides incr.
    if ( i == BL) {
      while ( xmax( *p) < v)
        ++p;
      CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
      CGAL_optimisation_postcondition( p == pp || xmax( *(p - 1)) < v);
    } else {
      // s[BT] is sorted acc. to left sides decr.
      while ( xmax( *--p) < v) {}
      ++p;
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition( p == end(i) || xmax( *p) < v);
    }
    return p;
  }
  
  const_iterator
  first_right_of( set_id i, FT v) const
  // PRE: i == BL or i == BT
  // i == BL:
  //   returns the maximal q from [ begin(i), end(i) ]
  //   such that all rectangles from [ begin(i), q )
  //   lie to the left of (not including) the line x = v.
  // i == BT:
  //   returns the minimal q from [ begin(i), end(i) ]
  //   such that all rectangles from [ q, end(i) )
  //   lie to the left of (not including) the line x = v.
  // TIME: O( log(#elements in set i))
  {
    CGAL_optimisation_precondition( i == BL || i == BT);
    const_iterator p;
    if ( i == BL) {
      p = lower_bound( begin( BL),
                       s[BL].end(),
                       v,
                       CGAL_compose2_2( less< FT >(),
                                        xmax,
                                        identity< FT >()));
      CGAL_optimisation_postcondition(
        p == begin(i) || xmax( *(p - 1)) < v);
    } else {
      p = lower_bound( begin( i),
                       s[i].end(),
                       v,
                       CGAL_compose2_2( greater_equal< FT >(),
                                        xmax,
                                        identity< FT >()));
      CGAL_optimisation_postcondition( p == end(i) || xmax( *p) < v);
    }
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    return p;
  }
  
  const_iterator
  first_left_of( set_id i, FT v, const_iterator p) const
  // PRE: i == BT, p in [ begin(i), end(i) ]
  // returns the maximal q from [ begin(i), p ]
  // such that all rectangles from [ begin(i), q )
  // lie to the right of (not including) the line x = v.
  // TIME: O( #elements in set i between p and q)
  {
    CGAL_optimisation_precondition(
      i == BT && p >= begin(i) && p <= end(i));
    while ( xmin( *p) <= v)
      --p;
    ++p;
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition(
      p == begin(i) || Xmin()( *(p - 1)) > v);
    CGAL_optimisation_postcondition( xmax( *p) <= v);
    return p;
  }
  
  const_iterator
  first_left_of( set_id i, FT v) const
  // PRE: i == BR or i == BT.
  // returns the maximal q from [ begin(i), end(i) ]
  // such that all rectangles from [ begin(i), q )
  // lie to the right of (not including) the line x = v.
  // TIME: O( log(#elements in set i))
  {
    CGAL_optimisation_precondition( i == BR || i == BT);
    const_iterator p =
      lower_bound( begin( i),
                   s[i].end(),
                   v,
                   CGAL_compose2_2( greater< FT >(),
                                    Xmin(),
                                    identity< FT >()));
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition(
      p == begin(i) || Xmin()( *(p - 1)) > v);
    return p;
  }
  
  const_iterator
  first_above_downwards( set_id i, FT v, const_iterator p) const
  // PRE: i == TL or i == LR, p in [ begin(i), end(i) ] and
  //      *p lies above (not including) the line y = v.
  //   returns the minimal q from [ begin(i), p ]
  //   such that all rectangles from [ q, end(i) )
  //   lie above (not including) the line y = v.
  // TIME: O( #elements in set i between p and q)
  {
    CGAL_optimisation_precondition( i == TL || i == LR);
    CGAL_optimisation_precondition(
      ymin(*p) > v && p >= begin(i) && p <= end(i));
    while ( ymin( *--p) > v) {}
    ++p;
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition( p == end(i) || Ymin()( *p) > v);
    return p;
  }
  
  const_iterator
  first_above_upwards( set_id i, FT v, const_iterator p) const
  // PRE: i == TR or i == LR and p in [ begin(i), end(i) ]
  //   returns the minimal q from [ p, end(i) ]
  //   such that all rectangles from [ q, end(i) )
  //   lie above (not including) the line y = v.
  // TIME: O( #elements in set i between p and q)
  {
    CGAL_optimisation_precondition(
      (i == TR || i == LR) && p >= begin(i) && p <= end(i));
    while ( ymin( *p) <= v)
      ++p;
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition( p == end(i) || Ymin()( *p) > v);
    return p;
  }
  
  const_iterator
  first_above( set_id i, FT v) const
  // PRE: i == TL or i == TR or i == LR
  //   returns the minimal q from [ begin(i), end(i) ]
  //   such that all rectangles from [ q, end(i) )
  //   lie above (not including) the line y = v.
  // TIME: O( log(#elements in set i))
  {
    CGAL_optimisation_precondition( i == TL || i == TR || i == LR);
    const_iterator p =
      lower_bound( begin( i),
                   s[i].end(),
                   v,
                   CGAL_compose2_2( less_equal< FT >(),
                                    Ymin(),
                                    identity< FT >()));
    CGAL_optimisation_postcondition( p >= begin(i) && p <= end(i));
    CGAL_optimisation_postcondition( p == end(i) || Ymin()( *p) > v);
    return p;
  }

private:
  void
  sort_set( set_id i)
  {
    CGAL_optimisation_precondition( i >= BL && i <= LR);
    CGAL_optimisation_precondition( s[i].begin() + 1 == begin( i));
  
    // NB: only used in constructor.
    //   using s[...].end() instead of end(...)
    //   since upper sentinels have not been inserted yet
  
    if ( i <= TL)
      sort( s[i].begin() + 1,
            s[i].end(),
            CGAL_compose2_2( less< FT >(), Xmax(), Xmax()));
    else if ( i <= TR)
      sort( s[i].begin() + 1,
            s[i].end(),
            CGAL_compose2_2( greater< FT >(), Xmin(), Xmin()));
    else if ( i == LR)
      sort( s[i].begin() + 1,
            s[i].end(),
            CGAL_compose2_2( less< FT >(), Ymax(), Ymax()));
  } // sort_set( set_id i)
  void
  remove_containing_rectangles()
  // make sure that in any partition set there are no
  // two rectangles such that if they are clipped
  // against the location domain, one contains the other.
  {
    iterator new_end;
  
    // NB: only used in constructor.
    //   using s[...].end() instead of end(...)
    //   since upper sentinels have not been inserted yet
  
    // ------------------------------------------------------
    // sorted according to their right sides (increasing):
  
    // s[BL]: discard second, if its top side is above
    CGAL_optimisation_assertion( s[BL].begin() + 1 == begin( BL));
    new_end =
      CGAL_remove_if_adjacent(
        s[BL].begin() + 1,
        s[BL].end(),
        CGAL_compose2_2( less_equal< FT >(), Ymax(), Ymax()));
    s[BL].erase( new_end, s[BL].end());
  
    // s[TL]: discard second, if its bottom side is below
    CGAL_optimisation_assertion( s[TL].begin() + 1 == begin( TL));
    new_end =
      CGAL_remove_if_adjacent(
        s[TL].begin() + 1,
        s[TL].end(),
        CGAL_compose2_2( greater_equal< FT >(), Ymin(), Ymin()));
    s[TL].erase( new_end, s[TL].end());
  
    // ------------------------------------------------------
    // sorted according to their left sides (decreasing):
  
    // s[BT]: discard second, if its right side is right
    CGAL_optimisation_assertion( s[BT].begin() + 1 == begin( BT));
    new_end =
      CGAL_remove_if_adjacent(
        s[BT].begin() + 1,
        s[BT].end(),
        CGAL_compose2_2( less_equal< FT >(), Xmax(), Xmax()));
    s[BT].erase( new_end, s[BT].end());
  
    // s[BR]: discard second, if its top side is above
    CGAL_optimisation_assertion( s[BR].begin() + 1 == begin( BR));
    new_end =
      CGAL_remove_if_adjacent(
        s[BR].begin() + 1,
        s[BR].end(),
        CGAL_compose2_2( less_equal< FT >(), Ymax(), Ymax()));
    s[BR].erase( new_end, s[BR].end());
  
    // s[TR]: discard second, if its bottom side is below
    CGAL_optimisation_assertion( s[TR].begin() + 1 == begin( TR));
    new_end =
      CGAL_remove_if_adjacent(
        s[TR].begin() + 1,
        s[TR].end(),
        CGAL_compose2_2( greater_equal< FT >(), Ymin(), Ymin()));
    s[TR].erase( new_end, s[TR].end());
  
    // ------------------------------------------------------
    // sorted according to their top sides (increasing):
  
    // s[LR]: discard second, if its bottom side is below
    CGAL_optimisation_assertion( s[LR].begin() + 1 == begin( LR));
    new_end =
      CGAL_remove_if_adjacent(
        s[LR].begin() + 1,
        s[LR].end(),
        CGAL_compose2_2( greater_equal< FT >(), Ymin(), Ymin()));
    s[LR].erase( new_end, s[LR].end());
  
  } // remove_containing_rectangles()

  // the partition sets:
  Container s[11];

  // coordinate accessors:
  Xmin xmin;
  Xmax xmax;
  Ymin ymin;
  Ymax ymax;

  // point builder:
  Build_point  build_point;

  // rectangle builder:
  Build_rectangle  build_rectangle;
};
template < class _Traits, class _RandomAccessIC >
CGAL__Rectangle_partition< _Traits, _RandomAccessIC>::
CGAL__Rectangle_partition(
  _RandomAccessIC f,
  _RandomAccessIC l,
  const CGAL__Loc_domain< _Traits, _RandomAccessIC >& d)
{
  //!!! reserves

  // sentinel rectangles (si <-> d.vertex(i)):
  FT rad = (f == l ? 0 : (*f).radius());

  Iso_rectangle_2 s0(
    build_rectangle(
      build_point( d.xmin() - FT(1) - rad, d.ymin() - FT(1) - rad)));
  Iso_rectangle_2 s1(
    build_rectangle(
      build_point( d.xmax() + FT(1) + rad, d.ymin() - FT(1) - rad)));
  Iso_rectangle_2 s2(
    build_rectangle(
      build_point( d.xmax() + FT(1) + rad, d.ymax() + FT(1) + rad)));
  Iso_rectangle_2 s3(
    build_rectangle(
      build_point( d.xmin() - FT(1) - rad, d.ymax() + FT(1) + rad)));

  // insert lower sentinels:
  s[BL].push_back( s3);
  s[BR].push_back( s2);
  s[BT].push_back( s2);
  s[TL].push_back( s0);
  s[TR].push_back( s1);
  s[LR].push_back( s1);

  while ( f != l) {
    // otherwise d would not be the location domain:
    CGAL_optimisation_expensive_assertion(
      (*f).xmax() >= d.xmin() && (*f).ymax() >= d.ymin() &&
      (*f).xmin() <= d.xmax() && (*f).ymin() <= d.ymax());

    if ( (*f).xmin() > d.xmin())
      if ( (*f).xmax() < d.xmax())
        // *f belongs to S_b, S_t, S_bt or S_no
        if ( (*f).ymin() > d.ymin())
          if ( (*f).ymax() < d.ymax())
            s[NO].push_back( *f);
          else
            s[T].push_back( *f);
        else
          if ( (*f).ymax() < d.ymax())
            s[B].push_back( *f);
          else
            s[BT].push_back( *f);
      else
        // *f belongs to S_r, S_tr, S_br or has to be discarded
        if ( (*f).ymin() > d.ymin())
          if ( (*f).ymax() < d.ymax())
            s[R].push_back( *f);
          else
            s[TR].push_back( *f);
        else {
          if ( (*f).ymax() < d.ymax())
            s[BR].push_back( *f);
        }
    else
      if ( (*f).xmax() < d.xmax())
        // *f belongs to S_l, S_tl, S_bl or has to be discarded
        if ( (*f).ymin() > d.ymin())
          if ( (*f).ymax() < d.ymax())
            s[L].push_back( *f);
          else
            s[TL].push_back( *f);
        else {
          if ( (*f).ymax() < d.ymax())
            s[BL].push_back( *f);
        }
      else {
        // *f belongs to S_lr or has to be discarded
        if ( (*f).ymin() > d.ymin() &&
             (*f).ymax() < d.ymax())
          s[LR].push_back( *f);
      }
    ++f;
  } // while ( f != l)

  // sort partition sets:
  sort_set( BL);
  sort_set( BR);
  sort_set( BT);
  sort_set( TL);
  sort_set( TR);
  sort_set( LR);

  // remove the rectangles containing a following rectangle:
  remove_containing_rectangles();

  // insert upper sentinels:
  s[BL].push_back( s1);
  s[BR].push_back( s0);
  s[BT].push_back( s3);
  s[TL].push_back( s2);
  s[TR].push_back( s3);
  s[LR].push_back( s3);

} // Rectangle_partition( f, l, d)

#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Pierce_rectangles_2_traits.h>
#endif // CGAL_REP_CLASS_DEFINED

template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_two_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
  bool& ok,
  const Traits&)
{
  CGAL_optimisation_precondition( f != l);

  // compute location domain:
  typedef CGAL__Loc_domain< Traits, RandomAccessIC > Loc_domain;
  Loc_domain d( f, l);

  return CGAL_two_pierce_rectangles( f, l, d, o, ok);
} // CGAL_two_pierce_rectangles( f, l, o, ok, i)
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_three_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
  bool& ok,
  const Traits&)
{
  CGAL_optimisation_precondition( f != l);

  // compute location domain:
  typedef CGAL__Loc_domain< Traits, RandomAccessIC > Loc_domain;
  Loc_domain d( f, l);

  return CGAL_three_pierce_rectangles( f, l, d, o, ok);
} // CGAL_three_pierce_rectangles( f, l, o, ok, i)
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_four_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
  bool& ok,
  const Traits&)
{
  CGAL_optimisation_precondition( f != l);

  // compute location domain:
  typedef CGAL__Loc_domain< Traits, RandomAccessIC > Loc_domain;
  Loc_domain d( f, l);

  return CGAL_four_pierce_rectangles( f, l, d, o, ok);
} // CGAL_four_pierce_rectangles( f, l, o, ok, i)
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_two_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  const CGAL__Loc_domain< Traits, RandomAccessIC >& d,
  OutputIterator o,
  bool& ok)
{
  CGAL_optimisation_precondition( f != l);

  // typedefs:
  typedef typename Traits::Point_2
    Point_2;
  typedef typename Traits::Iso_rectangle_2
    Iso_rectangle_2;
  typedef CGAL_Has_on_unbounded_side< Iso_rectangle_2, Point_2 >
    Has_on_unbounded_side;

  #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 2
  cerr << " ++ 2 pierce start ++" << endl;
  for ( RandomAccessIC u = f; u != l; ++u)
    cerr << &(*u) << " : " << *u << endl;
  cerr << "LocDomain is " << d.vertex( 0) << " --> " <<
    d.vertex( 2) << endl;
  /*
  CGAL_Window_stream W;
  W.init( -1024 / 5, 1024 + 1024 / 5, -1024 / 5);
  W << CGAL_RED;
  WindowOutput( W, f, l);
  W << CGAL_GREEN;
  W << Iso_rectangle_2( d.vertex( 0), d.vertex( 2));
  double x, y;
  while ( W.read_mouse( x, y) != -1) {}
  */
  #endif

  if ( !(*d[0]).has_on_unbounded_side( d.vertex( 0)) &&
       !(*d[1]).has_on_unbounded_side( d.vertex( 0)) &&
       !(*d[2]).has_on_unbounded_side( d.vertex( 0)) &&
       !(*d[3]).has_on_unbounded_side( d.vertex( 0))) {
    // the location domain is degenerate
    // and [f,l) is one-pierceable
    *o++ = d.vertex( 0);
    ok = true;
    return o;
  }
  // check, if {d.vertex( 0), d.vertex( 2)}
  // form a piercing set
  if ( l == find_if(
    f,
    l,
    compose2( logical_and< bool >(),
              bind2nd( Has_on_unbounded_side(),
                       d.vertex( 0)),
              bind2nd( Has_on_unbounded_side(),
                       d.vertex( 2)))))
  {
    // {d.vertex( 0), d.vertex( 2)} is a piercing set
    *o++ = d.vertex( 0);
    *o++ = d.vertex( 2);
    ok = true;
    return o;
  }
  // check, if {d.vertex( 1), d.vertex( 3)}
  // form a piercing set
  if ( l == find_if(
    f,
    l,
    compose2( logical_and< bool >(),
              bind2nd( Has_on_unbounded_side(),
                       d.vertex( 1)),
              bind2nd( Has_on_unbounded_side(),
                       d.vertex( 3)))))
  {
    // {d.vertex( 1), d.vertex( 3)} is a piercing set
    *o++ = d.vertex( 1);
    *o++ = d.vertex( 3);
    ok = true;
    return o;
  }

  // no piercing set exists:
  ok = false;
  return o;
} // CGAL_two_pierce_rectangles( f, l, d, o, ok)
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_three_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  CGAL__Loc_domain< Traits, RandomAccessIC >& d,
  OutputIterator o,
  bool& ok)
{
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef ptrdiff_t difference_type;
#else  // CGAL_CFG_NO_ITERATOR_TRAITS //
  typedef typename iterator_traits< RandomAccessIC >::difference_type
    difference_type;
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
  difference_type number_of_points( CGAL_iterator_distance( f, l));
  CGAL_optimisation_precondition( number_of_points > 0);

  // typedefs:
  typedef typename Traits::Point_2          Point_2;
  typedef typename Traits::Iso_rectangle_2  Iso_rectangle_2;
  typedef CGAL_Has_on_unbounded_side< Iso_rectangle_2, Point_2 >
    Has_on_unbounded_side;

#if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 3
  cerr << " ++ 3 pierce start ++" << endl;
  for ( RandomAccessIC u = f; u != l; ++u)
    cerr << &(*u) << " : " << *u << endl;
  cerr << "LocDomain is " << d.vertex( 0) << " --> "
       << d.vertex( 2) << endl;
#endif

  // define container to store disjoint rectangles:
  // typedef vector< Iso_rectangle_2 >  Rectangle_cont;
  typedef vector< Iso_rectangle_2 >   Rectangle_cont;
  typedef Rectangle_cont::size_type   size_type;
  typedef Rectangle_cont::iterator    Rectangle_iterator;
  Rectangle_cont disjoint;
  disjoint.reserve( number_of_points);
  
  // test the four corners:
  for ( int k( 0); k < 4; ++k) {
    disjoint.erase( disjoint.begin(), disjoint.end());
    #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 3
    cerr << "check corner nr. " << k << " ("
         << d.vertex(k) << ")" << endl;
    #endif
    
    // save changing sides of d:
    RandomAccessIC save_side1( d[k]);
    RandomAccessIC save_side2( d[( k + 1) % 4]);
    
    // extract all rectangles not containing d.vertex( k)
    // into rc and update d:
    typedef CGAL_Has_on_unbounded_side< Iso_rectangle_2, Point_2 >
      Does_not_contain;
    
    /*
    typedef CGAL_Get_address< Iso_rectangle_2 > Rect_address;
    typedef CGAL_Transform_iterator<
      back_insert_iterator< Rectangle_cont >,
      Rect_address >
    Address_iterator;
    */
    
    // find first rectangle not containing d.vertex( k):
    RandomAccessIC f_not(
      find_if( f, l, bind2nd( Does_not_contain(), d.vertex( k))));
    
    // do all rectangles contain the corner?
    if ( f_not == l) {
      CGAL_optimisation_assertion( k == 0);
      *o++ = d.vertex( k);
      ok = true;
      return o;
    } // if ( f_not == l)
    
    // (*f_not) does not contain the corner, so it can be used
    // to initialize the two (possibly) changing sides of d:
    Point_2 test_vertex( d.vertex( k));
    d.init( k, f_not);
    d.init( ( k + 1) % 4, f_not);
    disjoint.push_back( *f_not);
    
    // now run through it:
    while ( ++f_not != l)
      if ( Does_not_contain()( *f_not, test_vertex)) {
        disjoint.push_back( *f_not);
        d.check_and_update_border( k, f_not);
        d.check_and_update_border( ( k + 1) % 4, f_not);
      } // if ( Does_not_contain()( *f_not, d.vertex( k)))
    
    // check disjoint for two-pierceability:
    CGAL_two_pierce_rectangles(
      disjoint.begin(),
      disjoint.end(),
      d,
      o,
      ok);
    
    // restore saved sides of d:
    d.init( k, save_side1);
    d.init( ( k + 1) % 4, save_side2);
    
    if ( ok) {
      // does any rectangle contain the corner?
      if ( disjoint.size() < CGAL_static_cast( size_type, number_of_points))
        *o++ = d.vertex( k);
      return o;
    } // if ( ok)
  } // for ( int k( 0); k < 4; ++k)

  // no piercing set exists:
  ok = false;
  return o;

} // CGAL_three_pierce_rectangles( f, l, d, o, ok)
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
inline
OutputIterator
CGAL_four_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  CGAL__Loc_domain< Traits, RandomAccessIC >& d,
  OutputIterator o,
  bool& ok)
{
  // construct partition:
  typedef CGAL__Rectangle_partition< Traits, RandomAccessIC > Partition;
  Partition p( f, l, d);

  return CGAL_four_pierce_rectangles( f, l, d, p, o, ok);
} // CGAL_four_pierce_rectangles( f, l, d, o, ok)
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_four_pierce_rectangles(
  RandomAccessIC f,
  RandomAccessIC l,
  CGAL__Loc_domain< Traits, RandomAccessIC >& d,
  CGAL__Rectangle_partition< Traits, RandomAccessIC >& p,
  OutputIterator o,
  bool& ok)
{
  int number_of_points( CGAL_iterator_distance( f, l));
  CGAL_optimisation_precondition( number_of_points > 0);

#if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
  cerr << " ** 4 pierce start **" << endl;
  for ( RandomAccessIC u = f; u != l; ++u)
    cerr << &(*u) << " : " << *u << endl;
  cerr << "LocDomain is " << d.vertex( 0) << " --> "
       << d.vertex( 2) << endl;
#endif

  typedef typename Traits::Iso_rectangle_2     Iso_rectangle_2;
  typedef typename Traits::Point_2             Point_2;
  typedef typename Traits::FT                  FT;
  typedef typename Traits::Xmin                Xmin;
  typedef typename Traits::Xmax                Xmax;
  typedef typename Traits::Ymin                Ymin;
  typedef typename Traits::Ymax                Ymax;
  typedef typename Traits::Build_point         Build_point;
  typedef typename Traits::Build_rectangle     Build_rectangle;
  typedef pair< FT, FT >                       Intervall;
  typedef CGAL__Rectangle_partition< Traits, RandomAccessIC >  RP;
  typedef typename RP::const_iterator          iterator;
  typedef pair< iterator, iterator >           iterator_pair;
  
  // point builder:
  Build_point build_point;
  
  // coordinate accessors:
  Xmin xmin;
  Xmax xmax;
  Ymin ymin;
  Ymax ymax;
  // define container to store disjoint rectangles:
  // typedef vector< Iso_rectangle_2 >  Rectangle_cont;
  typedef vector< Iso_rectangle_2 >   Rectangle_cont;
  typedef Rectangle_cont::size_type   size_type;
  typedef Rectangle_cont::iterator    Rectangle_iterator;
  Rectangle_cont disjoint;
  disjoint.reserve( number_of_points);
  
  // test the four corners:
  for ( int k( 0); k < 4; ++k) {
    disjoint.erase( disjoint.begin(), disjoint.end());
    #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 3
    cerr << "check corner nr. " << k << " ("
         << d.vertex(k) << ")" << endl;
    #endif
    
    // save changing sides of d:
    RandomAccessIC save_side1( d[k]);
    RandomAccessIC save_side2( d[( k + 1) % 4]);
    
    // extract all rectangles not containing d.vertex( k)
    // into rc and update d:
    typedef CGAL_Has_on_unbounded_side< Iso_rectangle_2, Point_2 >
      Does_not_contain;
    
    /*
    typedef CGAL_Get_address< Iso_rectangle_2 > Rect_address;
    typedef CGAL_Transform_iterator<
      back_insert_iterator< Rectangle_cont >,
      Rect_address >
    Address_iterator;
    */
    
    // find first rectangle not containing d.vertex( k):
    RandomAccessIC f_not(
      find_if( f, l, bind2nd( Does_not_contain(), d.vertex( k))));
    
    // do all rectangles contain the corner?
    if ( f_not == l) {
      CGAL_optimisation_assertion( k == 0);
      *o++ = d.vertex( k);
      ok = true;
      return o;
    } // if ( f_not == l)
    
    // (*f_not) does not contain the corner, so it can be used
    // to initialize the two (possibly) changing sides of d:
    Point_2 test_vertex( d.vertex( k));
    d.init( k, f_not);
    d.init( ( k + 1) % 4, f_not);
    disjoint.push_back( *f_not);
    
    // now run through it:
    while ( ++f_not != l)
      if ( Does_not_contain()( *f_not, test_vertex)) {
        disjoint.push_back( *f_not);
        d.check_and_update_border( k, f_not);
        d.check_and_update_border( ( k + 1) % 4, f_not);
      } // if ( Does_not_contain()( *f_not, d.vertex( k)))
    
    // check disjoint for two-pierceability:
    CGAL_three_pierce_rectangles(
      disjoint.begin(),
      disjoint.end(),
      d,
      o,
      ok);
    
    // restore saved sides of d:
    d.init( k, save_side1);
    d.init( ( k + 1) % 4, save_side2);
    
    if ( ok) {
      // does any rectangle contain the corner?
      if ( disjoint.size() < CGAL_static_cast( size_type, number_of_points))
        *o++ = d.vertex( k);
      return o;
    } // if ( ok)
  } // for ( int k( 0); k < 4; ++k)
  #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
  cerr << "traversing boundary ..." << endl;
  #endif
  
  // test if four p'points can be placed on the boundary of d,
  // one on each side
  
  // if there is any rectangle completely contained in d,
  // there is no hope for a four-piercing in this way
  if ( p.is_empty( RP::NO)) {
  
  #if defined(CGAL_PCENTER_WINDOW_TRACE)
    // graphic debug window:
    CGAL_Window_stream Wd;
    const int XYSIZE( 1024);
    double dummy_x, dummy_y;
    CGAL_Ostream_iterator< Point_2, CGAL_Window_stream > Witer;
    Witer wout;
  #endif
  
    #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
    cerr << "compute one-piercing intervalls ..." << endl;
    #endif
    
    iterator_pair ip;
    // compute the intervall I_B where s[B] is (potentially)
    // one-pierced:
    // I_B is defined by the smallest maxpoint and the
    // largest minpoint of the rectangles in s[B]
    
    Intervall I_B;
    
    if ( p.is_empty( RP::B))
      I_B = Intervall( d.xmin(), d.xmax());
    else {
      ip = CGAL_min_max_element( p.begin( RP::B),
                                 p.end( RP::B),
                                 CGAL_compose2_2(
                                   less< FT >(), Xmax(), Xmax()),
                                 CGAL_compose2_2(
                                   greater< FT >(), Xmin(), Xmin()));
    
      I_B = Intervall( Xmin()(*(ip.second)), Xmax()(*(ip.first)));
      CGAL_optimisation_assertion(
        I_B.first >= d.xmin() && I_B.second <= d.xmax());
    }
    
    
    // compute the intervall I_L where s[L] is (potentially)
    // one-pierced:
    // I_L is defined by the smallest maxpoint and the
    // largest minpoint of the rectangles in s[L]
    
    Intervall I_L;
    
    if ( p.is_empty( RP::L))
      I_L = Intervall( d.ymin(), d.ymax());
    else {
      ip = CGAL_min_max_element( p.begin( RP::L),
                                 p.end( RP::L),
                                 CGAL_compose2_2(
                                   less< FT >(), Ymax(), Ymax()),
                                 CGAL_compose2_2(
                                   greater< FT >(), Ymin(), Ymin()));
    
      I_L = Intervall( Ymin()(*(ip.second)), Ymax()(*(ip.first)));
      CGAL_optimisation_assertion(
        I_L.first >= d.ymin() && I_L.second <= d.ymax());
    }
    
    
    // compute the intervall I_T where s[T] is (potentially)
    // one-pierced:
    // I_T is defined by the smallest maxpoint and the
    // largest minpoint of the rectangles in s[T]
    
    Intervall I_T;
    
    if ( p.is_empty( RP::T))
      I_T = Intervall( d.xmin(), d.xmax());
    else {
      ip = CGAL_min_max_element( p.begin( RP::T),
                                 p.end( RP::T),
                                 CGAL_compose2_2(
                                   less< FT >(), Xmax(), Xmax()),
                                 CGAL_compose2_2(
                                   greater< FT >(), Xmin(), Xmin()));
    
      I_T = Intervall( Xmin()(*(ip.second)), Xmax()(*(ip.first)));
      CGAL_optimisation_assertion(
        I_T.first >= d.xmin() && I_T.second <= d.xmax());
    }
    
    
    // compute the intervall I_R where s[R] is (potentially)
    // one-pierced:
    // I_R is defined by the smallest maxpoint and the
    // largest minpoint of the rectangles in s[R]
    
    Intervall I_R;
    
    if ( p.is_empty( RP::R))
      I_R = Intervall( d.ymin(), d.ymax());
    else {
      ip = CGAL_min_max_element( p.begin( RP::R),
                                 p.end( RP::R),
                                 CGAL_compose2_2(
                                   less< FT >(), Ymax(), Ymax()),
                                 CGAL_compose2_2(
                                   greater< FT >(), Ymin(), Ymin()));
    
      I_R = Intervall( Ymin()(*(ip.second)), Ymax()(*(ip.first)));
      CGAL_optimisation_assertion(
        I_R.first >= d.ymin() && I_R.second <= d.ymax());
    }
    
    
    
    #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
    cerr << "one piercing intervalls pierce ..." << endl;
    #endif
    
    #if defined(CGAL_PCENTER_WINDOW_TRACE)
    Wd.init( -XYSIZE / 5, XYSIZE + XYSIZE / 5, -XYSIZE / 5);
    Wd << CGAL_RED;
    copy( f, l, wout);
    Wd << CGAL_GREEN
       << CGAL_Iso_rectangleC2<FT>( d.min(), d.max())
       << CGAL_ORANGE
       << build_point( I_B.first, d.ymin())
       << build_point( I_B.second, d.ymin())
       << build_point( d.xmin(), I_L.first)
       << build_point( d.xmin(), I_L.second)
       << build_point( I_T.first, d.ymax())
       << build_point( I_T.second, d.ymax())
       << build_point( d.xmax(), I_R.first)
       << build_point( d.xmax(), I_R.second);
    while ( Wd.read_mouse( dummy_x, dummy_y) != -1) {}
    #endif
  
    // now try to position the bottom piercing point in each
    // of the intervalls formed by S_bt and S_br
    // (no need to consider S_bl, since we move from left
    // to right and leaving a rectangle won't make piercing easier)
  
    #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
    cerr << "round trip ..." << endl;
    #endif
    
    // x coordinate of bottom piercing point:
    FT bot( I_B.first);
    
    // first rectangle in s[BL] pierced by bot
    // ( ==> predecessor gives an upper bound for lef):
    iterator first_pierced_in_BL( p.first_right_of( RP::BL, bot));
    
    // upper bound for y coordinate of left piercing point:
    FT lef( min( ymax(*(first_pierced_in_BL - 1)), I_L.second));
    
    // first rectangle in s[BR] pierced by bot
    // ( ==> predecessor gives an upper bound for rig):
    iterator first_pierced_in_BR( p.first_left_of( RP::BR, bot));
    
    // upper bound for y coordinate of right piercing point:
    FT rig( min( ymax(*(first_pierced_in_BR - 1)), I_R.second));
    
    // first rectangle in s[BT] pierced by bot
    // ( ==> predecessor gives an upper bound for top):
    iterator first_pierced_in_BT( p.first_left_of( RP::BT, bot));
    
    // first rectangle in s[BT] after first_pierced_in_BT
    // not pierced by bot ( ==> gives a lower bound for top):
    iterator first_unpierced_in_BT( p.first_right_of( RP::BT, bot));
    
    // x coordinate of top piercing point:
    FT top;
    
    // first rectangle in s[TL] not pierced by lef
    // ( ==> gives an upper bound for top):
    iterator first_unpierced_in_TL( p.first_above( RP::TL, lef));
    
    // first rectangle in s[TR] not pierced by rig
    // ( ==> gives a lower bound for top):
    iterator first_unpierced_in_TR( p.first_above( RP::TR, rig));
    
    // first rectangle in s[LR] not pierced by lef:
    iterator first_unpierced_by_lef_in_LR( p.first_above( RP::LR, lef));
    
    // corresponding upper bound for top:
    iterator top_upper_bound_by_lef_in_LR(
      p.first_above( RP::TL,
                     min( ymax( *first_unpierced_by_lef_in_LR),
                          d.ymax())));
    
    // first rectangle in s[LR] not pierced by rig:
    iterator first_unpierced_by_rig_in_LR( p.first_above( RP::LR, rig));
    
    // corresponding lower bound for top:
    iterator top_lower_bound_by_rig_in_LR(
      p.first_above( RP::TR,
                     min( ymax( *first_unpierced_by_rig_in_LR),
                          d.ymax())));
    
    // top side of the first rectangle in s[LR]
    // ( > d.ymax, iff s[LR] is empty)
    FT top_side_of_lr( ymax( *(p.begin( RP::LR))));
    CGAL_optimisation_assertion(
      p.is_empty( RP::LR) && top_side_of_lr > d.ymax() ||
        !p.is_empty( RP::LR) && top_side_of_lr <= d.ymax());
    
    // top side of the
    // first rectangle in s[LR] above (not including)
    // the line y = top_side_of_lr
    // ( possibly > d.ymax())
    FT top_side_above_top_side_of_lr(
      ymax( *p.first_above( RP::LR,
                            min( top_side_of_lr, d.ymax()))));
    CGAL_optimisation_assertion( top_side_above_top_side_of_lr >= d.ymin());
    
    // right side of the
    // first rectangle in s[TL] above (not including)
    // the line y = top_side_of_lr
    // ( possibly > d.xmax())
    FT upper_bound_top_from_slr(
      xmax( *p.first_above( RP::TL,
                            min( top_side_above_top_side_of_lr, d.ymax()))));
    CGAL_optimisation_assertion( upper_bound_top_from_slr >= d.xmin());
    
    // left side of the
    // first rectangle in s[TR] above (not including)
    // the line y = top_side_of_lr
    // ( possibly < d.xmin())
    FT lower_bound_top_from_slr(
      xmin( *p.first_above( RP::TR,
                            min( top_side_above_top_side_of_lr, d.ymax()))));
    CGAL_optimisation_assertion( top_side_above_top_side_of_lr >= d.ymin());
    
    
    
    // one-piercing intervall of s[BT]:
    Intervall I_BT;
    
    // piercing intervall of s[T] with respect to lef and rig:
    Intervall I_TLR;
    
    // make sure, lef and rig are in I_L (I_R resp.):
    if ( lef < I_L.first || rig < I_R.first)
      goto next_iteration;
    
    for (;;) {
      #if defined(CGAL_PCENTER_WINDOW_TRACE)
      Wd << CGAL_RED;
      copy( f, l, wout);
      Wd << CGAL_GREEN
         << CGAL_Iso_rectangleC2<FT>( d.min(), d.max())
         << CGAL_BLUE
         << build_point( bot, d.ymin()));
      cerr << "** show bottom point" << endl;
      while ( Wd.read_mouse( dummy_x, dummy_y) != -1) {}
      #endif
      
      // ----------------------------------------------------------
      // compute subintervall I_BT of top side where s[BT] is pierced
      // (this is determined by the choice of bot)
      
      if ( p.begin( RP::BT) == first_pierced_in_BT) {
      
        // bot pierces the first rectangle from s[BT]
        // or no rectangle from s[BT] at all
        CGAL_optimisation_assertion(
          p.is_empty( RP::BT) ||
          !(*(first_pierced_in_BT)).has_on_unbounded_side(
            build_point( bot, d.ymin())) ||
            Xmax()( *first_pierced_in_BT) < bot);
      
        if ( p.end( RP::BT) == first_unpierced_in_BT) {
      
          // bot pierces the last rectangle from s[BT]
          // (if s[BT] is not empty)
          CGAL_optimisation_assertion(
            p.is_empty( RP::BT) ||
            !(*(p.end( RP::BT) - 1)).has_on_unbounded_side(
              build_point( bot, d.ymin())));
      
          // there is nothing left to pierce:
          I_BT = I_T;
      
        } else {
      
          // we have to pierce the first unpierced rectangle
          // as well as the last one:
          I_BT =
            Intervall( max( xmin( *(first_unpierced_in_BT)),
                            xmin( *(p.end( RP::BT) - 1))),
                       min( xmax( *(first_unpierced_in_BT)),
                            xmax( *(p.end( RP::BT) - 1))));
        }
      } else {
      
        // bot does not pierce the first rectangle from s[BT]:
        CGAL_optimisation_assertion(
          ( *(p.begin( RP::BT))).has_on_unbounded_side(
            build_point( bot, d.ymin())));
      
        if ( p.end( RP::BT) != first_unpierced_in_BT) {
          // now way to pierce both with a single point:
          goto next_iteration;
        }
      
        // we have to pierce the first rectangle
        // as well as the last unpierced one:
        /*
        I_BT =
          Intervall( max( xmin( *--(first_pierced_in_BT)),
                          xmin( *(p.begin( RP::BT)))),
                     min( xmax( *--(first_pierced_in_BT)),
                          xmax( *(p.begin( RP::BT)))));
        */
        I_BT =
          Intervall( max( xmin( *(first_pierced_in_BT - 1)),
                          xmin( *(p.begin( RP::BT)))),
                     min( xmax( *(first_pierced_in_BT - 1)),
                          xmax( *(p.begin( RP::BT)))));
      }
      
      // ----------------------------------------------------------
      // test if the intersection of I_T and I_BT is not empty
      // and contains a point that pierces the remaining rectangles
      // from s[TL] and s[TR]:
      
      I_BT.first = max( I_T.first, I_BT.first);
      I_BT.second = min( I_T.second, I_BT.second);
      #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
      cerr << "update bounds for top" << endl;
      #endif
      
      first_unpierced_in_TL =
        p.first_above_downwards( RP::TL, lef, first_unpierced_in_TL);
      
      first_unpierced_in_TR =
        p.first_above_upwards( RP::TR, rig, first_unpierced_in_TR);
      
      I_TLR = Intervall( xmin( *first_unpierced_in_TR),
                         xmax( *first_unpierced_in_TL));
      
      I_TLR.first = max( I_BT.first, I_TLR.first);
      I_TLR.second = min( I_BT.second, I_TLR.second);
      if ( I_TLR.first > I_TLR.second)
        goto next_iteration;
      
      // ----------------------------------------------------------
      // to pierce s[LR] there are at most two possible ways:
      // the lowest rectangle of s[LR] can either be pierced by
      // Point_2( ... lef ...) or Point_2( ... rig ...)
      // so we just try both:
      
      if ( !p.is_empty( RP::LR)) {
        #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
        cerr << "update bounds for s[LR]" << endl;
        #endif
        
        first_unpierced_by_lef_in_LR =
          p.first_above_downwards( RP::LR, lef, first_unpierced_by_lef_in_LR);
        
        top_upper_bound_by_lef_in_LR =
          p.first_above_downwards( RP::TL,
                                   ymax( *first_unpierced_by_lef_in_LR),
                                   top_upper_bound_by_lef_in_LR);
        
        first_unpierced_by_rig_in_LR =
          p.first_above_upwards( RP::LR, rig, first_unpierced_by_rig_in_LR);
        
        top_lower_bound_by_rig_in_LR =
          p.first_above_upwards( RP::TR,
                                 ymax( *first_unpierced_by_rig_in_LR),
                                 top_lower_bound_by_rig_in_LR);
        
        
        // ------------------------------------------------------
        // try to pierce the first rectangle from s[LR] with lef:
        
        if ( lef >= ymin( *p.begin( RP::LR)) &&
             I_L.first <= ymax( *p.begin( RP::LR))) {
        
          FT lef_save( lef);
          FT rig_save( rig);
        
          for (;;) {
            // top side of the first rectangle in s[LR] not pierced by lef:
            FT first_unpierced;
            // its bound for top:
            FT top_bound_first_unpierced;
            
            if ( lef > top_side_of_lr) {
              // lef has to be moved down
              // in order to pierce the first rectangle from s[LR]:
            
              lef = top_side_of_lr;
              I_TLR.second = upper_bound_top_from_slr;
            
              if ( I_TLR.first > I_TLR.second)
                // cannot move down lef without making it impossible
                // to place the top piercing point
                break;
            
              first_unpierced = top_side_above_top_side_of_lr;
              top_bound_first_unpierced = upper_bound_top_from_slr;
            
            } else {
              // no need to move lef down
            
              first_unpierced = ymax( *first_unpierced_by_lef_in_LR);
              top_bound_first_unpierced =
                Xmax()( *top_upper_bound_by_lef_in_LR);
            }
            // try to position rig:
            
            if ( rig > first_unpierced) {
              // rig has to be moved down
              rig = first_unpierced;
              if ( I_R.first > rig ||
                   top_bound_first_unpierced > I_TLR.second)
                // cannot move down rig without either leaving I_R
                // or making it impossible
                // to place the top piercing point
                break;
            } // if ( rig > first_unpierced)
            
            // check, if the last rectangle from s[LR] is pierced
            // hack: sentinel <=> first_unpierced > d.ymax()
            // (in this case we do not have to check anything
            if ( first_unpierced <= d.ymax() &&
                 ymin( *(p.end( RP::LR) - 1)) > rig)
              break;
            top = I_TLR.second;
            // ----------------------------------------------------------
            // We finally found a piercing set!
            
            *o++ = build_point( bot, d.ymin());
            *o++ = build_point( top, d.ymax());
            *o++ = build_point( d.xmin(), lef);
            *o++ = build_point( d.xmax(), rig);
            
            #if defined(CGAL_PCENTER_WINDOW_TRACE)
            Wd << CGAL_RED;
            copy( f, l, wout);
            Wd << CGAL_GREEN
               << CGAL_Iso_rectangle_2< R >( d.min(), d.max())
               << CGAL_BLACK
               << build_point( bot, d.ymin())
               << build_point( top, d.ymax())
               << build_point( d.xmin(), lef)
               << build_point( d.xmax(), rig);
            cerr << "** show piercing set" << endl;
            while ( Wd.read_mouse( dummy_x, dummy_y) != -1) {}
            #endif
            
            ok = true;
            return o;
          } // for (;;)
        
          lef = lef_save;
          rig = rig_save;
        
        } // if ( lef >= ymin( *p.begin( RP::LR)) && ... )
        
        // ------------------------------------------------------
        // try to pierce the first rectangle from s[LR] with rig:
        
        if ( rig >= ymin( *p.begin( RP::LR)) &&
             I_R.first <= ymax( *p.begin( RP::LR))) {
        
          // top side of the first rectangle in s[LR] not pierced by rig:
          FT first_unpierced;
          // its bound for top:
          FT top_bound_first_unpierced;
          
          if ( rig > top_side_of_lr) {
            // rig has to be moved down
            // in order to pierce the first rectangle from s[LR]:
          
            rig = top_side_of_lr;
            I_TLR.first = lower_bound_top_from_slr;
          
            if ( I_TLR.first > I_TLR.second)
              // cannot move down rig without making it impossible
              // to place the top piercing point
              goto next_iteration;
          
            first_unpierced = top_side_above_top_side_of_lr;
            top_bound_first_unpierced = lower_bound_top_from_slr;
          
          } else {
            // no need to move rig down
          
            first_unpierced = ymax( *first_unpierced_by_rig_in_LR);
            top_bound_first_unpierced =
              Xmin()( *top_lower_bound_by_rig_in_LR);
          }
          // try to position lef:
          
          if ( lef > first_unpierced) {
            // lef has to be moved down
            lef = first_unpierced;
            if ( I_L.first > lef ||
                 top_bound_first_unpierced < I_TLR.first)
              // cannot move down lef without either leaving I_L
              // or making it impossible
              // to place the top piercing point
              goto next_iteration;
          } // if ( lef > first_unpierced)
          
          // check, if the last rectangle from s[LR] is pierced
          // hack: sentinel <=> first_unpierced > d.ymax()
          // (in this case we do not have to check anything
          if ( first_unpierced <= d.ymax() &&
               ymin( *(p.end( RP::LR) - 1)) > lef)
            goto next_iteration;
          top = I_TLR.first;
        
        } // if ( rig >= ymin( *p.begin( RP::LR)) && ... )
        
        
      } // if ( !p.is_empty( RP::LR))
      else
        top = I_TLR.first;
      
      // ----------------------------------------------------------
      // We finally found a piercing set!
      
      *o++ = build_point( bot, d.ymin());
      *o++ = build_point( top, d.ymax());
      *o++ = build_point( d.xmin(), lef);
      *o++ = build_point( d.xmax(), rig);
      
      #if defined(CGAL_PCENTER_WINDOW_TRACE)
      Wd << CGAL_RED;
      copy( f, l, wout);
      Wd << CGAL_GREEN
         << CGAL_Iso_rectangle_2< R >( d.min(), d.max())
         << CGAL_BLACK
         << build_point( bot, d.ymin())
         << build_point( top, d.ymax())
         << build_point( d.xmin(), lef)
         << build_point( d.xmax(), rig);
      cerr << "** show piercing set" << endl;
      while ( Wd.read_mouse( dummy_x, dummy_y) != -1) {}
      #endif
      
      ok = true;
      return o;
    
      // ----------------------------------------------------------
      // test the next intervall:
    
      next_iteration:
    #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
      cerr << "!!! next iteration" << endl;
    #endif
    
      if ( xmin( *(first_pierced_in_BR - 1)) <
           xmin( *(first_pierced_in_BT - 1)))
        bot = xmin( *(--first_pierced_in_BR));
      else
        bot = xmin( *(--first_pierced_in_BT));
    
      if ( bot > I_B.second)
        break;
    
      CGAL_optimisation_assertion(
        first_pierced_in_BR >= p.begin( RP::BR) &&
        first_pierced_in_BT >= p.begin( RP::BT));
    
      #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
      cerr << "compute rig" << endl;
      #endif
      
      CGAL_optimisation_assertion( first_pierced_in_BR >= p.begin( RP::BR));
      
      rig = min( ymax( *(first_pierced_in_BR - 1)), I_R.second);
      
      // if that does not fit together with I_R --> continue
      // !!! in case of lef we can abort here
      if ( rig < I_R.first)
        goto next_iteration;
      
      // update first_pierced_in_BL:
      first_pierced_in_BL =
        p.first_right_of( RP::BL, bot, first_pierced_in_BL);
      
      #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
      cerr << "compute lef" << endl;
      #endif
      
      CGAL_optimisation_assertion( first_pierced_in_BL >= p.begin( RP::BL));
      
      lef = min( ymax( *(first_pierced_in_BL - 1)), I_L.second);
      
      // if that does not fit together with I_L --> continue
      // !!! in case of lef we can abort here
      if ( lef < I_L.first)
        goto next_iteration;
      
      #if defined(CGAL_PCENTER_WINDOW_TRACE)
      Wd << CGAL_BLUE
         << build_point( d.xmin(), lef)
         << build_point( d.xmax(), rig);
      cerr << "** show lef and rig" << endl;
      while ( Wd.read_mouse( dummy_x, dummy_y) != -1) {}
      #endif
      #if defined(CGAL_PCENTER_TRACE) && CGAL_PCENTER_TRACE <= 4
      cerr << "update bounds for s[BT]" << endl;
      #endif
      
      first_unpierced_in_BT =
        p.first_right_of( RP::BT, bot, first_unpierced_in_BT);
    
    } // for (;;)
  
  } // if ( p.is_empty( RP::NO))

  ok = false;
  return o;

} // CGAL_four_pierce_rectangles( f, l, d, p, o, ok)

template < class _Traits >
class CGAL_Two_piercing_algorithm {
public:
  // don't touch these typedefs ;)
  // Traits is not enough for sunpro ...
  typedef _Traits                            Traits;
  typedef typename _Traits::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename _Traits::Point_2          Point_2;
  typedef typename _Traits::FT               FT;
  typedef typename _Traits::X                X;
  typedef typename _Traits::Y                Y;
  typedef typename _Traits::Xmin             Xmin;
  typedef typename _Traits::Xmax             Xmax;
  typedef typename _Traits::Ymin             Ymin;
  typedef typename _Traits::Ymax             Ymax;
  typedef typename _Traits::Build_point      Build_point;
  typedef typename _Traits::Build_rectangle  Build_rectangle;

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#else
  typedef vector< Iso_rectangle_2 >::iterator
    RandomAccessIC;
  typedef back_insert_iterator< vector< Point_2 > >
    OutputIterator;

  CGAL_Wastebasket< Point_2 >
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    CGAL_Wastebasket< Point_2 > o,
    bool& ok)
  
  {
    return CGAL_two_pierce_rectangles( f, l, o, ok, Traits());
  }

#endif // CGAL_CFG_NO_MEMBER_TEMPLATES

  OutputIterator
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    OutputIterator o,
    bool& ok)
  
  {
    return CGAL_two_pierce_rectangles( f, l, o, ok, Traits());
  }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES
  OutputIterator
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    const CGAL__Loc_domain< Traits, RandomAccessIC >& d,
    OutputIterator o,
    bool& ok)
  { return CGAL_two_pierce_rectangles( f, l, d, o, ok); }

}; // class CGAL_Two_piercing_algorithm
template < class _Traits >
class CGAL_Three_piercing_algorithm {
public:
  // don't touch these typedefs ;)
  // Traits is not enough for sunpro ...
  typedef _Traits                            Traits;
  typedef typename _Traits::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename _Traits::Point_2          Point_2;
  typedef typename _Traits::FT               FT;
  typedef typename _Traits::X                X;
  typedef typename _Traits::Y                Y;
  typedef typename _Traits::Xmin             Xmin;
  typedef typename _Traits::Xmax             Xmax;
  typedef typename _Traits::Ymin             Ymin;
  typedef typename _Traits::Ymax             Ymax;
  typedef typename _Traits::Build_point      Build_point;
  typedef typename _Traits::Build_rectangle  Build_rectangle;

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#else
  typedef vector< Iso_rectangle_2 >::iterator
    RandomAccessIC;
  typedef back_insert_iterator< vector< Point_2 > >
    OutputIterator;

  CGAL_Wastebasket< Point_2 >
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    CGAL_Wastebasket< Point_2 > o,
    bool& ok)
  
  {
    return CGAL_three_pierce_rectangles( f, l, o, ok, Traits());
  }

#endif // CGAL_CFG_NO_MEMBER_TEMPLATES

  OutputIterator
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    OutputIterator o,
    bool& ok)
  
  {
    return CGAL_three_pierce_rectangles( f, l, o, ok, Traits());
  }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES
  OutputIterator
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    const CGAL__Loc_domain< Traits, RandomAccessIC >& d,
    OutputIterator o,
    bool& ok)
  { return CGAL_three_pierce_rectangles( f, l, d, o, ok); }

}; // class CGAL_Three_piercing_algorithm
template < class _Traits >
class CGAL_Four_piercing_algorithm {
public:
  // don't touch these typedefs ;)
  // Traits is not enough for sunpro ...
  typedef _Traits                            Traits;
  typedef typename _Traits::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename _Traits::Point_2          Point_2;
  typedef typename _Traits::FT               FT;
  typedef typename _Traits::X                X;
  typedef typename _Traits::Y                Y;
  typedef typename _Traits::Xmin             Xmin;
  typedef typename _Traits::Xmax             Xmax;
  typedef typename _Traits::Ymin             Ymin;
  typedef typename _Traits::Ymax             Ymax;
  typedef typename _Traits::Build_point      Build_point;
  typedef typename _Traits::Build_rectangle  Build_rectangle;

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#else
  typedef vector< Iso_rectangle_2 >::iterator
    RandomAccessIC;
  typedef back_insert_iterator< vector< Point_2 > >
    OutputIterator;

  CGAL_Wastebasket< Point_2 >
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    CGAL_Wastebasket< Point_2 > o,
    bool& ok)
  
  {
    return CGAL_four_pierce_rectangles( f, l, o, ok, Traits());
  }

#endif // CGAL_CFG_NO_MEMBER_TEMPLATES

  OutputIterator
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    OutputIterator o,
    bool& ok)
  
  {
    return CGAL_four_pierce_rectangles( f, l, o, ok, Traits());
  }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES
  OutputIterator
  operator()(
    RandomAccessIC f,
    RandomAccessIC l,
    const CGAL__Loc_domain< Traits, RandomAccessIC >& d,
    OutputIterator o,
    bool& ok)
  { return CGAL_four_pierce_rectangles( f, l, d, o, ok); }

}; // class CGAL_Four_piercing_algorithm

#endif // ! (CGAL_PIERCE_RECTANGLES_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

