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
// file          : include/CGAL/Polygon_2_edge_iterator.h
// source        : 
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_2_EDGE_ITERATOR_H
#define CGAL_POLYGON_2_EDGE_ITERATOR_H

#include <iterator.h>
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif // CGAL_CIRCULATOR_H

//-----------------------------------------------------------------------//
//                          CGAL_Polygon_2_edge_iterator
//-----------------------------------------------------------------------//
// Ideally the class CGAL_Polygon_2_edge_iterator would be a nested class of
// CGAL_Polygon_2, but this leads to compiler problems with SGI C++ 4.0
// with the iterator_category() function

template <class _Traits, class _Container> class CGAL_Polygon_2;

template <class _Traits, class _Container>
class CGAL_Polygon_2_edge_iterator {
  public:
    typedef _Traits Traits;
    typedef typename _Traits::Segment_2 Segment_2;
    typedef _Container Container;
    typedef typename _Container::const_iterator const_iterator;
    typedef typename _Container::difference_type difference_type;
  private:
    const _Container* container;   // needed for dereferencing the last edge
    const_iterator first_vertex;   // points to the first vertex of the edge
  public:
    CGAL_Polygon_2_edge_iterator() {}
    CGAL_Polygon_2_edge_iterator(const _Container* c, const_iterator f)
      : container(c), first_vertex(f) {}

    bool operator==(
      const CGAL_Polygon_2_edge_iterator<_Traits, _Container>& x) const
    {
      return first_vertex == x.first_vertex;
    }

    Segment_2 operator*() {
      const_iterator second_vertex = first_vertex;
      ++second_vertex;
      if (second_vertex == container->end())
        second_vertex = container->begin();
      return Segment_2(*first_vertex, *second_vertex);
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container>& operator++() {
      ++first_vertex;
      return *this;
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container> operator++(int) {
      CGAL_Polygon_2_edge_iterator<_Traits, _Container> tmp = *this;
      ++*this;
      return tmp;
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container>& operator--() {
      --first_vertex;
      return *this;
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container> operator--(int) {
      CGAL_Polygon_2_edge_iterator<_Traits, _Container> tmp = *this;
      --*this;
      return tmp;
    }

#ifndef CGAL_CFG_NO_LAZY_INSTANTIATION
// random access iterator requirements
    CGAL_Polygon_2_edge_iterator<_Traits, _Container>&
    operator+=(difference_type n) {
      first_vertex += n;
      return *this;
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container>
    operator+(difference_type n) {
      return CGAL_Polygon_2_edge_iterator<_Traits, _Container>(
        container, first_vertex + n);
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container>&
    operator-=(difference_type n) {
      return (*this) -= n;
    }

    CGAL_Polygon_2_edge_iterator<_Traits, _Container>
    operator-(difference_type n) {
      return CGAL_Polygon_2_edge_iterator<_Traits, _Container>(
        container, first_vertex - n);
    }

    difference_type
    operator-(const CGAL_Polygon_2_edge_iterator<_Traits, _Container>& a) {
      return first_vertex - a.first_vertex;
    }

    Segment_2 operator[](int n) {
      return *CGAL_Polygon_2_edge_iterator<_Traits, _Container>(
        container, first_vertex+n);
    }

    bool operator<(const CGAL_Polygon_2_edge_iterator<_Traits, _Container>& a)
    {
      return first_vertex < a.first_vertex;
    }

    bool operator>(const CGAL_Polygon_2_edge_iterator<_Traits, _Container>& a)
    {
      return first_vertex > a.first_vertex;
    }

    bool operator<=(const CGAL_Polygon_2_edge_iterator<_Traits, _Container>& a)
    {
      return first_vertex <= a.first_vertex;
    }

    bool operator>=(const CGAL_Polygon_2_edge_iterator<_Traits, _Container>& a)
    {
      return first_vertex >= a.first_vertex;
    }
#endif

    // the global function distance_type can only be defined as a friend
    // function within the class due to compiler problems with g++
    friend difference_type
    distance_type(const CGAL_Polygon_2_edge_iterator<_Traits,_Container>&)
      { return difference_type(); }

    friend
    Segment_2* value_type(
      const CGAL_Polygon_2_edge_iterator<_Traits,_Container>&)
    { return (Segment_2*)(0); }
};

//-----------------------------------------------------------------------//
//                          implementation
//-----------------------------------------------------------------------//

//--------------------------------------------------------------------//
// Comment: the iterator category of a CGAL_Polygon_2_edge_iterator should
// be equal to the iterator category of the corresponding container, but this
// cannot be implemented (???).
//--------------------------------------------------------------------//

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS 
template <class _Traits, class _Container>
inline
bidirectional_iterator_tag
iterator_category(const CGAL_Polygon_2_edge_iterator<_Traits, _Container>&)
{
  return bidirectional_iterator_tag();
}

#endif // CGAL_CFG_NO_ITERATOR_TRAITS

//--------------------------------------------------------------------//
// I don't know how to implement the following function:
//
// template <class _Traits, class _Container>
// inline
// CGAL_Polygon_2_edge_iterator<_Traits, _Container>
// operator+(_Container::difference_type n,
//           CGAL_Polygon_2_edge_iterator<_Traits, _Container>& a)
// { return a+n; }
//--------------------------------------------------------------------//

#endif

