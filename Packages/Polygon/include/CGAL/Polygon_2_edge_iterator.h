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
// author(s)     : Wieger Wesselink, Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_2_EDGE_ITERATOR_H
#define CGAL_POLYGON_2_EDGE_ITERATOR_H

#include <iterator>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Polygon_2_edge_iterator
//-----------------------------------------------------------------------//
// Ideally the class Polygon_2_edge_iterator would be a nested class of
// Polygon_2, but this leads to compiler problems with SGI C++ 4.0
// with the iterator_category() function

template <class Traits_, class Container_> class Polygon_2;

template <class Segment_>
class Polygon_2__Segment_ptr
{
public:
    typedef Segment_ Segment;
    Polygon_2__Segment_ptr(Segment const &seg) :m_seg(seg){}
    Segment* operator->() {return &m_seg;}
private:
    Segment m_seg;
};

template <class Traits_, class Container_>
class Polygon_2_edge_iterator {
  public:
#ifdef CGAL_CFG_NO_LAZY_INSTANTIATION
    typedef std::forward_iterator_tag iterator_category;
#else
    typedef std::random_access_iterator_tag iterator_category;
#endif
    typedef typename Traits_::Segment_2 Segment_2;
    typedef typename Traits_::Segment_2 value_type;
    typedef Container_ Container;
    typedef typename Container_::const_iterator const_iterator;
    typedef typename Container_::difference_type difference_type;
    typedef Segment_2*           pointer;
    typedef Segment_2&           reference;
  private:
    const Container_* container;   // needed for dereferencing the last edge
    const_iterator first_vertex;   // points to the first vertex of the edge
  public:
    Polygon_2_edge_iterator() {}
    Polygon_2_edge_iterator(const Container_* c, const_iterator f)
      : container(c), first_vertex(f) {}

    bool operator==(
      const Polygon_2_edge_iterator<Traits_, Container_>& x) const
    {
      return first_vertex == x.first_vertex;
    }
    
    bool operator!=(
      const Polygon_2_edge_iterator<Traits_, Container_>& x) const
    {
      return !(first_vertex == x.first_vertex);
    }

    Segment_2 operator*() {
      const_iterator second_vertex = first_vertex;
      ++second_vertex;
      if (second_vertex == container->end())
        second_vertex = container->begin();
      typename Traits_::Construct_segment_2 construct_segment_2 = 
            Traits_().construct_segment_2_object();
      return construct_segment_2(*first_vertex, *second_vertex);
    }
    
    Polygon_2__Segment_ptr<Segment_2> operator->()
        {return Polygon_2__Segment_ptr<Segment_2>(operator*());}

    Polygon_2_edge_iterator<Traits_, Container_>& operator++() {
      ++first_vertex;
      return *this;
    }

    Polygon_2_edge_iterator<Traits_, Container_> operator++(int) {
      Polygon_2_edge_iterator<Traits_, Container_> tmp = *this;
      ++*this;
      return tmp;
    }

    Polygon_2_edge_iterator<Traits_, Container_>& operator--() {
      --first_vertex;
      return *this;
    }

    Polygon_2_edge_iterator<Traits_, Container_> operator--(int) {
      Polygon_2_edge_iterator<Traits_, Container_> tmp = *this;
      --*this;
      return tmp;
    }

#ifndef CGAL_CFG_NO_LAZY_INSTANTIATION
// random access iterator requirements
    Polygon_2_edge_iterator<Traits_, Container_>&
    operator+=(difference_type n) {
      first_vertex += n;
      return *this;
    }

    Polygon_2_edge_iterator<Traits_, Container_>
    operator+(difference_type n) {
      return Polygon_2_edge_iterator<Traits_, Container_>(
        container, first_vertex + n);
    }

    Polygon_2_edge_iterator<Traits_, Container_>&
    operator-=(difference_type n) {
      return (*this) -= n;
    }

    Polygon_2_edge_iterator<Traits_, Container_>
    operator-(difference_type n) {
      return Polygon_2_edge_iterator<Traits_, Container_>(
        container, first_vertex - n);
    }

    difference_type
    operator-(const Polygon_2_edge_iterator<Traits_, Container_>& a) {
      return first_vertex - a.first_vertex;
    }

    Segment_2 operator[](int n) {
      return *Polygon_2_edge_iterator<Traits_, Container_>(
        container, first_vertex+n);
    }

    bool operator<(const Polygon_2_edge_iterator<Traits_, Container_>& a)
    {
      return first_vertex < a.first_vertex;
    }

    bool operator>(const Polygon_2_edge_iterator<Traits_, Container_>& a)
    {
      return first_vertex > a.first_vertex;
    }

    bool operator<=(const Polygon_2_edge_iterator<Traits_, Container_>& a)
    {
      return first_vertex <= a.first_vertex;
    }

    bool operator>=(const Polygon_2_edge_iterator<Traits_, Container_>& a)
    {
      return first_vertex >= a.first_vertex;
    }
#endif

};


template <class Traits_,  class Container_>
typename Container_::difference_type
distance_type(const Polygon_2_edge_iterator<Traits_,Container_>&)
{ return Container_::difference_type(); }

template <class Traits_,  class Container_>
typename Traits_::Segment_2*
value_type(const Polygon_2_edge_iterator<Traits_,Container_>&)
{ return (typename Traits_::Segment_2 *)(0); }


//-----------------------------------------------------------------------//
//                          implementation
//-----------------------------------------------------------------------//

//--------------------------------------------------------------------//
// Comment: the iterator category of a Polygon_2_edge_iterator should
// be equal to the iterator category of the corresponding container, but this
// cannot be implemented (???).
//--------------------------------------------------------------------//

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS 
template <class Traits_, class Container_>
inline
bidirectional_iterator_tag
iterator_category(const Polygon_2_edge_iterator<Traits_, Container_>&)
{
  return bidirectional_iterator_tag();
}

#endif // CGAL_CFG_NO_ITERATOR_TRAITS

//--------------------------------------------------------------------//
// I don't know how to implement the following function:
//
// template <class Traits_, class Container_>
// inline
// Polygon_2_edge_iterator<Traits_, Container_>
// operator+(Container_::difference_type n,
//           Polygon_2_edge_iterator<Traits_, Container_>& a)
// { return a+n; }
//--------------------------------------------------------------------//

CGAL_END_NAMESPACE

#endif

