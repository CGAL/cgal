// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 03
//
// file          : include/CGAL/IO/Tee_for_output_iterator.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_IO_TEE_FOR_OUTPUT_ITERATOR_H
#define CGAL_IO_TEE_FOR_OUTPUT_ITERATOR_H

#include <iterator>
#include <vector>
#include <CGAL/Handle.h>

CGAL_BEGIN_NAMESPACE
template <class T> class _Tee_for_output_iterator_rep;

template <class OutputIterator, class T>
class Tee_for_output_iterator 
  : public Handle
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  , public std::output_iterator
#endif // CGAL_CFG_NO_ITERATOR_TRAITS
{
  typedef std::vector<T>                             container;
  typedef typename container::iterator               iterator;
  typedef T                                          value_type;
  typedef std::output_iterator_tag                   iterator_category;
#ifndef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef std::iterator_traits<iterator>             iter_traits;
  typedef typename iter_traits::pointer              pointer;
  typedef typename iter_traits::reference            reference;
#endif // CGAL_CFG_NO_ITERATOR_TRAITS

public:
  Tee_for_output_iterator(const OutputIterator& o) : o_it(o) 
  {  PTR = (Rep*) new _Tee_for_output_iterator_rep<T>(); }

  Tee_for_output_iterator<OutputIterator,T>& 
  operator=(const T& value) 
  { 
    ptr()->output_so_far.push_back(value);
    *o_it = value;
    return *this;
  }

  Tee_for_output_iterator<OutputIterator,T>& 
  operator*() 
  { return *this; }

  Tee_for_output_iterator<OutputIterator,T>& 
  operator++() 
  { 
    ++o_it; 
    return *this; 
  } 

  Tee_for_output_iterator<OutputIterator,T> 
  operator++(int) 
  { 
    Tee_for_output_iterator<OutputIterator,T> tmp = *this;
    o_it++; 
    return tmp; 
  } 

  iterator
  output_so_far_begin()
  { return ptr()->output_so_far.begin(); }

  iterator
  output_so_far_end()
  { return ptr()->output_so_far.end(); }

  OutputIterator&
  to_output_iterator()
  { return o_it; }

  _Tee_for_output_iterator_rep<T>*
  ptr()
  { return (_Tee_for_output_iterator_rep<T>*)PTR; }

protected:
  OutputIterator o_it;
};

template <class T>
class _Tee_for_output_iterator_rep : public Rep
{
public:
  std::vector<T> output_so_far;
};

template <class OutputIterator, class T>
inline 
T*
value_type(const Tee_for_output_iterator<OutputIterator,T>&)
{ return (T*)0; }
CGAL_END_NAMESPACE

template <class OutputIterator, class T>
inline 
std::output_iterator_tag
iterator_category(const CGAL::Tee_for_output_iterator<OutputIterator,T>&)
{ return std::output_iterator_tag(); }


#endif // CGAL_IO_TEE_FOR_OUTPUT_ITERATOR_H
