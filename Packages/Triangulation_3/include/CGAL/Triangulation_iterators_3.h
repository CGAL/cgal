// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_iterators_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_ITERATORS_3_H
#define CGAL_TRIANGULATION_ITERATORS_3_H

#include <CGAL/basic.h>
#include <utility>
#include <iterator>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/triangulation_assertions.h>

// This file contains the Finite_*_iterators.

CGAL_BEGIN_NAMESPACE

// Wraps an iterator, but skips infinite elements.
// This one works for Vertex and Cell.
template <class Tr, class Iterator_base>
class Triangulation_finite_iterator_3
{
  typedef Triangulation_finite_iterator_3<Tr, Iterator_base> Iterator;
public:
  typedef typename Iterator_base::value_type  value_type;
  typedef typename Iterator_base::pointer     pointer;
  typedef typename Iterator_base::reference   reference;
  typedef std::size_t                         size_type;
  typedef std::ptrdiff_t                      difference_type;
  typedef std::bidirectional_iterator_tag     iterator_category;

  Triangulation_finite_iterator_3()
    : _ib(), _end(), _tr(NULL)
  {}

  Triangulation_finite_iterator_3(const Tr *tr,
	                          Iterator_base it,
	                          Iterator_base end)
    : _ib(it), _end(end), _tr(tr)
  { 
    while ( _tr->is_infinite(&*_ib) )
      ++_ib;
  }

  // Past the end iterator.
  Triangulation_finite_iterator_3(const Tr *tr, Iterator_base end)
    : _ib(end), _end(end), _tr(tr)
  {}

  bool
  operator==(const Iterator & it) const
  {
      return _ib == it._ib && _end == it._end && _tr == it._tr;
  }

  bool
  operator!=(const Iterator & it) const
  {
      return !(*this == it);
  }

  Iterator &
  operator++()
  {
      do {
	++_ib; 
      } while ( _ib != _end && _tr->is_infinite(&*_ib) );
      return *this;
  }

  Iterator &
  operator--()
  {
      do {
	--_ib;
      } while ( _ib != _end && _tr->is_infinite(&*_ib) );
      return *this;
  }

  Iterator
  operator++(int)
  {
    Iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Iterator
  operator--(int)
  {
    Iterator tmp(*this);
    --(*this);
    return tmp;
  }

  reference operator*() const
  {
    return *_ib;
  }

  pointer operator->() const
  {
    return &*_ib;
  }
     
private:
  Iterator_base _ib;
  Iterator_base _end;
  const Tr * _tr;
};

// This second version is for Facet/Edge : doesn't support operator->().
template <class Tr, class Iterator_base>
class Triangulation_finite_iterator2_3
{
  typedef Triangulation_finite_iterator2_3<Tr, Iterator_base> Iterator;
public:
  typedef typename Iterator_base::value_type  value_type;
  typedef typename Iterator_base::pointer     pointer;
  typedef typename Iterator_base::reference   reference;
  typedef std::size_t                         size_type;
  typedef std::ptrdiff_t                      difference_type;
  typedef std::bidirectional_iterator_tag     iterator_category;

  Triangulation_finite_iterator2_3()
    : _ib(), _end(), _tr(NULL)
  {}

  Triangulation_finite_iterator2_3(const Tr *tr,
	                          Iterator_base it,
	                          Iterator_base end)
    : _ib(it), _end(end), _tr(tr)
  { 
    while ( _tr->is_infinite(*_ib) )
      ++_ib;
  }

  // Past the end iterator.
  Triangulation_finite_iterator2_3(const Tr *tr, Iterator_base end)
    : _ib(end), _end(end), _tr(tr)
  {}

  bool
  operator==(const Iterator & it) const
  {
      return _ib == it._ib && _end == it._end && _tr == it._tr;
  }

  bool
  operator!=(const Iterator & it) const
  {
      return !(*this == it);
  }

  Iterator &
  operator++()
  {
      do {
	++_ib; 
      } while ( _ib != _end && _tr->is_infinite(*_ib) );
      return *this;
  }

  Iterator &
  operator--()
  {
      do {
	--_ib;
      } while ( _ib != _end && _tr->is_infinite(*_ib) );
      return *this;
  }

  Iterator
  operator++(int)
  {
    Iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Iterator
  operator--(int)
  {
    Iterator tmp(*this);
    --(*this);
    return tmp;
  }

  value_type operator*() const
  {
    return *_ib;
  }

private:
  Iterator_base _ib;
  Iterator_base _end;
  const Tr * _tr;
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_ITERATORS_3_H
