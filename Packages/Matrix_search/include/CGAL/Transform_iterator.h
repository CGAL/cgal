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
// file          : Transform_iterator.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// An OutputIterator Adaptor applying an unary function
// ============================================================================

#if ! (CGAL_TRANSFORM_ITERATOR_H)
#define CGAL_TRANSFORM_ITERATOR_H 1

#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H
#ifndef CGAL_CIRCULATOR_BASES_H
#include <CGAL/circulator_bases.h>
#endif // CGAL_CIRCULATOR_BASES_H
#ifndef CGAL_PROTECT_ITERATOR_H
#include <iterator.h>
#define CGAL_PROTECT_ITERATOR_H
#endif // CGAL_PROTECT_ITERATOR_H

template < class OutputIterator, class Operation >
class CGAL_Transform_iterator : public output_iterator {
public:
  typedef CGAL_Transform_iterator< OutputIterator, Operation >
    self;
  typedef typename Operation::argument_type
    argument_type;

  CGAL_Transform_iterator( const OutputIterator& o,
                           const Operation& op)
    : _o( o), _op( op)
  {}

  operator OutputIterator()
  { return _o; }

  self& operator*()
  { return *this; }

  self& operator++()
  { return *this; }

  self& operator++( int)
  { return *this; }

  self& operator=( const argument_type& a)
  {
    *(_o++) = _op( a);
    return *this;
  }

private:
  OutputIterator _o;
  Operation      _op;
};

template < class OutputIterator, class Operation > inline
CGAL_Transform_iterator< OutputIterator, Operation >
CGAL_transform_iterator( const OutputIterator& o,
                         const Operation& op)
{ return CGAL_Transform_iterator< OutputIterator, Operation >( o, op); }

template < class OutputIterator, class Operation > inline
output_iterator_tag
iterator_category(
  const CGAL_Transform_iterator< OutputIterator, Operation >&)
{ return output_iterator_tag(); }

template < class OutputIterator, class Operation > inline
CGAL_Iterator_tag
CGAL_query_circulator_or_iterator(
  const CGAL_Transform_iterator< OutputIterator, Operation >&)
{ return CGAL_Iterator_tag(); }    

#endif // ! (CGAL_TRANSFORM_ITERATOR_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

