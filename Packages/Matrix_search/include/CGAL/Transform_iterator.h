// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
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
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// An OutputIterator Adaptor applying an unary function
// ============================================================================

#if ! (CGAL_TRANSFORM_ITERATOR_H)
#define CGAL_TRANSFORM_ITERATOR_H 1

#include <CGAL/Optimisation/assertions.h>
#include <CGAL/circulator_bases.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class OutputIterator, class Operation >
struct Transform_iterator {
  typedef std::output_iterator_tag             iterator_category;
  typedef Transform_iterator< OutputIterator, Operation >   self;
  typedef typename Operation::argument_type        argument_type;

  Transform_iterator( const OutputIterator& o,
                      const Operation& op)
    : o_( o), op_( op)
  {}

  operator OutputIterator() { return o_; }

  self& operator*() { return *this; }

  self& operator++() { return *this; }

  self& operator++( int) { return *this; }

  self& operator=( const argument_type& a) {
    *(o_++) = op_( a);
    return *this;
  }

private:
  OutputIterator o_;
  Operation      op_;
};

template < class OutputIterator, class Operation > inline
Transform_iterator< OutputIterator, Operation >
transform_iterator( const OutputIterator& o,
                         const Operation& op)
{ return Transform_iterator< OutputIterator, Operation >( o, op); }

template < class OutputIterator, class Operation > inline
Iterator_tag
query_circulator_or_iterator(
  const Transform_iterator< OutputIterator, Operation >&)
{ return Iterator_tag(); }

CGAL_END_NAMESPACE

#endif // ! (CGAL_TRANSFORM_ITERATOR_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

