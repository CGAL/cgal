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

#include <CGAL/Optimisation/assertions.h>
#include <CGAL/circulator_bases.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class OutputIterator, class Operation >
class Transform_iterator : public CGAL_STD::output_iterator {
public:
  typedef Transform_iterator< OutputIterator, Operation >
    self;
  typedef typename Operation::argument_type
    argument_type;

  Transform_iterator( const OutputIterator& o,
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
Transform_iterator< OutputIterator, Operation >
transform_iterator( const OutputIterator& o,
                         const Operation& op)
{ return Transform_iterator< OutputIterator, Operation >( o, op); }

template < class OutputIterator, class Operation > inline
std::output_iterator_tag
iterator_category(
  const Transform_iterator< OutputIterator, Operation >&)
{ return output_iterator_tag(); }

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

