#line 1354 "mon_search.aw"
#line 18 "code_formatting.awi"
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

#line 1358 "mon_search.aw"
#line 54 "code_formatting.awi"
#if ! (TRANSFORM_ITERATOR_H)
#define TRANSFORM_ITERATOR_H 1

#line 269 "mon_search.aw"
#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H
#ifndef CGAL_CIRCULATOR_BASES_H
#include <CGAL/circulator_bases.h>
#endif // CGAL_CIRCULATOR_BASES_H
#ifndef CGAL_PROTECT_ITERATOR
#include <iterator>
#define CGAL_PROTECT_ITERATOR
#endif

#line 46 "code_formatting.awi"
CGAL_BEGIN_NAMESPACE
#line 279 "mon_search.aw"

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

#line 50 "code_formatting.awi"
CGAL_END_NAMESPACE
#line 335 "mon_search.aw"

#endif // ! (TRANSFORM_ITERATOR_H)

#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

