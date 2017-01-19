// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_TRANSFORM_ITERATOR_H
#define CGAL_TRANSFORM_ITERATOR_H

#include <CGAL/license/Matrix_search.h>


#include <CGAL/Optimisation/assertions.h>
#include <CGAL/circulator_bases.h>
#include <iterator>

namespace std {
  struct _Unchecked_iterator_tag;
}


namespace CGAL {

template < class OutputIterator, class Operation >
struct Transform_iterator {
  // Workaround. Added this non standard iterator category for VC8.
  // Strange that no other iterator complains about this "feature" missing  
  typedef std::_Unchecked_iterator_tag _Checked_iterator_category;
  typedef std::output_iterator_tag             iterator_category;
  typedef Transform_iterator< OutputIterator, Operation >   self;
  typedef typename Operation::argument_type        argument_type;

  typedef typename std::iterator_traits<OutputIterator>::difference_type difference_type;
  typedef typename std::iterator_traits<OutputIterator>::value_type      value_type;
  typedef typename std::iterator_traits<OutputIterator>::pointer         pointer;
  typedef typename std::iterator_traits<OutputIterator>::reference       reference;

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

} //namespace CGAL

#endif // CGAL_TRANSFORM_ITERATOR_H
