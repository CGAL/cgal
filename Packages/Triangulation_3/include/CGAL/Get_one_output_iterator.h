// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Get_one_output_iterator.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : INRIA Sophia Antipolis <Mariette.Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_GET_ONE_OUTPUT_ITERATOR_H
#define CGAL_GET_ONE_OUTPUT_ITERATOR_H

CGAL_BEGIN_NAMESPACE

// This output iterator stores a reference to an object T, which will be
// affected by operator*().  Its operator++() does nothing.

template <class T>
class Get_one_output_iterator
{
  T & _t;
public:
  Get_one_output_iterator(T & tt) : _t(tt) {}

  Get_one_output_iterator& operator++() {return  *this;}
  Get_one_output_iterator& operator++(int) {return  *this;}
  
  T & operator*() { return _t;} 
};

CGAL_END_NAMESPACE

#endif // CGAL_GET_ONE_OUTPUT_ITERATOR_H
