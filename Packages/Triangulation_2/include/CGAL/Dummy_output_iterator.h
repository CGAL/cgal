// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Dummy_output_iterator.h
// package       : Triangulation 
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_DUMMY_OUTPUT_ITERATOR_H
#define CGAL_DUMMY_OUTPUT_ITERATOR_H

CGAL_BEGIN_NAMESPACE

class Dummy_output_iterator
{
public:
  Dummy_output_iterator() {}
  Dummy_output_iterator( const Dummy_output_iterator&) {}

  template<class T>
  Dummy_output_iterator& operator=(const T&)    {return  *this; }

  Dummy_output_iterator&  operator++() {return  *this;}
  Dummy_output_iterator&  operator++(int) {return  *this;}
  
 
   Dummy_output_iterator&  operator*() { return *this ;} 
};


CGAL_END_NAMESPACE
#endif //CGAL_DUMMY_OUTPUT_ITERATOR_H

