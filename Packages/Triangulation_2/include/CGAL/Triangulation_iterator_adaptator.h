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
// file          : include/CGAL/Triangulation_iterator_adaptator.h
// package       : Triangulation_2 
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_TRIANGULATION_ITERATOR_ADAPTATOR_H
#define CGAL_TRIANGULATION_ITERATOR_ADAPTATOR_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE

template <class Base, class Handle>
struct Triangulation_iterator_handle_adaptor
  : public Base
{
  typedef Triangulation_iterator_handle_adaptor<Base, Handle>  Self;
  Triangulation_iterator_handle_adaptor() : Base() {}
   
  Triangulation_iterator_handle_adaptor(const Base & b) 
    : Base(b) {}

  operator Handle() const {return (*this)->handle();}
  Self&  operator++() { Base::operator++(); return *this;}
  Self&  operator--() { Base::operator--(); return *this;}
  Self   operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self   operator--(int) { Self tmp(*this); --(*this); return tmp; }
};



CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_ITERATOR_ADAPTATOR_H
