// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Filter_iterator_with_handle.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_FILTER_ITERATOR_WITH_HANDLE_H
#define CGAL_FILTER_ITERATOR_WITH_HANDLE_H

#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE 


template <class It, class Filter, class Handle>
class Filter_iterator_with_handle
  : public Filter_iterator<It, Filter>
{
public:
  typedef Filter_iterator<It, Filter> Base;

  Filter_iterator_with_handle()
    : Base()
  {}

  Filter_iterator_with_handle(It b, It e, Filter f)
    : Base(b, e, f)
  {}
  Filter_iterator_with_handle(It b, It e, Filter f, It c)
    : Base(b, e, f, c)
  {}


  operator Handle() const {return (*this)->handle();}  
  Self&  operator++() { Base::operator++(); return *this;}
  Self&  operator--() { Base::operator--(); return *this;}
  Self   operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self   operator--(int) { Self tmp(*this); --(*this); return tmp; }
};

template < class I, class P, class H >
inline Filter_iterator_with_handle< I, P, H >
filter_iterator_with_handle(I b, I e, const P& p, H h)
{ return Filter_iterator_with_handle< I, P, H >(b, e, p); }

template < class I, class P, class H >
inline Filter_iterator_with_handle< I, P, H >
filter_iterator_with_handle(I b, I e, const P& p, H h, I c)
{ return Filter_iterator_with_handle< I, P, H >(b, e, p, c); }




CGAL_END_NAMESPACE 

#endif // CGAL_FILTER_ITERATOR_WITH_HANDLE_H
