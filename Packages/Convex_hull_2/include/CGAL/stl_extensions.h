// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 03
//
// file          : include/CGAL/stl_extensions.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra 
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_STL_EXTENSIONS_H
#define CGAL_STL_EXTENSIONS_H

#include <iterator>
#include <vector>
#include <utility>

CGAL_BEGIN_NAMESPACE
template <class ForwardIterator>
inline
ForwardIterator
successor( ForwardIterator it )
{
 return ++it;
}

template <class BidirectionalIterator>
inline
BidirectionalIterator
predecessor( BidirectionalIterator it )
{
 return --it;
}

template <class InputIterator, class OutputIterator, class UnaryPredicate>
OutputIterator
copy_if( InputIterator first, InputIterator last,
              OutputIterator  result,
              UnaryPredicate  pred );

template <class InputIterator, class OutputIterator1, 
          class OutputIterator2, class UnaryPredicate>
std::pair<OutputIterator1,OutputIterator2>
copy_if_else( InputIterator first, InputIterator last,
                   OutputIterator1 result1,
                   OutputIterator2 result2,
                   UnaryPredicate  pred );

CGAL_END_NAMESPACE
#include <CGAL/IO/Tee_for_output_iterator.h>

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/stl_extensions.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION


#endif // CGAL_STL_EXTENSIONS_H

