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
// file          : stl_extensions.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : stl_extensions.lw
// revision      : 3.2.1
// revision_date : 19 Apr 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_STL_EXTENSIONS_C
#define CGAL_STL_EXTENSIONS_C

#include <CGAL/stl_extensions.h>

CGAL_BEGIN_NAMESPACE
template <class InputIterator, class OutputIterator, class UnaryPredicate>
OutputIterator
copy_if( InputIterator first, InputIterator last,
              OutputIterator  result,
              UnaryPredicate  pred )
{
 for ( ; first != last; first++ )
 {
    if ( pred(*first) )  *result++ = *first;
 }
 return result;
}


template <class InputIterator, class OutputIterator1, 
          class OutputIterator2, class UnaryPredicate>
std::pair<OutputIterator1,OutputIterator2>
copy_if_else( InputIterator first, InputIterator last,
                   OutputIterator1 result1,
                   OutputIterator2 result2,
                   UnaryPredicate  pred )
{
 for ( ; first != last; first++ )
 {
    if ( pred(*first) )
    {
        *result1++ = *first;
    }
    else
    {
     
        *result2++ = *first;
    }
 }
 return std::pair<OutputIterator1,OutputIterator2>(result1,result2);
}

CGAL_END_NAMESPACE

#endif // CGAL_STL_EXTENSIONS_C

