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
// file          : include/CGAL/ch_melkman.h
// package       : Convex_hull_2 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_MELKMAN_H
#define CGAL_CH_MELKMAN_H

#include <CGAL/basic.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_melkman( InputIterator first, InputIterator last,
            OutputIterator result, const Traits& ch_traits);


template <class InputIterator, class OutputIterator>
OutputIterator
ch_melkman( InputIterator first, InputIterator last,  OutputIterator result)
{ 
    typedef std::iterator_traits<InputIterator>   ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_melkman( first, last, result, Kernel()); 
}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_melkman.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_MELKMAN_H

