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
// file          : include/CGAL/ch_akl_toussaint.h
// package       : Convex_hull_2 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_AKL_TOUSSAINT_H
#define CGAL_CH_AKL_TOUSSAINT_H

#include <CGAL/basic.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

// same as |convex_hull_2(first,last,result)|.
// {\sc traits}: operates on |Traits::Point_2| using |Traits::Less_xy_2|, 
// |Traits::Less_yx_2|, |Traits::Equal_2| and |Traits::Left_turn_2|.
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator last, 
                 OutputIterator  result,
                 const Traits&   ch_traits);

template <class ForwardIterator, class OutputIterator>
inline
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator last, 
                 OutputIterator  result)
{
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_akl_toussaint( first, last, result, Kernel()); 
}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_akl_toussaint.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_AKL_TOUSSAINT_H

