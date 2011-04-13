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
// release_date  : 
//
// file          : include/CGAL/ch_bykat.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_BYKAT_H
#define CGAL_CH_BYKAT_H

#include <CGAL/ch_utils.h>
#include <CGAL/ch_value_type.h>

#ifdef CGAL_REP_CLASS_DEFINED
#ifdef STL_GCC
#ifndef GNU_ISTREAM_ITERATOR_VALUE_TYPE_FIX_H
#include <CGAL/gnu_istream_iterator_value_type_fix.h>
#endif // GNU_ISTREAM_ITERATOR_VALUE_TYPE_FIX_H
#endif // STL_GCC
#endif // CGAL_REP_CLASS_DEFINED

#ifndef CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_2.h>
#endif // CH_NO_POSTCONDITIONS


#include <CGAL/ch_selected_extreme_points_2.h>
#include <list>
#include <algorithm>
#include <CGAL/stl_extensions.h>
#include <CGAL/ch_graham_andrew.h>

CGAL_BEGIN_NAMESPACE
/*{\Moptions
outfile=cgal_ch_I_b.man
}*/

/*{\Mtext [[\#include <CGAL/ch_bykat.h>]]
}*/

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat(InputIterator first, InputIterator last, 
             OutputIterator  result,
             const Traits& ch_traits);
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.\\
{\sc traits}: uses |Traits::Point_2|, |Traits::Less_signed_distance_to_line_2|,
|Traits::Leftturn_2|, and |Traits::Less_xy_2|.
}*/

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat_with_threshold(InputIterator first, InputIterator last, 
                             OutputIterator  result,
                             const Traits& ch_traits);
#ifdef CGAL_POINT_2_H
/*{\Moptions
outfile=cgal_ch_b.man
}*/

template <class InputIterator, class OutputIterator, class R>
inline
OutputIterator
ch__bykat(InputIterator first, InputIterator last, 
              OutputIterator  result,
              Point_2<R>* )
{
  return ch_bykat(first, last, result, R() );
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator
ch_bykat(InputIterator first, InputIterator last, OutputIterator  result)
{
  return ch__bykat( first, last, result, ch_value_type(first) );
}
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.
}*/

template <class InputIterator, class OutputIterator, class R>
inline
OutputIterator
ch__bykat_with_threshold(InputIterator first, InputIterator last, 
                              OutputIterator  result,
                              Point_2<R>* )
{
  return ch_bykat_with_threshold(first, last, result, 
                                 R() );
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator
ch_bykat_with_threshold(InputIterator first, InputIterator last, 
                             OutputIterator  result)
{
  return ch__bykat_with_threshold( first, last, result, 
                                   ch_value_type(first) );
}
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.
}*/
#endif // CGAL_POINT_2_H
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_bykat.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_BYKAT_H

