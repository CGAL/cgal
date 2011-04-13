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
// file          : include/CGAL/ch_eddy.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_EDDY_H
#define CGAL_CH_EDDY_H

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

CGAL_BEGIN_NAMESPACE
/*{\Moptions
outfile=cgal_ch_I_e.man
}*/

/*{\Mtext [[\#include <CGAL/ch_eddy.h>]]
}*/

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_eddy(InputIterator first, InputIterator last, 
             OutputIterator  result,
             const Traits& ch_traits);
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.\\
{\sc traits}: uses |Traits::Point_2|, |Traits::Less_signed_distance_to_line_2|,
|Traits::Leftturn_2|, and |Traits::Less_xy_2|.
}*/

#ifdef CGAL_POINT_2_H
/*{\Moptions
outfile=cgal_ch_e.man
}*/

template <class InputIterator, class OutputIterator, class R>
inline
OutputIterator
ch__eddy(InputIterator first, InputIterator last, 
              OutputIterator  result,
              Point_2<R>* )
{
  return ch_eddy(first, last, result, R() );
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator
ch_eddy(InputIterator first, InputIterator last, OutputIterator  result)
{
  return ch__eddy( first, last, result, ch_value_type(first) );
}
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.
}*/

#endif // CGAL_POINT_2_H
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_eddy.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_EDDY_H

