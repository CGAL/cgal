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
// file          : include/CGAL/ch_graham_andrew.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_GRAHAM_ANDREW_H
#define CGAL_CH_GRAHAM_ANDREW_H

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


#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE
/*{\Moptions
outfile=cgal_ch_I_gas.man
}*/

/*{\Mtext
\settowidth{\typewidth}{|OutputIterator|}
\addtolength{\typewidth}{\colsep}
\settowidth{\callwidth}{|ch_|}
\computewidths
}*/

/*{\Mtext [[\#include <CGAL/ch_graham_andrew.h>]]
}*/

template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
                            BidirectionalIterator last,
                            OutputIterator        result,
                            const Traits& ch_traits );
/*{\Mfuncl 
computes the sorted sequence of extreme points which are not left 
of $pq$ and reports this sequence in a range starting at |result|,
where $p$ is the value of |first| and $q$ is the value of |last| $-1$.
The sequence reported starts with $p$, point $q$ is omitted.\\
{\it Precondition:} The points in [|first|,|last|) are sorted with respect
to $pq$ and the range [|first|,|last|) contains at least two different 
points.\\
{\sc traits}: uses |Traits::Leftturn_2| operating on the 
point type |Traits::Point_2|.
}*/

template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch__ref_graham_andrew_scan( BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 OutputIterator&       result,
                                 const Traits&         ch_traits);

/*{\Moptions
outfile=cgal_ch_I_ga.man
}*/

/*{\Mtext [[\#include <CGAL/ch_graham_andrew.h>]]
}*/

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew( InputIterator  first,
                       InputIterator  last,
                       OutputIterator result,
                       const Traits&  ch_traits );
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.\\
{\sc traits}: uses |Traits::Point_2|, |Traits::Leftturn_2|
and |Traits::Less_xy_2|.
}*/

#ifdef CGAL_POINT_2_H

/*{\Moptions
outfile=cgal_ch_gas.man
}*/

/*{\Mtext 
\settowidth{\typewidth}{|OutputIterator|}
\addtolength{\typewidth}{\colsep}
\settowidth{\callwidth}{|ch_|}
\computewidths
}*/

template <class BidirectionalIterator, class OutputIterator, class R>
inline
OutputIterator
ch__graham_andrew_scan( BidirectionalIterator first,
                             BidirectionalIterator last,
                             OutputIterator        result,
                             Point_2<R>* )
{
  return ch_graham_andrew_scan( first, last, result, 
                                R() );
}

template <class BidirectionalIterator, class OutputIterator>
inline
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
                            BidirectionalIterator last,
                            OutputIterator        result )
{ 
  return ch__graham_andrew_scan( first, last, result, 
                                 ch_value_type(first) ); 
}
/*{\Mfuncl 
computes the sorted sequence of extreme points which are not left 
of $pq$ and reports this sequence in a range starting at |result|,
where $p$ is the value of |first| and $q$ is the value of |last| $-1$.
The sequence reported starts with $p$, point $q$ is omitted.\\
}*/

/*{\Moptions
outfile=cgal_ch_ga.man
}*/

template <class InputIterator, class OutputIterator, class R>
inline
OutputIterator
ch__graham_andrew( InputIterator     first,
                        InputIterator     last,
                        OutputIterator    result,
                        Point_2<R>* )
{
  return ch_graham_andrew(first, last, result, 
                               R());
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator
ch_graham_andrew( InputIterator  first,
                       InputIterator  last,
                       OutputIterator result )
{ return ch__graham_andrew( first, last, result, ch_value_type(first) ); }
/*{\Mfuncl 
same as |convex_hull_2(first,last,result)|.
}*/

#endif // CGAL_POINT_2_H
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_lower_hull_scan( InputIterator  first,
                         InputIterator  last,
                         OutputIterator result,
                         const Traits&  ch_traits);
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_upper_hull_scan( InputIterator  first,
                         InputIterator  last,
                         OutputIterator result,
                         const Traits&  ch_traits);
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_graham_andrew.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_GRAHAM_ANDREW_H

