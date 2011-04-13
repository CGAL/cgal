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
// file          : include/CGAL/ch_selected_extreme_points_2.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_SELECTED_EXTREME_POINTS_2_H
#define CGAL_CH_SELECTED_EXTREME_POINTS_2_H

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



CGAL_BEGIN_NAMESPACE

/*{\Moptions
outfile=cgal_ch_I_sep.man
}*/

/*{\Mtext
\settowidth{\typewidth}{|OutputIterator|}
\addtolength{\typewidth}{\colsep}
\settowidth{\callwidth}{|ch_|}
\computewidths
}*/

/*{\Mtext [[\#include <CGAL/ch_selected_extreme_points_2.h>]]
}*/

template <class ForwardIterator, class Traits>
void
ch_nswe_point( ForwardIterator first, ForwardIterator last,
                    ForwardIterator& n,
                    ForwardIterator& s,
                    ForwardIterator& w,
                    ForwardIterator& e,
                    const Traits& ch_traits);
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |n| is
an iterator in the range such that |*n| $\ge_{\rm yx}$ |*it| for
all iterators |it| in the range. Similarly, for |s|, |w|, and |e| 
the inequations |*s| $\le_{\rm yx}$ |*it|, |*w| $\le_{\rm xy}$ |*it|,
and |*e| $\ge_{\rm yx}$ |*it| hold respectively for all iterators 
|it| in the range.\\
{\sc traits}: uses |Traits::Less_xy_2| and |Traits::Less_yx_2|.
}*/


/*{\Moptions
outfile=cgal_ch_I_sep2.man
}*/

/*{\Mtext
\settowidth{\typewidth}{|OutputIterator|}
\addtolength{\typewidth}{\colsep}
\settowidth{\callwidth}{|ch_|}
\computewidths
}*/

template <class ForwardIterator, class Traits>
void
ch_ns_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& n,
                  ForwardIterator& s,
                  const Traits& ch_traits );
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |n| is
an iterator in the range such that |*n| $\ge_{\rm yx}$ |*it| for
all iterators |it| in the range. Similarly, for |s| the inequation
|*s| $\le_{\rm yx}$ |*it| holds for all iterators |it| in the range.\\
{\sc traits}: uses function object type |Traits::Less_yx_2|.
}*/

template <class ForwardIterator, class Traits>
void
ch_we_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& w,
                  ForwardIterator& e,
                  const Traits& ch_traits );
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |w| is
an iterator in the range such that |*w| $\le_{\rm xy}$ |*it| for
all iterators |it| in the range. Similarly, for |e| the inequation
|*e| $\ge_{\rm yx}$ |*it| holds for all iterators |it| in the range.\\
{\sc traits}: uses function object type |Traits::Less_xy_2|.
}*/


template <class ForwardIterator, class Traits>
void
ch_n_point( ForwardIterator first, ForwardIterator last,
                 ForwardIterator& n,
                 const Traits& ch_traits );
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |n| is
an iterator in the range such that |*n| $\ge_{\rm yx}$ |*it| for
all iterators |it| in the range.\\ 
{\sc traits}: uses |Traits::Less_yx_2|.
}*/

template <class ForwardIterator, class Traits>
void
ch_s_point( ForwardIterator first, ForwardIterator last,
                 ForwardIterator& s,
                 const Traits& ch_traits );
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |s| is
an iterator in the range such that |*s| $\le_{\rm yx}$ |*it| for
all iterators |it| in the range.\\ 
{\sc traits}: uses |Traits::Less_yx_2|.
}*/

template <class ForwardIterator, class Traits>
void
ch__e_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& e,
                  const Traits& ch_traits );
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |e| is
an iterator in the range such that |*e| $\ge_{\rm yx}$ |*it| for
for all iterators |it| in the range.\\
{\sc traits}: uses |Traits::Less_xy_2|.
}*/

template <class ForwardIterator, class Traits>
void
ch_w_point( ForwardIterator first, ForwardIterator last,
                 ForwardIterator& w,
                 const Traits& ch_traits );
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |w| is
an iterator in the range such that |*w| $\le_{\rm yx}$ |*it| for
for all iterators |it| in the range.\\
{\sc traits}: uses |Traits::Less_yx_2|.
}*/
#ifdef CGAL_POINT_2_H
/*{\Moptions
outfile=cgal_ch_sep.man
}*/

/*{\Mtext 
\settowidth{\typewidth}{|OutputIterator|}
\addtolength{\typewidth}{\colsep}
\settowidth{\callwidth}{|ch_|}
\computewidths
}*/

template <class ForwardIterator, class R>
inline
void
ch__nswe_point( ForwardIterator first, ForwardIterator last,
                     ForwardIterator& n,
                     ForwardIterator& s,
                     ForwardIterator& w,
                     ForwardIterator& e,
                     Point_2<R>* )
{ ch_nswe_point(first, last, n, s, w, e, R()); }

template <class ForwardIterator>
inline
void
ch_nswe_point( ForwardIterator first, ForwardIterator last,
                    ForwardIterator& n,
                    ForwardIterator& s,
                    ForwardIterator& w,
                    ForwardIterator& e)
{ ch__nswe_point(first, last, n, s, w, e, ch_value_type(first) ); }
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |n| is
an iterator in the range such that |*n| $\ge_{\rm yx}$ |*it| for
all iterators |it| in the range. Similarly, for |s|, |w|, and |e| 
the inequations |*s| $\le_{\rm yx}$ |*it|, |*w| $\le_{\rm xy}$ |*it|,
and |*e| $\ge_{\rm yx}$ |*it| hold respectively for all iterators 
|it| in the range.}*/


/*{\Moptions
outfile=cgal_ch_sep.man
}*/

/*{\Mtext 
\settowidth{\typewidth}{|OutputIterator|}
\addtolength{\typewidth}{\colsep}
\settowidth{\callwidth}{|ch_|}
\computewidths
}*/

template <class ForwardIterator, class R>
inline
void
ch__n_point( ForwardIterator first, ForwardIterator last,
                 ForwardIterator& n, 
                 Point_2<R>* )
{
  ch_n_point(first, last, n, R() );
}

template <class ForwardIterator>
inline
void
ch_n_point( ForwardIterator first, ForwardIterator last,
                 ForwardIterator& n)
{
  ch__n_point(first, last, n, ch_value_type(first) );
}
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |n| is
an iterator in the range such that |*n| $\ge_{\rm yx}$ |*it| for
all iterators |it| in the range. 
}*/

template <class ForwardIterator, class R>
inline
void
ch__s_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& s,
                  Point_2<R>* )
{
  ch_s_point(first, last, s, R() );
}

template <class ForwardIterator>
inline
void
ch_s_point( ForwardIterator first, ForwardIterator last,
                 ForwardIterator& s)
{
  ch__s_point(first, last, s, ch_value_type(first) );
}
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |s| is
an iterator in the range such that |*s| $\le_{\rm yx}$ |*it| for
all iterators |it| in the range. 
}*/

template <class ForwardIterator, class R>
inline
void
ch__e_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& e,
                  Point_2<R>* )
{
  ch_e_point(first, last, e, R() );
}

template <class ForwardIterator>
inline
void
ch_e_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& e)
{
  ch__e_point(first, last, e, ch_value_type(first) );
}
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |e| is
an iterator in the range such that |*e| $\ge_{\rm xy}$ |*it| for
all iterators |it| in the range. 
}*/

template <class ForwardIterator, class R>
inline
void
ch__w_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& w,
                  Point_2<R>* )
{
  ch_w_point(first, last, w, R() );
}

template <class ForwardIterator>
inline
void
ch_w_point( ForwardIterator first, ForwardIterator last,
                      ForwardIterator& w)
{
 ch__w_point(first, last, w, ch_value_type(first) );
}
/*{\Mfuncl traverses the range [|first|,|last|). After execution, 
the value of |w| is
an iterator in the range such that |*w| $\le_{\rm xy}$ |*it| for
all iterators |it| in the range. 
}*/

template <class ForwardIterator, class R>
inline
void
ch__ns_point( ForwardIterator first, ForwardIterator last,
                   ForwardIterator& n,
                   ForwardIterator& s,
                   Point_2<R>* )
{
  ch_ns_point(first, last, n, s, R() );
}

template <class ForwardIterator>
inline
void
ch_ns_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& n,
                  ForwardIterator& s)
{
  ch__ns_point(first, last, n, s, ch_value_type(first) );
}
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |n| is
an iterator in the range such that |*n| $\ge_{\rm yx}$ |*it| for
all iterators |it| in the range. Similarly, for |s| the inequation
|*s| $\le_{\rm yx}$ |*it| holds for all iterators |it| in the range.}*/

template <class ForwardIterator, class R>
inline
void
ch__we_point( ForwardIterator first, ForwardIterator last,
                   ForwardIterator& w,
                   ForwardIterator& e,
                   Point_2<R>* )
{
  ch_we_point(first, last, w, e, R() );
}

template <class ForwardIterator>
inline
void
ch_we_point( ForwardIterator first, ForwardIterator last,
                  ForwardIterator& w,
                  ForwardIterator& e)
{
  ch__we_point(first, last, w, e, ch_value_type(first) );
}
/*{\Mfuncl traverses the range [|first|,|last|). 
After execution, the value of |w| is
an iterator in the range such that |*w| $\le_{\rm xy}$ |*it| for
all iterators |it| in the range. Similarly, for |e| the inequation
|*e| $\ge_{\rm yx}$ |*it| holds for all iterators |it| in the range.}*/
#endif // CGAL_POINT_2_H
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_selected_extreme_points_2.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_SELECTED_EXTREME_POINTS_2_H

