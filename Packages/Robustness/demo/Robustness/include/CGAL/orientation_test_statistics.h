// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : demo/Robustness/orientation_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/Interval_arithmetic.h>

// #include <sstream> // doesn't work with GCC 2.95.2
#include <strstream>

template <class ForwardIterator, class Traits>
void
orientation_statistics( ForwardIterator first, ForwardIterator last,
                        CGAL::Orientation* C,
                        leda_string& /*s1*/, leda_string& /*s2*/,
                        const Traits& T )
{
    typename Traits::Orientation_2  orientation = T.orientation_2_object();
    int c = 0;
    int success = 0;
    for( ForwardIterator i = first; i != last; ++i)
        for( ForwardIterator j = first; j != last; ++j)
            for( ForwardIterator k = first; k != last; ++k)
                if ( C[c++] == orientation(*i, *j, *k)) ++success;

    //s1 = leda_string("Out of %d orientation tests, %d", c, success);
    //s2 = leda_string("( %2.2f %%) give the correct result.",
    //                    (double)success/c * 100);
}

template <class ForwardIterator, class Traits>
void
fill_control_field( ForwardIterator first, ForwardIterator last,
                    CGAL::Orientation* C,
                    const Traits& T )
{
    typename Traits::Orientation_2  orientation = T.orientation_2_object();
    int c = 0;
    for( ForwardIterator i = first; i != last; ++i)
        for( ForwardIterator j = first; j != last; ++j)
            for( ForwardIterator k = first; k != last; ++k)
                C[c++] = orientation(*i, *j, *k);
}

template <class ForwardIterator, class Traits>
void
orientation_statistics_IA( ForwardIterator first, ForwardIterator last,
                           leda_string& s1, leda_string& s2,
                           const Traits& T )
{
    typename Traits::Orientation_2  orientation = T.orientation_2_object();
    int c = 0;
    int success = 0;
    for( ForwardIterator i = first; i != last; ++i)
        for( ForwardIterator j = first; j != last; ++j)
            for( ForwardIterator k = first; k != last; ++k)
                {
                  try
                  {
                    (void) orientation(*i, *j, *k);
                    ++success;
                  }
                  catch ( CGAL::Interval_base::unsafe_comparison )
                  {}
                  ++c;
                }

#if defined(CGAL_USE_CGAL_WINDOW)
    // std::ostringstream OS;
    std::ostrstream OS;
    OS << "Out of " << c << " orientation tests, " << c;
    OS << std::ends;
    s1 = OS.str();
    s2 = std::string("do not throw an exception.");
#else
    s1 = leda_string("Out of %d orientation tests, %d", c, success);
    s2 = leda_string("( %2.2f %%) do not throw an exception.",
                        (double)success/c * 100);
#endif
}
