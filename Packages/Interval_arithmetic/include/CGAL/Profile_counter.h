// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Profile_counter.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_PROFILE_COUNTER_H
#define CGAL_PROFILE_COUNTER_H

// This file contains the class Profile_counter which is able to keep track
// of a number, and prints a message in the destructor.
// Typically, it can be used as a profile counter in a static variable.

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

struct Profile_counter
{
    Profile_counter(const char *ss)
	: i(0), s(ss) {}

    void operator++() { ++i; }

    ~Profile_counter()
    {
	std::cerr << "Profile counter : " << s << " = " << i << std::endl;
    }

private:
    unsigned int i;
    const char *s;
};

#ifdef CGAL_PROFILE
#  define CGAL_PROFILER(X, Y) static CGAL::Profile_counter X(Y); ++X;
#else
#  define CGAL_PROFILER(X, Y)
#endif

CGAL_END_NAMESPACE

#endif // CGAL_PROFILE_COUNTER_H
