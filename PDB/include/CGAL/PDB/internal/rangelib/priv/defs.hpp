// Boost. Iterable Range Library (rangelib)
//
// Copyright 2003-2004 John Torjo (john@torjo.com) and Matthew Wilson (matthew@synesis.com.au)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.


#ifndef CGAL_PDB_BOOST_RTL_DEFS_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_DEFS_HPP_INCLUDED

#include <boost/config.hpp>

#if defined(BOOST_MSVC) || \
	defined(BOOST_INTEL)
# if _MSC_VER <= 1200
#  define CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
# endif
#endif


/*
    Does the following compile?

    template< class T> void f( T val) {}
    template< class Res, class Arg> void f( Res (*func)(Arg) ) { }

    void test( int) {}

    int main(int argc, char* argv[]) {
        f( &test);
    }
*/
#define CGAL_PDB_BOOST_RTL_ALLOWS_OVERLOAD_BY_FUNC_PARAM


#ifdef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
// what else is new? ;)
#undef CGAL_PDB_BOOST_RTL_ALLOWS_OVERLOAD_BY_FUNC_PARAM
#endif


#endif
