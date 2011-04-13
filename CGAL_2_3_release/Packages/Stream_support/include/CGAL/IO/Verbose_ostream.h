// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  :
//
// file          : IO/Verbose_ostream.h
// package       : Stream_support (2.4)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// A stream like output class for verbose output.
// ======================================================================

#ifndef CGAL_IO_VERBOSE_OSTREAM_H
#define CGAL_IO_VERBOSE_OSTREAM_H

#include <iostream>

CGAL_BEGIN_NAMESPACE

#define CGAL__VERB(x) if (b) *o << x; return *this

class Verbose_ostream {
    bool          b;
    std::ostream* o;
public:
    Verbose_ostream( bool active = false, std::ostream& out = std::cerr)
        : b(active), o(&out){}

    bool          verbose()           const { return b; }
    void          set_verbose(bool active)  { b = active; }
    std::ostream& out()                     { return *o; }

    template < class T >
    Verbose_ostream& operator<<(const T& t)
    { CGAL__VERB(t); }

    Verbose_ostream& operator<<( std::ostream& (*f)(std::ostream&))
    { CGAL__VERB(f); }

    Verbose_ostream& operator<<( std::ios& (*f)(std::ios&))
    { CGAL__VERB(f); }

    Verbose_ostream& flush()
    {
        if (b)
            o->flush();
        return *this;
    }

    Verbose_ostream& put(char c)
    {
        if (b)
            o->put(c);
        return *this;
    }

    Verbose_ostream& write(const char* s, int n)
    {
        if (b)
            o->write(s, n);
        return *this;
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_IO_VERBOSE_OSTREAM_H
