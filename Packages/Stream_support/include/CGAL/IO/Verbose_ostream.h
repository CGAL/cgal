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
    void          set_verbose( bool active) { b = active; }
    std::ostream& out()                     { return *o; }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
    template < class T >
    Verbose_ostream&  operator<<( const T& t)            { CGAL__VERB(t);}
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES //
    // The following specialisations avoid the & for their small args.
    Verbose_ostream&  operator<<( char c)                { CGAL__VERB(c);}
    Verbose_ostream&  operator<<( const char* s)         { CGAL__VERB(s);}
    Verbose_ostream&  operator<<( int a)                 { CGAL__VERB(a);}
    Verbose_ostream&  operator<<( long l)                { CGAL__VERB(l);}
    Verbose_ostream&  operator<<( double d)              { CGAL__VERB(d);}
    Verbose_ostream&  operator<<( float f)               { CGAL__VERB(f);}
    Verbose_ostream&  operator<<( unsigned int a)        { CGAL__VERB(a);}
    Verbose_ostream&  operator<<( unsigned long l)       { CGAL__VERB(l);}
#ifdef CGAL_HAS_INTEGER64
    Verbose_ostream&  operator<<( Integer64 l)           { CGAL__VERB(l);}
    Verbose_ostream&  operator<<( UInteger64 l)          { CGAL__VERB(l);}
#endif
    Verbose_ostream&  operator<<( void* p)               { CGAL__VERB(p);}
    Verbose_ostream&  operator<<( short i)               { CGAL__VERB(i);}
    Verbose_ostream&  operator<<( unsigned short i)      { CGAL__VERB(i);}

    Verbose_ostream&  operator<<( std::ostream& (*f)(std::ostream&))
                                                         { CGAL__VERB(f);}
    Verbose_ostream&  operator<<( std::ios& (*f)(std::ios&))
                                                         { CGAL__VERB(f);}
    Verbose_ostream&  flush() {
        if (b)
            o->flush();
        return *this;
    }
    Verbose_ostream&  put(char c) {
        if (b)
            o->put(c);
        return *this;
    }
    Verbose_ostream&  write(const char*  s,int n) {
        if (b)
            o->write( s, n);
        return *this;
    }
};

#undef CGAL__VERB

CGAL_END_NAMESPACE

#endif // CGAL_IO_VERBOSE_OSTREAM_H
