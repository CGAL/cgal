// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//

// release       : $CGAL_Revision: CGAL-0.9-I-04 $
// release_date  : $CGAL_Date: 1997/12/15 $
//
// file          : include/CGAL/IO/io.h
// source        : web/io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_IO_H
#define CGAL_IO_H

#include <iostream.h>
#include <CGAL/IO/io_tags.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Object.h>

class CGAL_IO {
public:
    static int mode;
    enum Mode {ASCII = 0, PRETTY, BINARY};
};

CGAL_IO::Mode
CGAL_get_mode(ios& i);

CGAL_IO::Mode
CGAL_set_ascii_mode(ios& i);

CGAL_IO::Mode
CGAL_set_binary_mode(ios& i);

CGAL_IO::Mode
CGAL_set_pretty_mode(ios& i);

CGAL_IO::Mode
CGAL_set_mode(ios& i, CGAL_IO::Mode m);
bool
CGAL_is_pretty(ios& i);

bool
CGAL_is_ascii(ios& i);

bool
CGAL_is_binary(ios& i);

inline CGAL_io_Read_write CGAL_io_tag(char){ return CGAL_io_Read_write(); }


template < class T >
inline
void
CGAL_write(ostream& os, const T& t, const CGAL_io_Read_write&)
{
    os.write((char*)&t, sizeof(t));
}


template < class T >
inline
void
CGAL_write(ostream& os, const T& t, const CGAL_io_Operator&)
{
    os << t;
}


template < class T >
inline
void
CGAL_write(ostream& os, const T& t, const CGAL_io_Extract_insert&)
{
    CGAL_insert(os, t);
}


template < class T >
inline
void
CGAL_write(ostream& os, const T& t)
{
    CGAL_write(os, t, CGAL_io_tag(t));
}


template < class T >
inline
void
CGAL_read(istream& is, T& t, const CGAL_io_Read_write&)
{
    is.read((char*)&t, sizeof(t));
}


template < class T >
inline
void
CGAL_read(istream& is, T& t, const CGAL_io_Operator&)
{
    is >> t;
}


template < class T >
inline
void
CGAL_read(istream& is, T& t, const CGAL_io_Extract_insert&)
{
    CGAL_extract(is, t);
}


template < class T >
inline
void
CGAL_read(istream& is, T& t)
{
    CGAL_read(is, t, CGAL_io_tag(t));
}


inline
ostream& operator<<( ostream& out, const CGAL_Color& col)
{
    switch(out.iword(CGAL_IO::mode)) {
    case CGAL_IO::ASCII :
        return out << col.red() << ' ' << col.green() << ' ' << col.blue();
    case CGAL_IO::BINARY :
        CGAL_write(out, col.red());
        CGAL_write(out, col.green());
        CGAL_write(out, col.blue());
        return out;
    default:
        return out << "Color(" << col.red() << ", " << col.green() << ", "
                   << col.blue() << ')';
    }
}

inline
istream &operator>>(istream &is, CGAL_Color& col)
{
    int r, g, b;
    switch(is.iword(CGAL_IO::mode)) {
    case CGAL_IO::ASCII :
        is >> r >> g >> b;
        break;
    case CGAL_IO::BINARY :
        CGAL_read(is, r);
        CGAL_read(is, g);
        CGAL_read(is, b);
        break;
    default:
        cerr << "" << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;
    }
    col = CGAL_Color(r,g,b);
    return is;
}


#endif // CGAL_IO_H
