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

// release       : 
// release_date  : 
//
// file          : include/CGAL/IO/io.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_IO_H
#define CGAL_IO_H

#include <iostream>
#include <CGAL/IO/io_tags.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

class IO {
public:
    static int mode;
    enum Mode {ASCII = 0, PRETTY, BINARY};
};

IO::Mode
get_mode(std::ios& i);

IO::Mode
set_ascii_mode(std::ios& i);

IO::Mode
set_binary_mode(std::ios& i);

IO::Mode
set_pretty_mode(std::ios& i);

IO::Mode
set_mode(std::ios& i, IO::Mode m);
bool
is_pretty(std::ios& i);

bool
is_ascii(std::ios& i);

bool
is_binary(std::ios& i);

inline io_Read_write io_tag(char){ return io_Read_write(); }


template < class T >
inline
void
write(std::ostream& os, const T& t, const io_Read_write&)
{
    os.write((char*)&t, sizeof(t));
}


template < class T >
inline
void
write(std::ostream& os, const T& t, const io_Operator&)
{
    os << t;
}


template < class T >
inline
void
write(std::ostream& os, const T& t, const io_Extract_insert&)
{
    insert(os, t);
}


template < class T >
inline
void
write(std::ostream& os, const T& t)
{
    write(os, t, io_tag(t));
}


template < class T >
inline
void
read(std::istream& is, T& t, const io_Read_write&)
{
    is.read((char*)&t, sizeof(t));
}


template < class T >
inline
void
read(std::istream& is, T& t, const io_Operator&)
{
    is >> t;
}


template < class T >
inline
void
read(std::istream& is, T& t, const io_Extract_insert&)
{
    extract(is, t);
}


template < class T >
inline
void
read(std::istream& is, T& t)
{
    read(is, t, io_tag(t));
}


inline
std::ostream& operator<<( std::ostream& out, const Color& col)
{
    switch(out.iword(IO::mode)) {
    case IO::ASCII :
        return out << col.red() << ' ' << col.green() << ' ' << col.blue();
    case IO::BINARY :
        write(out, col.red());
        write(out, col.green());
        write(out, col.blue());
        return out;
    default:
        return out << "Color(" << col.red() << ", " << col.green() << ", "
                   << col.blue() << ')';
    }
}

inline
std::istream &operator>>(std::istream &is, Color& col)
{
    int r, g, b;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> r >> g >> b;
        break;
    case IO::BINARY :
        read(is, r);
        read(is, g);
        read(is, b);
        break;
    default:
        cerr << "" << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;
    }
    col = Color(r,g,b);
    return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_IO_H
