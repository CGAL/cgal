// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_IO_H
#define CGAL_IO_H

#include <iostream>
#include <CGAL/IO/io_tags.h>
#include <CGAL/IO/Color.h>


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
    os.write(reinterpret_cast<const char*>(&t), sizeof(t));
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
    is.read(reinterpret_cast<char*>(&t), sizeof(t));
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
        return out << static_cast<int>(col.red())   << ' ' 
		   << static_cast<int>(col.green()) << ' ' 
		   << static_cast<int>(col.blue());
    case IO::BINARY :
        write(out, static_cast<int>(col.red()));
        write(out, static_cast<int>(col.green()));
        write(out, static_cast<int>(col.blue()));
        return out;
    default:
        return out << "Color(" << static_cast<int>(col.red()) << ", " 
		   << static_cast<int>(col.green()) << ", "
                   << static_cast<int>(col.blue()) << ')';
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
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    col = Color(r,g,b);
    return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_IO_H
