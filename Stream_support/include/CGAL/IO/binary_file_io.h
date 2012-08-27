// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>


#ifndef CGAL_IO_BINARY_FILE_IO_H
#define CGAL_IO_BINARY_FILE_IO_H

#include <CGAL/basic.h>
#include <iostream>
#include <boost/cstdint.hpp>

namespace CGAL {

inline void
I_Binary_write_integer32(std::ostream& out, boost::int32_t i) {
    out.write( (char*)(&i), 4);
}
inline void
I_Binary_write_float32(std::ostream& out, float f) {
    out.write( (char*)(&f), 4);
}

inline void
I_Binary_read_integer32(std::istream& in, boost::int32_t& i) {
    in.read( (char*)(&i), 4);
}
inline void
I_Binary_read_float32(std::istream& in, float& f) {
    in.read( (char*)(&f), 4);
}

inline void
I_swap_to_big_endian( boost::uint32_t& u) {
    (void) u;
#ifdef CGAL_LITTLE_ENDIAN
    u = ((u >> 24) | (u << 24) | ((u >> 8) & 0xff00) | ((u << 8) & 0xff0000));
#endif
}

inline void
I_swap_to_big_endian( boost::int32_t& i) {
    // We need to use a union instead of the 2 lines below,
    // otherwise we get aliasing issues.
    // boost::uint32_t& u = (boost::uint32_t&)i;
    // I_swap_to_big_endian( u);
    union {
      boost::int32_t  in;
      boost::uint32_t ui;
    } u;
    u.in = i;
    I_swap_to_big_endian(u.ui);
    i = u.in;
}

inline void
I_swap_to_big_endian( float& f) {
    // We need to use a union instead of the 2 lines below,
    // otherwise we get aliasing issues.
    // boost::uint32_t& u = (boost::uint32_t&)f;
    // I_swap_to_big_endian( u);
    union {
      boost::uint32_t ui;
      float           fl;
    } u;
    u.fl = f;
    I_swap_to_big_endian(u.ui);
    f = u.fl;
}

inline void
I_Binary_write_big_endian_integer32(std::ostream& out, boost::int32_t i) {
    I_swap_to_big_endian( i);
    out.write( (char*)(&i), 4);
}
inline void
I_Binary_write_big_endian_float32(std::ostream& out, float f) {
    I_swap_to_big_endian( f);
    out.write( (char*)(&f), 4);
}

inline void
I_Binary_read_big_endian_integer32(std::istream& in, boost::int32_t& i) {
    in.read( (char*)(&i), 4);
    I_swap_to_big_endian( i);
}
inline void
I_Binary_read_big_endian_float32(std::istream& in, float& f) {
    in.read( (char*)(&f), 4);
    I_swap_to_big_endian( f);
}

} //namespace CGAL

#endif // CGAL_IO_BINARY_FILE_IO_H
