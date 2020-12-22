// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>


#ifndef CGAL_IO_BINARY_FILE_IO_H
#define CGAL_IO_BINARY_FILE_IO_H

#include <CGAL/config.h>
#include <CGAL/assertions.h>
#include <iostream>
#include <cstdint>
#include <limits>

namespace CGAL {

inline void
I_Binary_write_uinteger32(std::ostream& out, std::uint32_t u) {
    out.write( (char*)(&u), 4);
}
// Special function to write size_t in 32b integer to ensure files
// written by 64b systems are still readable by 32b ones
inline void
I_Binary_write_size_t_into_uinteger32 (std::ostream& out, std::size_t s) {
    CGAL_assertion_msg
      (s <= static_cast<std::size_t>((std::numeric_limits<std::uint32_t>::max)()),
       "Trying to write size_t that does not fit in uint32_t");
    I_Binary_write_uinteger32 (out, static_cast<std::uint32_t>(s));
}
inline void
I_Binary_write_integer32(std::ostream& out, std::int32_t i) {
    out.write( (char*)(&i), 4);
}
inline void
I_Binary_write_float32(std::ostream& out, float f) {
    out.write( (char*)(&f), 4);
}
inline void
I_Binary_write_bool(std::ostream& out, bool b) {
    char c = (b ? 1 : 0);
    out.write(&c, 1);
}

// Special function to read size_t from 32b integer to ensure files
inline void
I_Binary_read_uinteger32(std::istream& in, std::uint32_t& u) {
    in.read( (char*)(&u), 4);
}
// written by 64b systems are still readable by 32b ones
inline void
I_Binary_read_size_t_from_uinteger32(std::istream& in, std::size_t& s) {
    std::uint32_t s32;
    I_Binary_read_uinteger32 (in, s32);
    s = static_cast<std::size_t>(s32);
}
inline void
I_Binary_read_integer32(std::istream& in, std::int32_t& i) {
    in.read( (char*)(&i), 4);
}
inline void
I_Binary_read_float32(std::istream& in, float& f) {
    in.read( (char*)(&f), 4);
}
inline void
I_Binary_read_bool(std::istream& in, bool& b) {
    char c;
    in.read(&c, 1);
    b = (c != 0);
}

inline void
I_swap_to_big_endian( std::uint32_t& u) {
    (void) u;
#ifdef CGAL_LITTLE_ENDIAN
    u = ((u >> 24) | (u << 24) | ((u >> 8) & 0xff00) | ((u << 8) & 0xff0000));
#endif
}

inline void
I_swap_to_big_endian( std::int32_t& i) {
    // We need to use a union instead of the 2 lines below,
    // otherwise we get aliasing issues.
    // std::uint32_t& u = (std::uint32_t&)i;
    // I_swap_to_big_endian( u);
    union {
      std::int32_t  in;
      std::uint32_t ui;
    } u;
    u.in = i;
    I_swap_to_big_endian(u.ui);
    i = u.in;
}

inline void
I_swap_to_big_endian( float& f) {
    // We need to use a union instead of the 2 lines below,
    // otherwise we get aliasing issues.
    // std::uint32_t& u = (std::uint32_t&)f;
    // I_swap_to_big_endian( u);
    union {
      std::uint32_t ui;
      float           fl;
    } u;
    u.fl = f;
    I_swap_to_big_endian(u.ui);
    f = u.fl;
}

inline void
I_Binary_write_big_endian_integer32(std::ostream& out, std::int32_t i) {
    I_swap_to_big_endian( i);
    out.write( (char*)(&i), 4);
}
inline void
I_Binary_write_big_endian_float32(std::ostream& out, float f) {
    I_swap_to_big_endian( f);
    out.write( (char*)(&f), 4);
}

inline void
I_Binary_read_big_endian_integer32(std::istream& in, std::int32_t& i) {
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
