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
// Author(s)     : Andreas Fabri

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/assertions.h>

#include <sstream>
#include <string>

namespace CGAL {

CGAL_INLINE_FUNCTION
IO::Mode
get_mode(std::ios& i)
{
  return static_cast<IO::Mode>(i.iword(IO::Static::get_mode()));
}

CGAL_INLINE_FUNCTION
IO::Mode
set_ascii_mode(std::ios& i)
{
    IO::Mode m = get_mode(i);
    i.iword(IO::Static::get_mode()) = IO::ASCII;
    return m;
}

CGAL_INLINE_FUNCTION
IO::Mode
set_binary_mode(std::ios& i)
{
    IO::Mode m = get_mode(i);
    i.iword(IO::Static::get_mode()) = IO::BINARY;
    return m;
}


CGAL_INLINE_FUNCTION
IO::Mode
set_pretty_mode(std::ios& i)
{
    IO::Mode m = get_mode(i);
    i.iword(IO::Static::get_mode()) = IO::PRETTY;
    return m;
}

CGAL_INLINE_FUNCTION
IO::Mode
set_mode(std::ios& i, IO::Mode m)
{
    IO::Mode old = get_mode(i);
    i.iword(IO::Static::get_mode()) = m;
    return old;
}

CGAL_INLINE_FUNCTION
bool
is_pretty(std::ios& i)
{
    return i.iword(IO::Static::get_mode()) == IO::PRETTY;
}

CGAL_INLINE_FUNCTION
bool
is_ascii(std::ios& i)
{
    return i.iword(IO::Static::get_mode()) == IO::ASCII;
}

CGAL_INLINE_FUNCTION
bool
is_binary(std::ios& i)
{
    return i.iword(IO::Static::get_mode()) == IO::BINARY;
}

CGAL_INLINE_FUNCTION
const char*
mode_name( IO::Mode m) {
    static const char* const names[] = {"ASCII", "PRETTY", "BINARY" };
    CGAL_assertion( IO::ASCII <= m && m <= IO::BINARY );
    return names[m];
}

CGAL_INLINE_FUNCTION
void
swallow(std::istream &is, char d) {
    char c;
    do is.get(c); while (isspace(c));
    if (c != d) {
      std::stringstream msg;
      msg << "input error: expected '" << d << "' but got '" << c << "'";
      CGAL_error_msg( msg.str().c_str() );
    }
}

CGAL_INLINE_FUNCTION
void
swallow(std::istream &is, const std::string& s ) {
    std::string t;
    is >> t;
    if (s != t) {
      std::stringstream msg;
      msg << "input error: expected '" << s << "' but got '" << t << "'";
      CGAL_error_msg( msg.str().c_str() );
    }
}

} //namespace CGAL
