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
//
// release       : 
// release_date  : 
//
// file          : src/io.C
// package       : iostream
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Mariette.Yvinec (Mariette.Yvinec@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_IO_C
#define CGAL_IO_C

#include <CGAL/basic.h>
#include <CGAL/IO/io.h>

CGAL_BEGIN_NAMESPACE

int IO::mode = std::ios::xalloc();


IO::Mode
get_mode(std::ios& i)
{
    return static_cast<IO::Mode>(i.iword(IO::mode));
}

IO::Mode
set_ascii_mode(std::ios& i)
{
    IO::Mode m = get_mode(i);
    i.iword(IO::mode) = IO::ASCII;
    return m;
}


IO::Mode
set_binary_mode(std::ios& i)
{
    IO::Mode m = get_mode(i);
    i.iword(IO::mode) = IO::BINARY;
    return m;
}


IO::Mode
set_pretty_mode(std::ios& i)
{
    IO::Mode m = get_mode(i);
    i.iword(IO::mode) = IO::PRETTY;
    return m;
}


IO::Mode
set_mode(std::ios& i, IO::Mode m)
{
    IO::Mode old = get_mode(i);
    i.iword(IO::mode) = m;
    return old;
}

bool
is_pretty(std::ios& i)
{
    return i.iword(IO::mode) == IO::PRETTY;
}

bool
is_ascii(std::ios& i)
{
    return i.iword(IO::mode) == IO::ASCII;
}


bool
is_binary(std::ios& i)
{
    return i.iword(IO::mode) == IO::BINARY;
}

CGAL_END_NAMESPACE

#endif // CGAL_IO_C
