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
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_IO_VERBOSE_OSTREAM_H
#define CGAL_IO_VERBOSE_OSTREAM_H

#include <iostream>

namespace CGAL {

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

} //namespace CGAL

#endif // CGAL_IO_VERBOSE_OSTREAM_H
