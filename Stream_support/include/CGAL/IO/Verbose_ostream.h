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

/*!
\ingroup PkgStreamSupportRef

The class `Verbose_ostream` can be used as an output stream. The stream
output operator `<<` is defined for any type. The class
`Verbose_ostream` stores in an internal state a stream and whether the
output is active or not. If the state is active, the stream output
operator `<<` uses the internal stream to output its argument. If
the state is inactive, nothing happens.

\cgalHeading{Example}

The class `Verbose_ostream` can be conveniently used to implement for
example the `is_valid()` member function for triangulations or
other complex data structures.

\code{.cpp}
bool is_valid( bool verbose = false, int level = 0) {
Verbose_ostream verr( verbose);
verr << "Triangulation::is_valid( level = " << level << ')' << endl;
verr << " Number of vertices = " << size_of_vertices() << endl;
// ...
}
\endcode
*/
class Verbose_ostream
{
  bool b;
  std::ostream* o;

public:

  /// \name Creation
  /// @{

  /*!
  creates an output stream with state set to `active` that writes to the stream `out`.
  */
  Verbose_ostream(bool active = false, std::ostream& out = std::cerr)
    : b(active), o(&out)
  {}

  /// @}

  bool verbose() const { return b; }
  void set_verbose(bool active) { b = active; }
  std::ostream& out() { return *o; }

  /// \name Operations
  /// @{

  /*!
  writes the object `t` into the stream `out`.
  */
  template < class T >
  Verbose_ostream& operator<<(const T& t) { CGAL__VERB(t); }

  /// @}

  Verbose_ostream& operator<<( std::ostream& (*f)(std::ostream&)) { CGAL__VERB(f); }
  Verbose_ostream& operator<<( std::ios& (*f)(std::ios&)) { CGAL__VERB(f); }

  Verbose_ostream& flush()
  {
    if(b)
      o->flush();
    return *this;
  }

  Verbose_ostream& put(char c)
  {
    if(b)
      o->put(c);
    return *this;
  }

  Verbose_ostream& write(const char* s, int n)
  {
    if(b)
      o->write(s, n);
    return *this;
  }
};

} // namespace CGAL

#endif // CGAL_IO_VERBOSE_OSTREAM_H
