// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Lutz Kettner <kettner@inf.ethz.ch>
//                 Herve Bronnimann <Herve.Bronnimann@sophia.inria.fr>
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_IO_VRML_INVENTOR_OSTREAM_H
#define CGAL_IO_VRML_INVENTOR_OSTREAM_H

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/boost/graph/named_params_helper.h>

// OpenInventor and VRML 1.0 are quite similar formats, so
// output operators could be shared if they use the following
// base class, which is common for both output streams.

namespace CGAL {

class Inventor_ostream_base
{
private:
  std::ostream* m_os = nullptr;

public:
  Inventor_ostream_base() : m_os(nullptr)  {}
  Inventor_ostream_base(std::ostream& o) : m_os(&o) {}
  ~Inventor_ostream_base() { close(); }
  void open(std::ostream& o) { m_os = &o; }
  void close()
  {
    if ( m_os)
      os() << std::endl;
    m_os = 0;
  }

  bool fail() const { return m_os->fail(); }
  bool good() const { return m_os->good(); }

  void precision(const int p) { m_os->precision(p); }
  std::streamsize precision() const { return m_os->precision(); }

  explicit operator bool()
  {
    return m_os && !m_os->fail();
  }

  std::ostream& os() const
  {
    // The behavior if m_os == nullptr could be changed to return
    // cerr or a file handle to /dev/null. The latter one would
    // mimic the behavior that one can still use a stream with
    // an invalid stream, but without producing any output.
    CGAL_assertion( m_os != nullptr );
    return *m_os;
  }
};

class Inventor_ostream
  : public Inventor_ostream_base
{
public:
  Inventor_ostream() {}
  Inventor_ostream(std::ostream& o) : Inventor_ostream_base(o) { header(); }

  void open(std::ostream& o)
  {
    Inventor_ostream_base::open(o);
    header();
  }

private:
  void header()
  {
    os() << "#Inventor V2.0 ascii" << std::endl;
    os() << "# File written with the help of the CGAL Library"
         << std::endl;
  }
};

template<typename NP>
void set_stream_precision_from_NP(Inventor_ostream_base& os, const NP& np)
{
  return set_stream_precision_from_NP(os.os(), np);
}

} // namespace CGAL

#endif // CGAL_IO_INVENTOR_OSTREAM_H

#ifdef CGAL_TETRAHEDRON_3_H
#ifndef CGAL_INVENTOR_TETRAHEDRON_3
#define CGAL_INVENTOR_TETRAHEDRON_3

namespace CGAL {

template <class R >
Inventor_ostream&
operator<<(Inventor_ostream& os,
           const Tetrahedron_3<R > &t)
{
  const char *Indent = "   ";
  os.os() << "\n Separator {";
  os.os() << "\n   Coordinate3 { \n"
          << Indent << "point [\n"
          << Indent << "  "
          << CGAL::to_double(t[0].x()) << " "
          << CGAL::to_double(t[0].y()) << " "
          << CGAL::to_double(t[0].z()) << " ,\n"
          << Indent << "  "
          << CGAL::to_double(t[1].x()) << " "
          << CGAL::to_double(t[1].y()) << " "
          << CGAL::to_double(t[1].z()) << " ,\n"
          << Indent << "  "
          << CGAL::to_double(t[2].x()) << " "
          << CGAL::to_double(t[2].y()) << " "
          << CGAL::to_double(t[2].z()) << " ,\n"
          << Indent << "  "
          << CGAL::to_double(t[3].x()) << " "
          << CGAL::to_double(t[3].y()) << " "
          << CGAL::to_double(t[3].z()) << " ]"
          << "\n   } #Coordinate3" ;
  os.os() << "\n   IndexedFaceSet {"
          << Indent << "coordIndex  [ 0,1,2,-1, 1,3,2,-1,\n"
          << Indent << "              0,2,3,-1, 0,3,1,-1 ]\n"
          << "\n   } #IndexedFaceSet"
          << "\n } #Separator\n";
  return os;
}

} //namespace CGAL

#endif // CGAL_TETRAHEDRON_3_H
#endif // CGAL_VRML_INVENTOR_TETRAHEDRON_3
