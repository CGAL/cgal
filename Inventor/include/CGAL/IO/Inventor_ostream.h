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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Andreas Fabri
//                 Lutz Kettner <kettner@inf.ethz.ch>
//                 Herve Bronnimann <Herve.Bronnimann@sophia.inria.fr>
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_IO_INVENTOR_OSTREAM_H
#define CGAL_IO_INVENTOR_OSTREAM_H

#include <CGAL/basic.h>
#include <iostream>

// OpenInventor and VRML 1.0 are quite similar formats, so
// output operators could be shared if they use the following
// base class, which is common for both output streams.

namespace CGAL {

class Inventor_ostream_base {
private:
    std::ostream*  m_os;
public:
    Inventor_ostream_base()           : m_os(0)  {}
    Inventor_ostream_base(std::ostream& o) : m_os(&o) {}
    ~Inventor_ostream_base()  { close(); }
    void open(std::ostream& o)        { m_os = &o; }
    void close() {
        if ( m_os)
            os() << std::endl;
        m_os = 0;
    }
    #if defined BOOST_NO_EXPLICIT_CONVERSION_OPERATORS \
        || defined BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
    typedef const void* Const_void_ptr;
    operator Const_void_ptr () const {
      if ( m_os->fail() )
        return 0;
      else
        return static_cast<Const_void_ptr>(m_os);
    }
    #else
    explicit operator bool ()
    {
      return m_os && !m_os->fail();
    }
    #endif
    std::ostream& os() {
        // The behaviour if m_os == 0 could be changed to return
        // cerr or a file handle to /dev/null. The latter one would
        // mimick the behaviour that one can still use a stream with
        // an invalid stream, but without producing any output.
        CGAL_assertion( m_os != 0 );
        return *m_os;
    }
};


class Inventor_ostream : public  Inventor_ostream_base
{
public:
    Inventor_ostream() {}
    Inventor_ostream(std::ostream& o) : Inventor_ostream_base(o) {
        header();
    }
    void open(std::ostream& o) {
        Inventor_ostream_base::open(o);
        header();
    }
private:
    void header() {
        os() << "#Inventor V2.0 ascii" << std::endl;
        os() << "# File written with the help of the CGAL Library" 
	     << std::endl;
    }
};

} //namespace CGAL
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

#endif // CGAL_INVENTOR_TETRAHEDRON_3
#endif // CGAL_TETRAHEDRON_3_H
