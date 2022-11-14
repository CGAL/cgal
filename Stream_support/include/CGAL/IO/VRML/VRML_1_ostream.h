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

#ifndef CGAL_IO_VRML_VRML_1_OSTREAM_H
#define CGAL_IO_VRML_VRML_1_OSTREAM_H

// Declare the common base class for OpenInventor and VRML 1.0 format.
#include <CGAL/IO/OI.h>

#include <iostream>

// OpenInventor and VRML 1.0 are quite similar formats, so
// output operators could be shared if they use the common
// base class Inventor_ostream_base, which is common for
// both output streams.

namespace CGAL {

class VRML_1_ostream
  : public Inventor_ostream_base
{
public:
  VRML_1_ostream() {}
  VRML_1_ostream(std::ostream& o) : Inventor_ostream_base(o) { header(); }
  void open(std::ostream& o) {
    Inventor_ostream_base::open(o);
    header();
  }

private:
  void header()
  {
    os() << "#VRML V1.0 ascii" << std::endl;
    os() << "# File written with the help of the CGAL Library"
         << std::endl;
  }
};

} //namespace CGAL

#endif // CGAL_IO_VRML_1_OSTREAM_H

#ifdef CGAL_TETRAHEDRON_3_H
#ifndef CGAL_IO_VRML_1_TETRAHEDRON_3
#define CGAL_IO_VRML_1_TETRAHEDRON_3

namespace CGAL {

template <class R >
VRML_1_ostream&
operator<<(VRML_1_ostream& os,
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
#endif // CGAL_IO_VRML_VRML_1_TETRAHEDRON_3
