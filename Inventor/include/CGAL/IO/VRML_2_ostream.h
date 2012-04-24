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
// Author(s)     : Andreas Fabri
//                 Lutz Kettner <kettner@inf.ethz.ch>
//                 Herve Bronnimann
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_IO_VRML_2_OSTREAM_H
#define CGAL_IO_VRML_2_OSTREAM_H

#include <CGAL/basic.h>
#include <iostream>

namespace CGAL {

class VRML_2_ostream {
public:
    VRML_2_ostream()           : m_os(0)  {}
    VRML_2_ostream(std::ostream& o) : m_os(&o) { header();}
    ~VRML_2_ostream()  { close(); }
    void open(std::ostream& o)   { m_os = &o; header(); }
    void close() {
        if ( m_os)
            footer();
        m_os = 0;
    }
    typedef const void* Const_void_ptr;
    operator Const_void_ptr () const {
        if ( m_os)
            return *m_os;
        return 0;
    }
    std::ostream& os() {
        // The behaviour if m_os == 0 could be changed to return
        // cerr or a file handle to /dev/null. The latter one would
        // mimick the behaviour that one can still use a stream with
        // an invalid stream, but without producing any output.
        CGAL_assertion( m_os != 0 );
        return *m_os;
    }
private:
    void header() {
        os() << "#VRML V2.0 utf8\n"
                "# File written with the help of the CGAL Library\n"
                "#-- Begin of file header\n"
                "Group {\n"
                "    children [\n"
                "        Shape {\n"
                "          appearance DEF A1 Appearance {\n"
                "            material Material {\n"
                "              diffuseColor .6 .5 .9\n"
                "            }\n         }\n"
                "            appearance\n"
                "                Appearance {\n"
                "                    material DEF Material Material {}\n"
                "                }\n"
                "            geometry NULL\n"
                "        }\n"
                "        #-- End of file header" << std::endl;
    }
    void footer() {
        os() << "        #-- Begin of file footer\n"
                "    ]\n"
                "}\n"
                "#-- End of file footer" << std::endl;
    }
    std::ostream*  m_os;
};

inline
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const char* s)
{
  os.os() << s;
  return os;
}

inline
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const double& d)
{
  os.os() << d;
  return os;
}

} //namespace CGAL

#endif // CGAL_IO_VRML_2_OSTREAM_H

#ifdef CGAL_TETRAHEDRON_3_H
#ifndef CGAL_IO_VRML_2_TETRAHEDRON_3
#define CGAL_IO_VRML_2_TETRAHEDRON_3

namespace CGAL {

template <class R >
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const Tetrahedron_3<R > &t)
{
  const char *Indent = "                                    ";
  os <<      "        Group {\n"
             "            children [\n"
             "                Shape {\n"
             "                    appearance\n"
             "                        Appearance {\n"
             "                            material USE Material\n"
             "                        } #Appearance\n"
             "                    geometry\n"
             "                        IndexedFaceSet {\n"
             "                            coord Coordinate {\n"
             "                                point [ \n"
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
          << CGAL::to_double(t[3].z()) <<
             "\n                                ]\n"
             "                            }\n"
             "                            solid   FALSE\n"
          << Indent << "coordIndex  [ 0,1,2,-1, 1,3,2,-1,\n"
          << Indent << "              0,2,3,-1, 0,3,1,-1 ]\n"
             "                        } #IndexedFaceSet\n"
             "                } #Shape\n"
             "            ] #children\n"
             "        } #Group\n";
  return os;
}

} //namespace CGAL

#endif // CGAL_IO_VRML_2_TETRAHEDRON_3
#endif // CGAL_TETRAHEDRON_3_H

#ifdef CGAL_POINT_3_H
#ifndef CGAL_IO_VRML_2_POINT_3
#define CGAL_IO_VRML_2_POINT_3

namespace CGAL {

template <class R >
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const Point_3<R > &p)
{
  const char *Indent = "                                    ";
  os <<      "        Group {\n"
             "            children [\n"
             "                Shape {\n"
             "                    appearance USE A1\n"
             "                    geometry\n"
             "                        PointSet {\n"
             "                            coord Coordinate {\n"
             "                                point [ ";
  os << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) 
     << " " << CGAL::to_double(p.z()) << " ]\n";
  os << Indent << "}\n";
  os << Indent << "} # PointSet\n";
  os << "                } #Shape\n"
        "            ] #children\n"
        "        } #Group\n";
  return os;
}

} //namespace CGAL

#endif // CGAL_IO_VRML_2_POINT_3
#endif // CGAL_POINT_3_H



#ifdef CGAL_TRIANGLE_3_H
#ifndef CGAL_IO_VRML_2_TRIANGLE_3
#define CGAL_IO_VRML_2_TRIANGLE_3

namespace CGAL {

template <class R >
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const Triangle_3<R > &t)
{
  const char *Indent = "                                    ";
  os <<      "        Group {\n"
             "            children [\n"
             "                Shape {\n"
             "                    appearance USE A1\n"
             "                    geometry\n"
             "                        IndexedLineSet {\n"
             "                            coord Coordinate {\n"
             "                                point [ \n";
  os << Indent ;
  os << CGAL::to_double(t[0].x()) << " " << CGAL::to_double(t[0].y()) 
     << " " << CGAL::to_double(t[0].z()) << ",\n";
  os << Indent;
  os << CGAL::to_double(t[1].x()) << " " << CGAL::to_double(t[1].y()) 
     << " " << CGAL::to_double(t[1].z()) << ",\n";
  os << Indent;
  os << CGAL::to_double(t[2].x()) << " " << CGAL::to_double(t[2].y()) 
     << " " << CGAL::to_double(t[2].z()) << " ]\n";
  os << Indent << "}\n" << Indent << "coordIndex [ 0 1, 1 2, 2 0 -1 ]\n";
  os << Indent << "} # IndexedLineSet\n";
  os << "                } #Shape\n"
        "            ] #children\n"
        "        } #Group\n";
  return os;
}

} //namespace CGAL

#endif // CGAL_IO_VRML_2_TRIANGLE_3
#endif // CGAL_TRIANGLE_3_H


#ifdef CGAL_SEGMENT_3_H
#ifndef CGAL_IO_VRML_2_SEGMENT_3
#define CGAL_IO_VRML_2_SEGMENT_3

namespace CGAL {

template <class R >
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const Segment_3<R > &s)
{
  const char *Indent = "                                    ";
  os <<      "        Group {\n"
             "            children [\n"
             "                Shape {\n"
             "                    appearance USE A1\n"
             "                    geometry\n"
             "                        IndexedLineSet {\n"
             "                            coord Coordinate {\n"
             "                                point [ \n";
  os << Indent << CGAL::to_double(s.source().x());
  os << " " << CGAL::to_double(s.source().y()) 
     << " " << CGAL::to_double(s.source().z()) << ",\n";
  os << Indent;
  os << CGAL::to_double(s.target().x())
     << " " << CGAL::to_double(s.target().y()) 
     << " " << CGAL::to_double(s.target().z()) << " ]\n";
  os << Indent << "}\n" << Indent << "coordIndex [ 0 1 -1 ]\n";
  os << Indent << "} # IndexedLineSet\n";
  os << "                } #Shape\n"
        "            ] #children\n"
        "        } #Group\n";

  return os;
}

} //namespace CGAL

#endif // CGAL_IO_VRML_2_SEGMENT_3
#endif // CGAL_SEGMENT_3_H

#ifdef CGAL_SPHERE_3_H
#ifndef CGAL_IO_VRML_2_SPHERE_3
#define CGAL_IO_VRML_2_SPHERE_3

namespace CGAL {

template <class R >
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           const Sphere_3<R > &s)
{
  os <<      "        Group {\n"
             "            children [\n"
             "              Transform {\n"
             "                translation ";
  os <<      CGAL::to_double(s.center().x()) << " " 
     <<      CGAL::to_double(s.center().y()) << " "
     <<      CGAL::to_double(s.center().z()) << "\n";
  os <<      "                children Shape {\n"
             "                    appearance USE A1\n"
             "                    geometry\n"
             "                        Sphere { "
             "radius ";
  os <<      std::sqrt(CGAL::to_double(s.squared_radius())) <<" }\n";
  os <<      "                } #children Shape\n"
             "              } # Transform\n"
             "            ] #children\n"
             "        } #Group\n";

  return os;
}

} //namespace CGAL

#endif // CGAL_IO_VRML_2_SEGMENT_3
#endif // CGAL_SPHERE_3_H
