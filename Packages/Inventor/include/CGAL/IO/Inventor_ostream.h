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
// release       :
// date          :
//
// file          : include/CGAL/IO/Inventor_ostream.h
// source        : web/Inventor_ostream.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Lutz Kettner <kettner@inf.ethz.ch>
//                 Herve bronnimann
//                 <Herve.Bronnimann@sophia.inria.fr>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_IO_INVENTOR_OSTREAM_H
#define CGAL_IO_INVENTOR_OSTREAM_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_IOSTREAM_H
#include <iostream.h>
#define CGAL_PROTECT_IOSTREAM_H
#endif // CGAL_PROTECT_IOSTREAM_H

// OpenInventor and VRML 1.0 are quite similar formats, so
// output operators could be shared if they use the following
// base class, which is common for both output streams.

class CGAL_Inventor_ostream_base {
private:
    ostream*  m_os;
public:
    CGAL_Inventor_ostream_base()           : m_os(0)  {}
    CGAL_Inventor_ostream_base(ostream& o) : m_os(&o) {}
    ~CGAL_Inventor_ostream_base()  { close(); }
    void open(ostream& o)        { m_os = &o; }
    void close() {
        if ( m_os)
            os() << endl;
        m_os = 0;
    }
    typedef const void* Const_void_ptr;
    operator Const_void_ptr () const {
        if ( m_os)
            return *m_os;
        return 0;
    }
    ostream& os() {
        // The behaviour if m_os == 0 could be changed to return
        // cerr or a file handle to /dev/null. The latter one would
        // mimick the behaviour that one can still use a stream with
        // an invalid stream, but without producing any output.
        CGAL_assertion( m_os);
        return *m_os;
    }
};


class CGAL_Inventor_ostream : public  CGAL_Inventor_ostream_base
{
public:
    CGAL_Inventor_ostream() {}
    CGAL_Inventor_ostream(ostream& o) : CGAL_Inventor_ostream_base(o) {
        header();
    }
    void open(ostream& o) {
        CGAL_Inventor_ostream_base::open(o);
        header();
    }
private:
    void header() {
        os() << "#Inventor V2.0 ascii" << endl;
        os() << "# File written with the help of the CGAL Library" << endl;
    }
};


#endif // CGAL_IO_INVENTOR_OSTREAM_H


#ifdef CGAL_TETRAHEDRON_3_H
#ifndef CGAL_INVENTOR_TETRAHEDRON_3
#define CGAL_INVENTOR_TETRAHEDRON_3


template <class R >
CGAL_Inventor_ostream&
operator<<(CGAL_Inventor_ostream& os,
           const CGAL_Tetrahedron_3<R > &t)
{
  const char *Indent = "   ";
  os.os() << "\n Separator {";
  os.os() << "\n   Coordinate3 { \n"
          << Indent << "point [\n"
          << Indent << "  "
          << CGAL_to_double(t[0].x()) << " "
          << CGAL_to_double(t[0].y()) << " "
          << CGAL_to_double(t[0].z()) << " ,\n"
          << Indent << "  "
          << CGAL_to_double(t[1].x()) << " "
          << CGAL_to_double(t[1].y()) << " "
          << CGAL_to_double(t[1].z()) << " ,\n"
          << Indent << "  "
          << CGAL_to_double(t[2].x()) << " "
          << CGAL_to_double(t[2].y()) << " "
          << CGAL_to_double(t[2].z()) << " ,\n"
          << Indent << "  "
          << CGAL_to_double(t[3].x()) << " "
          << CGAL_to_double(t[3].y()) << " "
          << CGAL_to_double(t[3].z()) << " ]"
          << "\n   } #Coordinate3" ;
  os.os() << "\n   IndexedFaceSet {"
          << Indent << "coordIndex  [ 0,1,2,-1, 1,3,2,-1,\n"
          << Indent << "              0,2,3,-1, 0,3,1,-1 ]\n"
          << "\n   } #IndexedFaceSet"
          << "\n } #Separator\n";
  return os;
}


#endif // CGAL_INVENTOR_TETRAHEDRON_3
#endif // CGAL_TETRAHEDRON_3_H
