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
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : File_scanner_OFF.h
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// File scanner for an object in an object file format (OFF) file
// ============================================================================

#ifndef CGAL_IO_FILE_SCANNER_OFF_H
#define CGAL_IO_FILE_SCANNER_OFF_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_KNOWN_BIT_SIZE_INTEGERS_H
#include <CGAL/known_bit_size_integers.h>
#endif
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H
#ifndef CGAL_IO_FILE_HEADER_OFF_H
#include <CGAL/IO/File_header_OFF.h>
#endif // CGAL_IO_FILE_HEADER_OFF_H
#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>

CGAL_BEGIN_NAMESPACE

class File_scanner_OFF : public File_header_OFF {
    std::istream&  m_in;
    bool           normals_read;
    void skip_comment() { m_in >> skip_comment_OFF; }
public:
    File_scanner_OFF( std::istream& in, bool verbose = false)
      : File_header_OFF(verbose), m_in(in), normals_read(false) {
        in >> static_cast<File_header_OFF&>( *this);
    }
    File_scanner_OFF( std::istream& in, const File_header_OFF& header)
      : File_header_OFF(header), m_in(in), normals_read(false) {}

    std::istream& in() { return m_in; }

    // The scan_vertex() routine is provided for multiple
    // coordinate types to support parameterized polytopes.
    void scan_vertex( float&  x, float&  y, float&  z) {
        if ( binary()) {
            I_Binary_read_big_endian_float32( m_in, x);
            I_Binary_read_big_endian_float32( m_in, y);
            I_Binary_read_big_endian_float32( m_in, z);
            if ( is_homogeneous()) {
                float w;
                I_Binary_read_big_endian_float32( m_in, w);
                x /= w;
                y /= w;
                z /= w;
            }
        } else {
            skip_comment();
            m_in >> x >> y >> z;
            if ( is_homogeneous()) {
                float w;
                m_in >> w;
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( double& x, double& y, double& z) {
        if ( binary()) {
            float f;
            I_Binary_read_big_endian_float32( m_in, f);
            x = f;
            I_Binary_read_big_endian_float32( m_in, f);
            y = f;
            I_Binary_read_big_endian_float32( m_in, f);
            z = f;
            if ( is_homogeneous()) {
                I_Binary_read_big_endian_float32( m_in, f);
                x /= f;
                y /= f;
                z /= f;
            }
        } else {
            skip_comment();
            m_in >> x >> y >> z;
            if ( is_homogeneous()) {
                double w;
                m_in >> w;
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( int& x, int& y, int& z) {
        if ( binary()) {
            float fx, fy, fz;
            I_Binary_read_big_endian_float32( m_in, fx);
            I_Binary_read_big_endian_float32( m_in, fy);
            I_Binary_read_big_endian_float32( m_in, fz);
            if ( is_homogeneous()) {
                float fw;
                I_Binary_read_big_endian_float32( m_in, fw);
                x = int( fx / fw);
                y = int( fy / fw);
                y = int( fz / fw);
            } else {
                x = int(fx);
                y = int(fy);
                z = int(fz);
            }
        } else {
            skip_comment();
            if ( is_homogeneous()) {
                double fx, fy, fz, fw;
                m_in >> fx >> fy >> fz >> fw;
                x = int( fx / fw);
                y = int( fy / fw);
                y = int( fz / fw);
            } else {
                double d;
                m_in >> d;
                x = int(d);
                m_in >> d;
                y = int(d);
                m_in >> d;
                z = int(d);
            }
        }
    }

    void scan_vertex( float&  x, float&  y, float&  z, float&  w) {
        w = 1;
        if ( binary()) {
            I_Binary_read_big_endian_float32( m_in, x);
            I_Binary_read_big_endian_float32( m_in, y);
            I_Binary_read_big_endian_float32( m_in, z);
            if ( is_homogeneous())
                I_Binary_read_big_endian_float32( m_in, w);
        } else {
            skip_comment();
            m_in >> x >> y >> z;
            if ( is_homogeneous())
                m_in >> w;
        }
    }
    void scan_vertex( double& x, double& y, double& z, double& w) {
        w = 1;
        if ( binary()) {
            float f;
            I_Binary_read_big_endian_float32( m_in, f);
            x = f;
            I_Binary_read_big_endian_float32( m_in, f);
            y = f;
            I_Binary_read_big_endian_float32( m_in, f);
            z = f;
            if ( is_homogeneous()) {
                I_Binary_read_big_endian_float32( m_in, f);
                w = f;
            }
        } else {
            skip_comment();
            m_in >> x >> y >> z;
            if ( is_homogeneous())
                m_in >> w;
        }
    }
    void scan_vertex( int& x, int& y, int& z, int& w) {
        w = 1;
        if ( binary()) {
            float f;
            I_Binary_read_big_endian_float32( m_in, f);
            x = int(f);
            I_Binary_read_big_endian_float32( m_in, f);
            y = int(f);
            I_Binary_read_big_endian_float32( m_in, f);
            z = int(f);
            if ( is_homogeneous()) {
                I_Binary_read_big_endian_float32( m_in, f);
                w = int(f);
            }
        } else {
            skip_comment();
            double d;
            m_in >> d;
            x = int(d);
            m_in >> d;
            y = int(d);
            m_in >> d;
            z = int(d);
            if ( is_homogeneous()) {
                m_in >> d;
                w = int(d);
            }
        }
    }

    void scan_normal( float&  x, float&  y, float&  z) {
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                I_Binary_read_big_endian_float32( m_in, x);
                I_Binary_read_big_endian_float32( m_in, y);
                I_Binary_read_big_endian_float32( m_in, z);
                if ( is_homogeneous()) {
                    float w;
                    I_Binary_read_big_endian_float32( m_in, w);
                    x /= w;
                    y /= w;
                    z /= w;
                }
            } else {
                m_in >> x >> y >> z;
                if ( is_homogeneous()) {
                    float w;
                    m_in >> w;
                    x /= w;
                    y /= w;
                    z /= w;
                }
            }
        }
    }
    void scan_normal( double& x, double& y, double& z) {
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float fx, fy, fz;
                I_Binary_read_big_endian_float32( m_in, fx);
                I_Binary_read_big_endian_float32( m_in, fy);
                I_Binary_read_big_endian_float32( m_in, fz);
                if ( is_homogeneous()) {
                    float fw;
                    I_Binary_read_big_endian_float32( m_in, fw);
                    x = fx / fw;
                    y = fy / fw;
                    y = fz / fw;
                } else {
                    x = fx;
                    y = fy;
                    z = fz;
                }
            } else {
                if ( is_homogeneous()) {
                    float fx, fy, fz, fw;
                    m_in >> fx >> fy >> fz >> fw;
                    x = fx / fw;
                    y = fy / fw;
                    y = fz / fw;
                } else
                    m_in >> x >> y >> z;
            }
        }
    }
    void scan_normal( int& x, int& y, int& z) {
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float fx, fy, fz;
                I_Binary_read_big_endian_float32( m_in, fx);
                I_Binary_read_big_endian_float32( m_in, fy);
                I_Binary_read_big_endian_float32( m_in, fz);
                if ( is_homogeneous()) {
                    float fw;
                    I_Binary_read_big_endian_float32( m_in, fw);
                    x = int( fx / fw);
                    y = int( fy / fw);
                    y = int( fz / fw);
                } else {
                    x = int(fx);
                    y = int(fy);
                    z = int(fz);
                }
            } else {
                if ( is_homogeneous()) {
                    float fx, fy, fz, fw;
                    m_in >> fx >> fy >> fz >> fw;
                    x = int( fx / fw);
                    y = int( fy / fw);
                    y = int( fz / fw);
                } else {
                    double d;
                    m_in >> d;
                    x = int(d);
                    m_in >> d;
                    y = int(d);
                    m_in >> d;
                    z = int(d);
                }
            }
        }
    }

    void scan_normal( float&  x, float&  y, float&  z, float&  w) {
        w = 1;
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                I_Binary_read_big_endian_float32( m_in, x);
                I_Binary_read_big_endian_float32( m_in, y);
                I_Binary_read_big_endian_float32( m_in, z);
                if ( is_homogeneous())
                    I_Binary_read_big_endian_float32( m_in, w);
            } else {
                m_in >> x >> y >> z;
                if ( is_homogeneous())
                    m_in >> w;
            }
        }
    }
    void scan_normal( double& x, double& y, double& z, double& w) {
        w = 1;
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float f;
                I_Binary_read_big_endian_float32( m_in, f);
                x = f;
                I_Binary_read_big_endian_float32( m_in, f);
                y = f;
                I_Binary_read_big_endian_float32( m_in, f);
                z = f;
                if ( is_homogeneous()) {
                    I_Binary_read_big_endian_float32( m_in, f);
                    w = f;
                }
            } else {
                m_in >> x >> y >> z;
                if ( is_homogeneous())
                    m_in >> w;
            }
        }
    }
    void scan_normal( int& x, int& y, int& z, int& w) {
        w = 1;
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float f;
                I_Binary_read_big_endian_float32( m_in, f);
                x = int(f);
                I_Binary_read_big_endian_float32( m_in, f);
                y = int(f);
                I_Binary_read_big_endian_float32( m_in, f);
                z = int(f);
                if ( is_homogeneous()) {
                    I_Binary_read_big_endian_float32( m_in, f);
                    w = int(f);
                }
            } else {
                double d;
                m_in >> d;
                x = int(d);
                m_in >> d;
                y = int(d);
                m_in >> d;
                z = int(d);
                if ( is_homogeneous()) {
                    m_in >> d;
                    w = int(d);
                }
            }
        }
    }

    void skip_to_next_vertex( int current_vertex);

    void scan_facet( Integer32& size, int current_facet) {
        CGAL_assertion( current_facet < size_of_facets());
        if ( binary())
            I_Binary_read_big_endian_integer32( m_in, size);
        else {
            skip_comment();
            m_in >> size;
        }
    }

    void scan_facet_vertex_index( Integer32& index,
                                  int current_facet) {
        if ( binary())
            I_Binary_read_big_endian_integer32( m_in, index);
        else
            m_in >> index;
        if( ! m_in) {
            if ( verbose()) {
                std::cerr << " " << std::endl;
                std::cerr << "File_scanner_OFF::" << std::endl;
                std::cerr << "scan_facet_vertex_index(): input error:  "
                             "cannot read OFF file beyond facet "
                          << current_facet << "." << std::endl;
            }
            set_off_header( false);
            return;
        }
        index -= index_offset();
        if( index < 0 || index >= size_of_vertices()) {
            m_in.clear( std::ios::badbit);
            if ( verbose()) {
                std::cerr << " " << std::endl;
                std::cerr << "File_scanner_OFF::" << std::endl;
                std::cerr << "scan_facet_vertex_index(): input error: "
                             "facet " << current_facet << ": vertex index "
                          << index + index_offset() << ": is out of range."
                          << std::endl;
            }
            set_off_header( false);
            return;
        }
    }

    void skip_to_next_facet( int current_facet);
};

template < class R> inline
Point_3<R>&
file_scan_vertex( File_scanner_OFF& scanner, Point_3<R>& p) {
    typedef typename R::RT  RT;
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    if ( w == 1)
        p = Point_3<R>( RT(x), RT(y), RT(z));
    else
        p = Point_3<R>( RT(x), RT(y), RT(z), RT(w));
    return p;
}

template < class R> inline
Vector_3<R>&
file_scan_normal( File_scanner_OFF& scanner, Vector_3<R>& v) {
    typedef typename R::RT  RT;
    double x, y, z, w;
    scanner.scan_normal( x, y, z, w);
    if ( w == 1)
        v = Vector_3<R>( RT(x), RT(y), RT(z));
    else
        v = Vector_3<R>( RT(x), RT(y), RT(z), RT(w));
    return v;
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_SCANNER_OFF_H //
// EOF //
