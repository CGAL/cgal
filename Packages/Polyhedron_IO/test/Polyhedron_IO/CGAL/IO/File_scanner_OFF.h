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
// file          : CGAL_File_scanner_OFF.h
// package       : $CGAL_Package: $
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

#ifndef CGAL_PROTECT_STDDEF_H
#include <stddef.h>
#define CGAL_PROTECT_STDDEF_H
#endif // CGAL_PROTECT_STDDEF_H

#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H

#ifndef CGAL_IO_FILE_HEADER_OFF_H
#include <CGAL/IO/File_header_OFF.h>
#endif // CGAL_IO_FILE_HEADER_OFF_H

// Forward declarations.
class istream;

class CGAL_File_scanner_OFF : public CGAL_File_header_OFF {
    istream&  _in;
    bool      normals_read;
    void skip_comment() {
        char c;
        while ( (_in >> c) && c == '#') {
            while ( _in.get(c) && c != '\n')
                ;
        }
        _in.putback(c);
    }
public:
    CGAL_File_scanner_OFF( istream& in, bool verbose = false)
        : _in( in), normals_read( false),
          CGAL_File_header_OFF( in, verbose) {}
    CGAL_File_scanner_OFF( istream& in,
                          const CGAL_File_header_OFF& header)
        : _in( in), normals_read( false), CGAL_File_header_OFF( header) {}
    // The scan_vertex() routine is provided for multiple
    // coordinate types to support parameterized polytopes.
    void scan_vertex( float&  x, float&  y, float&  z) {
        if ( _binary) {
            CGAL__Binary_read_float32( _in, x);
            CGAL__Binary_read_float32( _in, y);
            CGAL__Binary_read_float32( _in, z);
            if ( _tag4) {
                float w;
                CGAL__Binary_read_float32( _in, w);
                x /= w;
                y /= w;
                z /= w;
            }
        } else {
            skip_comment();
            _in >> x >> y >> z;
            if ( _tag4) {
                float w;
                _in >> w;
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( double& x, double& y, double& z) {
        if ( _binary) {
            float f;
            CGAL__Binary_read_float32( _in, f);
            x = f;
            CGAL__Binary_read_float32( _in, f);
            y = f;
            CGAL__Binary_read_float32( _in, f);
            z = f;
            if ( _tag4) {
                CGAL__Binary_read_float32( _in, f);
                x /= f;
                y /= f;
                z /= f;
            }
        } else {
            skip_comment();
            _in >> x >> y >> z;
            if ( _tag4) {
                double w;
                _in >> w;
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( int&    x, int&    y, int&    z) {
        if ( _binary) {
            float fx, fy, fz;
            CGAL__Binary_read_float32( _in, fx);
            CGAL__Binary_read_float32( _in, fy);
            CGAL__Binary_read_float32( _in, fz);
            if ( _tag4) {
                float fw;
                CGAL__Binary_read_float32( _in, fw);
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
            if ( _tag4) {
                double fx, fy, fz, fw;
                _in >> fx >> fy >> fz >> fw;
                x = int( fx / fw);
                y = int( fy / fw);
                y = int( fz / fw);
            } else {
                double d;
                _in >> d;
                x = int(d);
                _in >> d;
                y = int(d);
                _in >> d;
                z = int(d);
            }
        }
    }

    void scan_vertex( float&  x, float&  y, float&  z, float&  w) {
        w = 1;
        if ( _binary) {
            CGAL__Binary_read_float32( _in, x);
            CGAL__Binary_read_float32( _in, y);
            CGAL__Binary_read_float32( _in, z);
            if ( _tag4)
                CGAL__Binary_read_float32( _in, w);
        } else {
            skip_comment();
            _in >> x >> y >> z;
            if ( _tag4)
                _in >> w;
        }
    }
    void scan_vertex( double& x, double& y, double& z, double& w) {
        w = 1;
        if ( _binary) {
            float f;
            CGAL__Binary_read_float32( _in, f);
            x = f;
            CGAL__Binary_read_float32( _in, f);
            y = f;
            CGAL__Binary_read_float32( _in, f);
            z = f;
            if ( _tag4) {
                CGAL__Binary_read_float32( _in, f);
                w = f;
            }
        } else {
            skip_comment();
            _in >> x >> y >> z;
            if ( _tag4)
                _in >> w;
        }
    }
    void scan_vertex( int&    x, int&    y, int&    z, int&    w) {
        w = 1;
        if ( _binary) {
            float f;
            CGAL__Binary_read_float32( _in, f);
            x = int(f);
            CGAL__Binary_read_float32( _in, f);
            y = int(f);
            CGAL__Binary_read_float32( _in, f);
            z = int(f);
            if ( _tag4) {
                CGAL__Binary_read_float32( _in, f);
                w = int(f);
            }
        } else {
            skip_comment();
            double d;
            _in >> d;
            x = int(d);
            _in >> d;
            y = int(d);
            _in >> d;
            z = int(d);
            if ( _tag4) {
                _in >> d;
                w = int(d);
            }
        }
    }

    void scan_normal( float&  x, float&  y, float&  z) {
        if ( _normals) {
            normals_read = true;
            if ( _binary) {
                CGAL__Binary_read_float32( _in, x);
                CGAL__Binary_read_float32( _in, y);
                CGAL__Binary_read_float32( _in, z);
                if ( _tag4) {
                    float w;
                    CGAL__Binary_read_float32( _in, w);
                    x /= w;
                    y /= w;
                    z /= w;
                }
            } else {
                _in >> x >> y >> z;
                if ( _tag4) {
                    float w;
                    _in >> w;
                    x /= w;
                    y /= w;
                    z /= w;
                }
            }
        }
    }
    void scan_normal( double& x, double& y, double& z) {
        if ( _normals) {
            normals_read = true;
            if ( _binary) {
                float fx, fy, fz;
                CGAL__Binary_read_float32( _in, fx);
                CGAL__Binary_read_float32( _in, fy);
                CGAL__Binary_read_float32( _in, fz);
                if ( _tag4) {
                    float fw;
                    CGAL__Binary_read_float32( _in, fw);
                    x = fx / fw;
                    y = fy / fw;
                    y = fz / fw;
                } else {
                    x = fx;
                    y = fy;
                    z = fz;
                }
            } else {
                if ( _tag4) {
                    float fx, fy, fz, fw;
                    _in >> fx >> fy >> fz >> fw;
                    x = fx / fw;
                    y = fy / fw;
                    y = fz / fw;
                } else
                    _in >> x >> y >> z;
            }
        }
    }
    void scan_normal( int&    x, int&    y, int&    z) {
        if ( _normals) {
            normals_read = true;
            if ( _binary) {
                float fx, fy, fz;
                CGAL__Binary_read_float32( _in, fx);
                CGAL__Binary_read_float32( _in, fy);
                CGAL__Binary_read_float32( _in, fz);
                if ( _tag4) {
                    float fw;
                    CGAL__Binary_read_float32( _in, fw);
                    x = int( fx / fw);
                    y = int( fy / fw);
                    y = int( fz / fw);
                } else {
                    x = int(fx);
                    y = int(fy);
                    z = int(fz);
                }
            } else {
                if ( _tag4) {
                    float fx, fy, fz, fw;
                    _in >> fx >> fy >> fz >> fw;
                    x = int( fx / fw);
                    y = int( fy / fw);
                    y = int( fz / fw);
                } else {
                    double d;
                    _in >> d;
                    x = int(d);
                    _in >> d;
                    y = int(d);
                    _in >> d;
                    z = int(d);
                }
            }
        }
    }

    void scan_normal( float&  x, float&  y, float&  z, float&  w) {
        w = 1;
        if ( _normals) {
            normals_read = true;
            if ( _binary) {
                CGAL__Binary_read_float32( _in, x);
                CGAL__Binary_read_float32( _in, y);
                CGAL__Binary_read_float32( _in, z);
                if ( _tag4)
                    CGAL__Binary_read_float32( _in, w);
            } else {
                _in >> x >> y >> z;
                if ( _tag4)
                    _in >> w;
            }
        }
    }
    void scan_normal( double& x, double& y, double& z, double& w) {
        w = 1;
        if ( _normals) {
            normals_read = true;
            if ( _binary) {
                float f;
                CGAL__Binary_read_float32( _in, f);
                x = f;
                CGAL__Binary_read_float32( _in, f);
                y = f;
                CGAL__Binary_read_float32( _in, f);
                z = f;
                if ( _tag4) {
                    CGAL__Binary_read_float32( _in, f);
                    w = f;
                }
            } else {
                _in >> x >> y >> z;
                if ( _tag4)
                    _in >> w;
            }
        }
    }
    void scan_normal( int&    x, int&    y, int&    z, int&    w) {
        w = 1;
        if ( _normals) {
            normals_read = true;
            if ( _binary) {
                float f;
                CGAL__Binary_read_float32( _in, f);
                x = int(f);
                CGAL__Binary_read_float32( _in, f);
                y = int(f);
                CGAL__Binary_read_float32( _in, f);
                z = int(f);
                if ( _tag4) {
                    CGAL__Binary_read_float32( _in, f);
                    w = int(f);
                }
            } else {
                double d;
                _in >> d;
                x = int(d);
                _in >> d;
                y = int(d);
                _in >> d;
                z = int(d);
                if ( _tag4) {
                    _in >> d;
                    w = int(d);
                }
            }
        }
    }

    void skip_to_next_vertex( int current_vertex);

    void scan_facet( CGAL_Integer32& size, int current_facet) {
        CGAL_assertion( current_facet < n_facets);
        if ( _binary)
            CGAL__Binary_read_integer32( _in, size);
        else {
            skip_comment();
            _in >> size;
        }
    }

    void scan_facet_vertex_index( CGAL_Integer32& index,
                                  int current_facet) {
        if ( _binary)
            CGAL__Binary_read_integer32( _in, index);
        else
            _in >> index;
        if( ! _in) {
            if ( _verbose) {
                cerr << " " << endl;
                cerr << "CGAL_File_scanner_OFF::" << endl;
                cerr << "scan_facet_vertex_index(): input error: cannot "
                        "read OFF file beyond facet " << current_facet
                     << "." << endl;
            }
            return;
        }
        index -= _offset;
        if( index < 0 || index >= n_vertices) {
            _in.clear( ios::badbit);
            if ( _verbose) {
                cerr << " " << endl;
                cerr << "CGAL_File_scanner_OFF::" << endl;
                cerr << "scan_facet_vertex_index(): input error: facet "
                     << current_facet << ": vertex index " << index+_offset
                     << ": is out of range." << endl;
            }
            return;
        }
    }

    void skip_to_next_facet( int current_facet);
};

#ifdef CGAL_IO_FILE_SCANNER_OFF_H

template < class R> inline
CGAL_Point_3<R>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<R>& p) {
    typedef typename R::RT  RT;
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    if ( w == 1)
        p = CGAL_Point_3<R>( RT(x), RT(y), RT(z));
    else
        p = CGAL_Point_3<R>( RT(x), RT(y), RT(z), RT(w));
    return p;
}

template < class R> inline
CGAL_Vector_3<R>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<R>& v) {
    typedef typename R::RT  RT;
    double x, y, z, w;
    scanner.scan_normal( x, y, z, w);
    if ( w == 1)
        v = CGAL_Vector_3<R>( RT(x), RT(y), RT(z));
    else
        v = CGAL_Vector_3<R>( RT(x), RT(y), RT(z), RT(w));
    return v;
}

#else  // CGAL_IO_FILE_SCANNER_OFF_H //
// This is code from the CEBAP project now useless in CGAL.

#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif
#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif
#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif
#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif

template < class Pt, class RT> inline
Pt&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner, Pt& p,
                      CGAL_Cartesian_tag, RT*) {
    double x, y, z;
    scanner.scan_vertex( x, y, z);
    p = Pt( RT(x), RT(y), RT(z));
    return p;
}

inline
CGAL_Point_3<float>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<float>& p,
                      CGAL_Cartesian_tag,
                      float*) {
    scanner.scan_vertex( p.x(), p.y(), p.z());
    return p;
}

inline
CGAL_Point_3<double>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<double>& p,
                      CGAL_Cartesian_tag,
                      double*) {
    scanner.scan_vertex( p.x(), p.y(), p.z());
    return p;
}

inline
CGAL_Point_3<int>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<int>& p,
                      CGAL_Cartesian_tag,
                      int*) {
    scanner.scan_vertex( p.x(), p.y(), p.z());
    return p;
}

template < class Pt, class RT> inline
Pt&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner, Pt& p,
                      CGAL_Rational_tag, RT*) {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    p = Pt( RT(x), RT(y), RT(z), RT(w));
    return p;
}

inline
CGAL_Point_3<float,float>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<float,float>& p,
                      CGAL_Rational_tag, float*) {
    scanner.scan_vertex( p.hx(), p.hy(), p.hz(), p.rw());
    return p;
}

inline
CGAL_Point_3<double,double>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<double,double>& p,
                      CGAL_Rational_tag, double*) {
    scanner.scan_vertex( p.hx(), p.hy(), p.hz(), p.rw());
    return p;
}

inline
CGAL_Point_3<double,int>&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner,
                      CGAL_Point_3<double,int>& p,
                      CGAL_Rational_tag, int*) {
    scanner.scan_vertex( p.hx(), p.hy(), p.hz(), p.rw());
    return p;
}

template < class Pt> inline
Pt&
CGAL_file_scan_vertex( CGAL_File_scanner_OFF& scanner, Pt& p) {
    typedef typename Pt::RT  RT;
    return CGAL_file_scan_vertex( scanner, p, CGAL_query_representation(p),
                                 (RT*)(0));
}





template < class V, class RT> inline
V&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner, V& v,
                      CGAL_Cartesian_tag, RT*) {
    double x, y, z;
    scanner.scan_normal( x, y, z);
    v = V( RT(x), RT(y), RT(z));
    return v;
}

inline
CGAL_Vector_3<float>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<float>& v,
                      CGAL_Cartesian_tag,
                      float*) {
    scanner.scan_normal( v.x(), v.y(), v.z());
    return v;
}

inline
CGAL_Vector_3<double>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<double>& v,
                      CGAL_Cartesian_tag, double*) {
    scanner.scan_normal( v.x(), v.y(), v.z());
    return v;
}

inline
CGAL_Vector_3<int>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<int>& v,
                      CGAL_Cartesian_tag,
                      int*) {
    scanner.scan_normal( v.x(), v.y(), v.z());
    return v;
}

template < class V, class RT> inline
V&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner, V& v,
                      CGAL_Rational_tag, RT*) {
    double x, y, z, w;
    scanner.scan_normal( x, y, z, w);
    v = V( RT(x), RT(y), RT(z), RT(w));
    return v;
}

inline
CGAL_Vector_3<float,float>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<float,float>& v,
                      CGAL_Rational_tag, float*) {
    scanner.scan_normal( v.hx(), v.hy(), v.hz(), v.rw());
    return v;
}

inline
CGAL_Vector_3<double,double>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<double,double>& v,
                      CGAL_Rational_tag, double*) {
    scanner.scan_normal( v.hx(), v.hy(), v.hz(), v.rw());
    return v;
}

inline
CGAL_Vector_3<double,int>&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner,
                      CGAL_Vector_3<double,int>& v,
                      CGAL_Rational_tag, int*) {
    scanner.scan_normal( v.hx(), v.hy(), v.hz(), v.rw());
    return v;
}

template < class V> inline
V&
CGAL_file_scan_normal( CGAL_File_scanner_OFF& scanner, V& v) {
    typedef typename V::RT  RT;
    return CGAL_file_scan_normal( scanner, v,
                                 CGAL_query_representation(v),
                                 (RT*)(0));
}

#endif // CGAL_IO_FILE_SCANNER_OFF_H //
#endif // CGAL_IO_FILE_SCANNER_OFF_H //
// EOF //
