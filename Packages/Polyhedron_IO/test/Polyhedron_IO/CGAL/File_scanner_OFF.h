#line 12 "cgal_header.fw"
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
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// File scanner for an object in an object file format (OFF) file
// ============================================================================
#line 120 "polyhedron_io.fw"

#ifndef CGAL_FILE_SCANNER_OFF_H
#define CGAL_FILE_SCANNER_OFF_H 1
#line 915 "polyhedron_io.fw"
#ifndef _STDDEF_H
#include <stddef.h>
#endif

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_BINARY_FILE_IO_H
#include <CGAL/binary_file_io.h>
#endif

#ifndef CGAL_FILE_HEADER_OFF_H
#include <CGAL/File_header_OFF.h>
#endif

// Forward declarations.
class istream;

class CGAL_File_scanner_OFF : public CGAL_File_header_OFF {
    istream&  _in;
    bool      normals_read;
public:
    CGAL_File_scanner_OFF( istream& in)
        : _in( in), normals_read( false), CGAL_File_header_OFF( in) {}
    CGAL_File_scanner_OFF( istream& in,
                          const CGAL_File_header_OFF& header)
        : _in( in), normals_read( false), CGAL_File_header_OFF( header) {}
    // The scan_vertex() routine is provided for multiple
    // coordinate types to support parameterized polytopes.
    void scan_vertex( float&  x, float&  y, float&  z) {
        if ( binary) {
            CGAL__Binary_read_float32( _in, x);
            CGAL__Binary_read_float32( _in, y);
            CGAL__Binary_read_float32( _in, z);
            if ( tag4) {
                float w;
                CGAL__Binary_read_float32( _in, w);
                x /= w;
                y /= w;
                z /= w;
            }
        } else {
            // Skip any possible comment.
            char c;
            while ( (_in >> c) && c == '#') {
                while ( _in.get(c) && c != '\n')
                    ;
            }
            _in.putback(c);
            _in >> x >> y >> z;
            if ( tag4) {
                float w;
                _in >> w;
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( double& x, double& y, double& z) {
        if ( binary) {
            float f;
            CGAL__Binary_read_float32( _in, f);
            x = f;
            CGAL__Binary_read_float32( _in, f);
            y = f;
            CGAL__Binary_read_float32( _in, f);
            z = f;
            if ( tag4) {
                CGAL__Binary_read_float32( _in, f);
                x /= f;
                y /= f;
                z /= f;
            }
        } else {
            // Skip any possible comment.
            char c;
            while ( (_in >> c) && c == '#') {
                while ( _in.get(c) && c != '\n')
                    ;
            }
            _in.putback(c);
            _in >> x >> y >> z;
            if ( tag4) {
                double w;
                _in >> w;
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( short&  x, short&  y, short&  z) {
        if ( binary) {
            float fx, fy, fz;
            CGAL__Binary_read_float32( _in, fx);
            CGAL__Binary_read_float32( _in, fy);
            CGAL__Binary_read_float32( _in, fz);
            if ( tag4) {
                float fw;
                CGAL__Binary_read_float32( _in, fw);
                x = short( fx / fw);
                y = short( fy / fw);
                y = short( fz / fw);
            } else {
                x = short(fx);
                y = short(fy);
                z = short(fz);
            }
        } else {
            // Skip any possible comment.
            char c;
            while ( (_in >> c) && c == '#') {
                while ( _in.get(c) && c != '\n')
                    ;
            }
            _in.putback(c);
            if ( tag4) {
                float fx, fy, fz, fw;
                _in >> fx >> fy >> fz >> fw;
                x = short( fx / fw);
                y = short( fy / fw);
                y = short( fz / fw);
            } else
                _in >> x >> y >> z;
        }
    }
    void scan_vertex( int&    x, int&    y, int&    z) {
        if ( binary) {
            float fx, fy, fz;
            CGAL__Binary_read_float32( _in, fx);
            CGAL__Binary_read_float32( _in, fy);
            CGAL__Binary_read_float32( _in, fz);
            if ( tag4) {
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
            // Skip any possible comment.
            char c;
            while ( (_in >> c) && c == '#') {
                while ( _in.get(c) && c != '\n')
                    ;
            }
            _in.putback(c);
            if ( tag4) {
                double fx, fy, fz, fw;
                _in >> fx >> fy >> fz >> fw;
                x = int( fx / fw);
                y = int( fy / fw);
                y = int( fz / fw);
            } else
                _in >> x >> y >> z;
        }
    }
    void scan_vertex( long&   x, long&   y, long&   z) {
        if ( binary) {
            float fx, fy, fz;
            CGAL__Binary_read_float32( _in, fx);
            CGAL__Binary_read_float32( _in, fy);
            CGAL__Binary_read_float32( _in, fz);
            if ( tag4) {
                float fw;
                CGAL__Binary_read_float32( _in, fw);
                x = long( fx / fw);
                y = long( fy / fw);
                y = long( fz / fw);
            } else {
                x = long(fx);
                y = long(fy);
                z = long(fz);
            }
        } else {
            // Skip any possible comment.
            char c;
            while ( (_in >> c) && c == '#') {
                while ( _in.get(c) && c != '\n')
                    ;
            }
            _in.putback(c);
            if ( tag4) {
                double fx, fy, fz, fw;
                _in >> fx >> fy >> fz >> fw;
                x = long( fx / fw);
                y = long( fy / fw);
                y = long( fz / fw);
            } else
                _in >> x >> y >> z;
        }
    }

    void scan_normal( float&  x, float&  y, float&  z) {
        if ( normals) {
            normals_read = true;
            if ( binary) {
                CGAL__Binary_read_float32( _in, x);
                CGAL__Binary_read_float32( _in, y);
                CGAL__Binary_read_float32( _in, z);
                if ( tag4) {
                    float w;
                    CGAL__Binary_read_float32( _in, w);
                    x /= w;
                    y /= w;
                    z /= w;
                }
            } else {
                _in >> x >> y >> z;
                if ( tag4) {
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
        if ( normals) {
            normals_read = true;
            if ( binary) {
                float fx, fy, fz;
                CGAL__Binary_read_float32( _in, fx);
                CGAL__Binary_read_float32( _in, fy);
                CGAL__Binary_read_float32( _in, fz);
                if ( tag4) {
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
                if ( tag4) {
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
    void scan_normal( short&  x, short&  y, short&  z) {
        if ( normals) {
            normals_read = true;
            if ( binary) {
                float fx, fy, fz;
                CGAL__Binary_read_float32( _in, fx);
                CGAL__Binary_read_float32( _in, fy);
                CGAL__Binary_read_float32( _in, fz);
                if ( tag4) {
                    float fw;
                    CGAL__Binary_read_float32( _in, fw);
                    x = short( fx / fw);
                    y = short( fy / fw);
                    y = short( fz / fw);
                } else {
                    x = short(fx);
                    y = short(fy);
                    z = short(fz);
                }
            } else {
                if ( tag4) {
                    float fx, fy, fz, fw;
                    _in >> fx >> fy >> fz >> fw;
                    x = short( fx / fw);
                    y = short( fy / fw);
                    y = short( fz / fw);
                } else
                    _in >> x >> y >> z;
            }
        }
    }
    void scan_normal( int&    x, int&    y, int&    z) {
        if ( normals) {
            normals_read = true;
            if ( binary) {
                float fx, fy, fz;
                CGAL__Binary_read_float32( _in, fx);
                CGAL__Binary_read_float32( _in, fy);
                CGAL__Binary_read_float32( _in, fz);
                if ( tag4) {
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
                if ( tag4) {
                    float fx, fy, fz, fw;
                    _in >> fx >> fy >> fz >> fw;
                    x = int( fx / fw);
                    y = int( fy / fw);
                    y = int( fz / fw);
                } else
                    _in >> x >> y >> z;
            }
        }
    }
    void scan_normal( long&   x, long&   y, long&   z) {
        if ( normals) {
            normals_read = true;
            if ( binary) {
                float fx, fy, fz;
                CGAL__Binary_read_float32( _in, fx);
                CGAL__Binary_read_float32( _in, fy);
                CGAL__Binary_read_float32( _in, fz);
                if ( tag4) {
                    float fw;
                    CGAL__Binary_read_float32( _in, fw);
                    x = long( fx / fw);
                    y = long( fy / fw);
                    y = long( fz / fw);
                } else {
                    x = long(fx);
                    y = long(fy);
                    z = long(fz);
                }
            } else {
                if ( tag4) {
                    float fx, fy, fz, fw;
                    _in >> fx >> fy >> fz >> fw;
                    x = long( fx / fw);
                    y = long( fy / fw);
                    y = long( fz / fw);
                } else
                    _in >> x >> y >> z;
            }
        }
    }

    void skip_to_next_vertex( int current_vertex);

    void scan_facet( CGAL_Integer32& size, int current_facet) {
        CGAL_assertion( current_facet < n_facets);
        if ( binary)
            CGAL__Binary_read_integer32( _in, size);
        else {
            // Skip any possible comment.
            char c;
            while ( (_in >> c) && c == '#') {
                while ( _in.get(c) && c != '\n')
                    ;
            }
            _in.putback(c);
            _in >> size;
        }
        if( ! skel && size < 3) {
            cerr << " " << endl;
            cerr << "CGAL_Polyhedron_scan_OFF<Traits>::" << endl;
            cerr << "build(): input error: facet " << current_facet
                 << " has less than 3 vertices." << endl;
            abort();
        }
    }

    void scan_facet_vertex_index( CGAL_Integer32& index, int current_facet) {
        if ( binary)
            CGAL__Binary_read_integer32( _in, index);
        else
            _in >> index;
        if( ! _in) {
            cerr << " " << endl;
            cerr << "CGAL_Polyhedron_scan_OFF<Traits>::" << endl;
            cerr << "build(): input error: cannot read OFF file "
                    "beyond facet " << current_facet << "." << endl;
            abort();
        }
        index -= offset;
        if( index < 0 || index >= n_vertices) {
            cerr << " " << endl;
            cerr << "CGAL_Polyhedron_scan_OFF<Traits>::" << endl;
            cerr << "build(): input error: facet " << current_facet
                 << ": vertex index " << index + offset
                 << ": is out of range." << endl;
            abort();
        }
    }

    void skip_to_next_facet( int current_facet);
};
#line 123 "polyhedron_io.fw"
  
#endif // CGAL_FILE_SCANNER_OFF_H //
// EOF //
