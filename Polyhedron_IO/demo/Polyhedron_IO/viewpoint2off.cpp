// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Lutz Kettner

// Copies from Viewpoint ASCII format into OFF format.
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/cstdint.hpp>

using namespace std;

typedef  CGAL::Cartesian<float>                 Kernel;
typedef  Kernel::Point_3                        Point;
typedef  vector<Point>                          Point_vector;
typedef  vector<std::size_t>                            Facet;
typedef  vector<Facet>                          Facet_vector;

bool  verbose      = false;
bool  binary       = false;
bool  skel         = false;
bool  noc          = false;
bool  binary_mesh  = false;
bool  ascii_mesh   = false;
bool  normals_file = false;
bool  no_normals   = false;

const char* normals_filename = 0;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
#if CGAL_CFG_INCOMPLETE_TYPE_BUG_5
    Point _work_around_incomplete_type;
#endif  // CGAL_CFG_INCOMPLETE_TYPE_BUG_5
    int n = 0;
    char *filename[3];
    bool help = false;
    for ( int i = 1; i < argc; i++) {
        // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
        else if ( strcmp( "-b", argv[i]) == 0)
            binary = true;
        else if ( strcmp( "-skel", argv[i]) == 0)
            skel = true;
        else if ( strcmp( "-noc", argv[i]) == 0)
            noc = true;
        else if ( strcmp( "-mesh", argv[i]) == 0)
            binary_mesh = true;
        else if ( strcmp( "-ascii", argv[i]) == 0)
            ascii_mesh = true;
        else if ( strcmp( "-no_normals", argv[i]) == 0)
            no_normals = true;
        else if ( strcmp( "-normals", argv[i]) == 0) {
            normals_file = true;
            i++;
            if ( i < argc) {
                normals_filename = argv[i];
            } else {
                cerr << argv[0] << ": error: -normals need a filename as "
                     "parameter." << endl;
                help = true;
            }
        } else if ( (strcmp( "-h", argv[i]) == 0) ||
                  (strcmp( "-help", argv[i]) == 0))
            help = true;
        else if ( n < 3 ) {
            filename[n ++] = argv[i];
        } else {
	    ++n;
            break;
	}
    }
    if ((n < 1) || (n > 3) || help) {
        if ( ! help)
            cerr << "Error: in parameter list" << endl;
        cerr << "Usage: " << argv[0]
             << " [<options>] <infile.coo> <infile.ele> [<outfile>]" << endl;
        cerr << "   or: " << argv[0]
             << " [<options>] <infile.msh> [<outfile>]" << endl;
        cerr << "       convert an object from Viewpoint formats into OFF."
             << endl;
        cerr << "       -mesh   forces input to be read from mesh file."
             << endl;
        cerr << "       -ascii  forces input to be read from ascii mesh file."
             << endl;
        cerr << "       -normals <filename>   reads normals file (in "
                "polygon format)." << endl;
        cerr << "       -no_normals           ignore normals." << endl;
        cerr << "       -b      binary output (default is ASCII)." << endl;
        cerr << "       -skel   Geomview SKEL format." << endl;
        cerr << "       -noc    no comments in file." << endl;
        cerr << "       -v      verbose." << endl;
        exit( ! help);
    }

    CGAL::Verbose_ostream vout( verbose);
    vout << argv[0] << ": verbosity on." << endl;

    ifstream     in( filename[0]);
    if ( ! in) {
        cerr << argv[0] << ": error: cannot open file '"<< filename[0]
         << "' for reading." <<endl;
        exit( 1);
    }

    Point_vector  points;
    Point_vector  normals;
    Facet_vector  facets;
    points.reserve(  64000);  // Avoid too much reallocation
    normals.reserve( 64000);
    facets.reserve(  64000);
    float x,y,z;
    char  c,d,e;
    in.get(c);
    if ( c < '0' && c > '9')
        binary_mesh = 1;
    in.putback(c);
    boost::int32_t    number;
    if ( ! binary_mesh) {
        in >> number;
        if ( number != 1)
            ascii_mesh = 1;
    }
    if ( binary_mesh) {
        vout << "Scanning Viewpoint binary mesh data file ..." << endl;
        CGAL::I_Binary_read_big_endian_integer32( in, number);
        int tmesh = 0;
        while( in) {
            tmesh ++;
            for ( int j = 0; j < number; j++) {
                if ( j > 1) {
                    facets.push_back( Facet());
                    Facet& facet = (*(facets.end() - 1));
                    facet.push_back( points.size());
                    if (j & 1) {
                        facet.push_back( points.size()-2);
                        facet.push_back( points.size()-1);
                    } else {
                        facet.push_back( points.size()-1);
                        facet.push_back( points.size()-2);
                    }
                }
                // Scan vertex coordinates.
                CGAL::I_Binary_read_big_endian_float32( in, x);
                CGAL::I_Binary_read_big_endian_float32( in, y);
                CGAL::I_Binary_read_big_endian_float32( in, z);
                points.push_back( Point( x, y, z));
            }
            for ( int k = 0; k < number; k++) {
                // Scan vertex normal.
                CGAL::I_Binary_read_big_endian_float32( in, x);
                CGAL::I_Binary_read_big_endian_float32( in, y);
                CGAL::I_Binary_read_big_endian_float32( in, z);
                if ( ! no_normals)
                    normals.push_back( Point( x, y, z));
            }
            CGAL::I_Binary_read_big_endian_integer32( in, number);
        }
        in.close();
        vout << points.size() << " vertex coordinates read." << endl;
        vout << facets.size() << " facets read." << endl;
        vout << tmesh << " triangle meshes read." << endl;
    } else if ( ascii_mesh) {
        vout << "Scanning Viewpoint ASCII mesh data file ..." << endl;
        int tmesh = 0;
        while( in) {
            tmesh ++;
            for ( int j = 0; j < number; j++) {
                if ( j > 1) {
                    facets.push_back( Facet());
                    Facet& facet = (*(facets.end() - 1));
                    facet.push_back( points.size());
                    if (j & 1) {
                        facet.push_back( points.size()-2);
                        facet.push_back( points.size()-1);
                    } else {
                        facet.push_back( points.size()-1);
                        facet.push_back( points.size()-2);
                    }
                }
                // Scan vertex coordinates.
                in >> x >> y >> z;
                points.push_back( Point( x, y, z));
            }
            for ( int k = 0; k < number; k++) {
                // Scan vertex normal.
                in >> x >> y >> z;
                if ( ! no_normals)
                    normals.push_back( Point( x, y, z));
            }
            in >> number;
        }
        in.close();
        vout << points.size() << " vertex coordinates read." << endl;
        vout << facets.size() << " facets read." << endl;
        vout << tmesh << " triangle meshes read." << endl;
    } else {
        vout << "Scanning Viewpoint polygon data files ..." << endl;
        vout << "    ... start with coordinates file ..." << endl;
        if ( n < 2) {
            cerr << argv[0] << ": error: two input filenames needed." << endl;
            exit( 1);
        }
        while( in) {
            in >> c >> x >> d >> y >> e >> z;
            if ( ! in || (size_t)number != (points.size() + 1) ||
                 c != ',' || d != ',' || e != ',') {
                cerr << argv[0] << ": error: cannot read line "
                     << points.size() + 1 << " from file '" << filename[0]
                     << "'." <<endl;
                exit( 1);
            }
            points.push_back( Point( x, y, z));
            in >> number;
        }
        in.close();
        vout << points.size() << " vertex coordinates read." << endl;

        if ( normals_file && ! no_normals) {
            vout << "    ... next is normals file ..." << endl;
            in.open( normals_filename);
            if ( ! in) {
                cerr << argv[0] << ": error: cannot open file '"
                     << normals_filename << "' for reading." <<endl;
                exit( 1);
            }
            in >> x;
            while( in) {
                in >> y >> z;
                if ( ! in ) {
                    cerr << argv[0] << ": error: cannot read line "
                         << points.size() + 1 << " from file '"
                         << normals_filename << "'." <<endl;
                    exit( 1);
                }
                normals.push_back( Point( x, y, z));
                in >> x;
            }
            in.close();
            vout << normals.size() << " normals read." << endl;
            if ( normals.size() != points.size()) {
                cerr << argv[0] << ": error: number of points and normals "
                    "differ." << endl;
                exit( 1);
            }
        }

        vout << "    ... next is element file ..." << endl;
        in.open( filename[1]);
        if ( ! in) {
            cerr << argv[0] << ": error: cannot open file '"<< filename[1]
             << "' for reading." <<endl;
            exit( 1);
        }
        char  part_name[16000];
        in >> part_name;
        while( in) {
            in.get(c);
            facets.push_back( Facet());
            Facet& facet = (*(facets.end() - 1));
            while( c == ' ') {
                in >> number;
                if ( number < 1 || (size_t)number > points.size()) {
                    cerr << argv[0] << ": error: parsing line "
                         << facets.size() << " from file '" << filename[1]
                         << "': index " << number << " is out of range."
                         << endl;
                    exit( 1);
                }
                facet.push_back( number - 1);
                in.get(c);
            }
            in >> part_name;
        }
        in.close();
        vout << facets.size() << " facets read." << endl;
    }

    const char*  oname    = "cout";
    ostream*     p_out = &cout;
    ofstream     out;
    if ( n > 2) {
        out.open( filename[2]);
        p_out = &out;
        oname = filename[2];
    }
    if ( (ascii_mesh || binary_mesh) && n > 1) {
        out.open( filename[1]);
        p_out = &out;
        oname = filename[1];
    }
    if ( !*p_out) {
        cerr << argv[0] << ": error: cannot open file '"<< oname
             << "' for writing." <<endl;
        exit( 1);
    }

    vout << "CGAL::File_writer_OFF( " << (binary ? ", binary" : ", ASCII")
         << ") ...." << endl;
    CGAL::File_header_OFF  header( binary, noc, skel, verbose);
    CGAL::File_writer_OFF  writer( header);
    writer.write_header( *p_out, points.size(), 0, facets.size(),
                         normals.size()>0);
    Point_vector::iterator norm = normals.begin();
    for ( Point_vector::iterator p = points.begin(); p != points.end(); ++p) {
        writer.write_vertex( (*p).x(), (*p).y(), (*p).z());
        if ( normals.size()>0) {
            writer.write_normal( (*norm).x(), (*norm).y(), (*norm).z());
            ++norm;
        }
    }
    writer.write_facet_header();
    for ( Facet_vector::iterator fc = facets.begin(); fc!=facets.end(); ++fc) {
        writer.write_facet_begin( (*fc).size());
        for ( Facet::iterator fi = (*fc).begin(); fi != (*fc).end(); ++fi) {
            writer.write_facet_vertex_index( *fi);
        }
        writer.write_facet_end();
    }
    writer.write_footer();
    vout << "    .... done." << endl;
    if ( !*p_out) {
        cerr << argv[0] << " write error: while writing file '"<< oname
             << "'."  << endl;
        exit( 1);
    }
    return 0;
}
// EOF //
