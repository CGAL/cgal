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
// file          : examples/Polyhedron_IO/terr_trian.C
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Delaunay Triangulation of a set of 3D points in the xy-plane. 
// (Terrain triangulation)
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_rational.h>
#endif
#include "triangulation_print_OFF.h"

using namespace std;
using CGAL::Point_3;

template <class K>
struct  Indexed_point: public Point_3<K> {
    int*  index;
    Indexed_point()                                 {}
    Indexed_point( Point_3<K> p) : Point_3<K>(p)    {}
    Indexed_point( double x, double y, double z, int* i) 
        : Point_3<K>(x,y,z), index(i)               {}
};

#ifdef CGAL_USE_LEDA
typedef  CGAL::Cartesian<leda_rational>                     Kernel;
#else
typedef  CGAL::Cartesian<double>                            Kernel;
#endif
typedef  Indexed_point<Kernel>                              IPoint;
typedef  CGAL::Triangulation_euclidean_traits_xy_3<Kernel>  Gtraits;

struct Gt : public Gtraits {
    typedef IPoint Point;
};

typedef  CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef  CGAL::Triangulation_face_base_2<Gt>  Fb;

typedef  CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef  CGAL::Triangulation_2<Gt,Tds>                 Triangulation;
typedef  CGAL::Delaunay_triangulation_2<Gt,Tds>        Delaunay_triangulation;

bool  verbose      = false;
bool  binary       = false;
bool  noc          = false;
bool  delaunay     = false;
bool  incr         = false;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
    int n = 0; // number of filenames
    char *filename[2];
    bool help = false;
    for (int i = 1; i < argc; i++) { // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
	else if ( strcmp( "-b", argv[i]) == 0)
	    binary = true;
	else if ( strcmp( "-noc", argv[i]) == 0)
	    noc = true;
	else if ( strcmp( "-delaunay", argv[i]) == 0)
	    delaunay = true;
	else if ( strcmp( "-incr", argv[i]) == 0)
	    incr = true;
        else if ( (strcmp( "-h", argv[i]) == 0) || 
                  (strcmp( "-help", argv[i]) == 0))
            help = true;
        else if ( n < 2 ) {
            filename[ n++] = argv[i];
        } else {
	    ++n;
            break;
	}
    }
    if ((n > 2) || help) {
        if ( ! help)
            cerr << "Error: in parameter list" << endl;
        cerr << "Usage: " << argv[0] << " [<options>] [<infile> [<outfile>]]"
             << endl;
	cerr << "       Terrain triangulation in the xy-plane." << endl;
	cerr << "       -delaunay  Delaunay triangulation (default)." << endl;
	cerr << "       -incr      Incremental insertion (no flips)." << endl;
	cerr << "       -b         binary output (default is ASCII)." << endl;
	cerr << "       -noc       no comments in file." << endl;
	cerr << "       -v         verbose." << endl;
	exit( ! help);
    }

    CGAL::Verbose_ostream vout( verbose);
    vout << argv[0] << ": verbosity on." << endl;

    const char*  iname = "cin";
    istream*     p_in  = &cin;
    ifstream     in;
    if ( n > 0) {
	in.open( filename[0]);
	p_in = &in;
	iname = filename[0];
    }
    if ( !*p_in) { 
	cerr << argv[0] << ": error: cannot open file '" << iname
	 << "' for reading." <<endl;
	exit( 1);
    }

    CGAL::File_scanner_OFF scanner( * p_in, true);
    if ( !*p_in)
	exit( 1);

    const char*  oname = "cout";
    ostream*     p_out = &cout;
    ofstream     out;
    if ( n > 1) {
	out.open( filename[1]);
	p_out = &out;
	oname = filename[1];
    }
    if ( !*p_out) { 
	cerr << argv[0] << ": error: cannot open file '"<< oname
	     << "' for writing." <<endl;
	exit( 1);
    }

    // index array.
    int* indices = new int[ scanner.size_of_vertices()];
    for ( int k = 0; k < scanner.size_of_vertices(); k++)
	indices[k] = -1;

    if ( delaunay || ! incr) {
	Delaunay_triangulation triang;
	vout << "Scanning and triangulating ..." << endl;
	for ( int j = 0; j < scanner.size_of_vertices(); j++) {
	    double x, y, z;
	    scanner.scan_vertex( x, y, z);
	    IPoint p( x, y, z, indices + j);
	    triang.insert( p);
	}
	vout << "    .... done." << endl;
	vout << "write_triangulation( " << oname 
	     << (binary ? ", binary" : ", ASCII") << ") ...." << endl;
	CGAL::triangulation_print_OFF( *p_out, triang, binary, noc, verbose);
	vout << "    .... done." << endl;
    } else {
        Triangulation triang;
	vout << "Scanning and triangulating ..." << endl;
	for ( int j = 0; j < scanner.size_of_vertices(); j++) {
	    double x, y, z;
	    scanner.scan_vertex( x, y, z);
	    IPoint p( x, y, z, indices + j);
	    triang.insert( p);
	}
	vout << "    .... done." << endl;
	vout << "write_triangulation( " << oname 
	     << (binary ? ", binary" : ", ASCII") << ") ...." << endl;
	CGAL::triangulation_print_OFF( *p_out, triang, binary, noc, verbose);
	vout << "    .... done." << endl;
    }
    if ( !*p_in) { 
	cerr << argv[0] << " read error: while reading file '"<< iname << "'." 
	     << endl;
	exit( 1);
    }
    if ( !*p_out) { 
	cerr << argv[0] << " write error: while writing file '"<< oname << "'." 
	     << endl;
	exit( 1);
    }
    delete[] indices;
    return 0;
}
// EOF //
