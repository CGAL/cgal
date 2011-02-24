// Applies translation and scaling to an OFF object.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

typedef  CGAL::Simple_cartesian<double>  Kernel;
typedef  Kernel::Point_3                 Point;
typedef  Kernel::Vector_3                Vector;

bool    verbose      = false;
bool    binary       = false;
bool    skel         = false;
bool    noc          = false;

double  transx       =  0.0;
double  transy       =  0.0;
double  transz       =  0.0;
double  scale        =  1.0;

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
        else if ( strcmp( "-skel", argv[i]) == 0)
            skel = true;
        else if ( strcmp( "-noc", argv[i]) == 0)
            noc = true;
        else if ( strcmp( "-scale", argv[i]) == 0) {
            i++;
            if ( i < argc) {
                scale = atof( argv[i]);
            } else {
                cerr << argv[0] << ": error: -scale needs a double parameter."
                     << endl;
                help = true;
            }
        } else if ( strcmp( "-trans", argv[i]) == 0) {
            i++;
            if ( i+2 < argc) {
                transx = atof( argv[i]);
                i++;
                transy = atof( argv[i]);
                i++;
                transz = atof( argv[i]);
            } else {
                cerr << argv[0] << ": error: -trans needs three double "
                        "parameters." << endl;
                help = true;
            }
        } else if ( (strcmp( "-h", argv[i]) == 0) ||
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
        cerr << "       transforms coordinate values of an OFF object."
             << endl;
        cerr << "       -trans <x> <y> <y>    translation." << endl;
        cerr << "       -scale <s>            uniform scaling." << endl;
        cerr << "       -b                    binary (default is ASCII)."
             << endl;
        cerr << "       -skel                 Geomview SKEL format." << endl;
        cerr << "       -noc                  no comments in file." << endl;
        cerr << "       -v                    verbose." << endl;
        exit( ! help);
    }

    CGAL::Verbose_ostream verr( verbose);
    verr << argv[0] << ": verbosity on." << endl;

    const char*  name = "cin";
    istream*     p_in = &cin;
    ifstream     in;
    if ( n > 0) {
        in.open( filename[0]);
        p_in = &in;
        name = filename[0];
    }
    if ( ! * p_in) {
        cerr << argv[0] << ": error: cannot open file '"<< name
         << "' for reading." <<endl;
        exit( 1);
    }

    verr << "CGAL::File_scanner_OFF( " << name << ") ...." << endl;
    CGAL::File_scanner_OFF scanner( * p_in);
    if ( ! * p_in) {
        cerr << argv[0] << ": error: file '"<< name
         << "' is not in OFF format." << endl;
        std::abort();
    }
    verr << "CGAL::File_writer_OFF( ..." << endl;
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

    CGAL::File_header_OFF  header( binary, noc, skel, verbose);
    CGAL::File_writer_OFF  writer( header);
    writer.write_header( * p_out,
                         scanner.size_of_vertices(),
                         scanner.size_of_halfedges(),
                         scanner.size_of_facets());

    Vector v( transx, transy, transz);
    for ( std::size_t k = 0; k < scanner.size_of_vertices(); ++k) {
        double x, y, z;
        scanner.scan_vertex( x, y, z);
        Point q( x, y, z);
        q = q + v;
        q = CGAL::ORIGIN + ( (q - CGAL::ORIGIN) * scale );
        scanner.skip_to_next_vertex( k);
        writer.write_vertex( q.x(), q.y(), q.z());
    }
    verr << "    .... done." << scanner.size_of_vertices() << " points read."
	 << endl;

    if ( ! *p_in) {
        cerr << argv[0] << " read error: while reading file '"<< name << "'."
             << endl;
        exit( 1);
    }
    writer.write_facet_header();
    * p_out << endl;
    char c;
    while ( (*p_in).get(c))
        * p_out << c;
    return 0;
}
// EOF //
