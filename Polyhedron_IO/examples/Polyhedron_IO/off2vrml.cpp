// Convert from OFF format to VRML (.wrl) 1.0 or 2.0 format.

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/VRML_1_ostream.h>
#include <CGAL/IO/VRML_2_ostream.h>
#include <CGAL/IO/File_writer_inventor.h>
#include <CGAL/IO/File_writer_VRML_2.h>
#include <CGAL/IO/generic_copy_OFF.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

bool  verbose      = false;
int   version      = 1;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
    int n = 0; // number of filenames
    char *filename[2];
    bool help = false;
    for (int i = 1; i < argc; i++) { // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
        else if ( strcmp( "-2", argv[i]) == 0)
            version = 2;
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
        cerr << "       convert a CGAL object (OFF) to VRML .wrl format."
             << endl;
        cerr << "       -2      VRML 2.0 (default is VRML 1.0)." << endl;
        cerr << "       -v      verbose." << endl;
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
        cerr << argv[0] << ": error: cannot open file '"<< iname
             << "' for reading." <<endl;
        exit( 1);
    }

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

    vout << "CGAL::generic_copy_OFF( " << iname << ", " << oname
         << ", V" << version << " ) ...." << endl;
    if ( version == 1) {
        CGAL::VRML_1_ostream os( *p_out);
        CGAL::File_writer_inventor  writer;
        CGAL::generic_copy_OFF( *p_in, *p_out, writer);
        os.close();
    } else {
        CGAL::VRML_2_ostream os( *p_out);
        CGAL::File_writer_VRML_2  writer;
        CGAL::generic_copy_OFF( *p_in, *p_out, writer);
        os.close();
    }
    vout << "    .... done." << endl;

    if ( !*p_in) {
        cerr << argv[0] << " read error: while reading file '"<< iname << "'."
             << endl;
        exit( 1);
    }
    if ( !*p_out) {
        cerr << argv[0] <<" write error: while writing file '"<< oname << "'."
             << endl;
        exit( 1);
    }
    return 0;
}
// EOF //
