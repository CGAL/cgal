// computes bbox of an OFF object.

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cfloat>


using namespace std;

bool  verbose   = false;
bool  unitcube  = false;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
    int n = 0; // number of filenames
    char *filename[1] = { NULL }; // stop compiler warning (too hard to rewrite the code to avoid it)
    bool help = false;
    for (int i = 1; i < argc; i++) { // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
        else if ( strcmp( "-unit", argv[i]) == 0)
            unitcube = true;
        else if ( (strcmp( "-h", argv[i]) == 0) ||
                  (strcmp( "-help", argv[i]) == 0))
            help = true;
        else if ( n < 1 ) {
            filename[ n++] = argv[i];
        } else {
            n++;
            break;
        }
    }
    if ((n > 1) || help) {
        if ( ! help)
            cerr << "Error: in parameter list" << endl;
        cerr << "Usage: " << argv[0] << " [<options>] [<infile> [<outfile>]]"
             << endl;
        cerr << "Usage: " << argv[0] << " [<options>] [<infile>]" << endl;
        cerr << "       computes the bbox of the coordinates of an OFF object."
             << endl;
        cerr << "       -unit     prints transformation to unit cube." << endl;
        cerr << "                 (can be used with 'off_transform')" << endl;
        cerr << "       -v      verbose." << endl;
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
    if ( scanner.size_of_vertices() <= 0) {
        cerr << argv[0] << ": error: file '"<< name
         << "' has no vertices." << endl;
        std::abort();
    }
    size_t  v = scanner.size_of_vertices();
    CGAL::Bbox_3 bbox;
    double x, y, z;
    scanner.scan_vertex( x, y, z);
    bbox = CGAL::Bbox_3( x,y,z, x,y,z);
    v--;
    while (v--) {
        scanner.scan_vertex( x, y, z);
        bbox = bbox + CGAL::Bbox_3( x,y,z, x,y,z);
        scanner.skip_to_next_vertex( scanner.size_of_vertices() - v - 1);
    }
    verr << ".... done." << scanner.size_of_vertices() << " points read."
         << endl;

    if ( !in) {
        cerr << argv[0] << " read error: while reading file '"<< name << "'."
             << endl;
        exit( 1);
    }
    if ( ! unitcube) {
        cout << bbox.xmin() << "  " << bbox.ymin() << "  " << bbox.zmin()
             << '\n';
        cout << bbox.xmax() << "  " << bbox.ymax() << "  " << bbox.zmax()
             << endl;
    } else {
        double s = DBL_MAX;
        double d = bbox.xmax() - bbox.xmin();
        if ( d > 0 && 2/d < s)
            s = 2/d;
        d = bbox.ymax() - bbox.ymin();
        if ( d > 0 && 2/d < s)
            s = 2/d;
        d = bbox.zmax() - bbox.zmin();
        if ( d > 0 && 2/d < s)
            s = 2/d;
        if ( s == DBL_MAX)
            s = 1;
        cout << "-trans  " << (-(bbox.xmin() + bbox.xmax())/2)
             << "  "       << (-(bbox.ymin() + bbox.ymax())/2)
             << "  "       << (-(bbox.zmin() + bbox.zmax())/2)
             << "  -scale  " << s << endl;
    }
    return 0;
}
// EOF //
