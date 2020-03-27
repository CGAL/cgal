// Converts Geometry Information in Inventor files (*.iv) into OFF format.
// The scanner triggers on Coordinate3, IndexedFaceSet, and IndexedLineSet
// keywords. It does not recognize transformations nor groups.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <list>

using namespace std;

typedef  CGAL::Simple_cartesian<double>  Kernel;
typedef  Kernel::Point_3                 Point;
typedef  list<Point>                     Point_list;
typedef  list<std::size_t>                       Facet;
typedef  list<Facet>                     Facet_list;

// Command line arguments
// ======================
bool  verbose      = false;
bool  binary       = false;
bool  skel         = false;
bool  noc          = false;
bool  block        = false;
int   n_block      = 0;


// Global data describing the object
// =================================
Point_list  points;
Facet_list  facets;

// iv File Scanner
// ===============
void iv_file_scanner( istream& in) {
  std::size_t offset    = 0;    // offset for the Index....Set
    char      c;                // one read character (comment, vertex, or facet)
    char      str[2000];        // temporal storage for rest of identifier
    int       blocks    = 0;    // Number of blocks found.

    CGAL::Verbose_ostream vout( verbose);
    in >> str;
    while( in) {
        if ( strcmp( str, "label") == 0) {
            in >> str;
            vout << "Label " << str << " found." << endl;
            in >> str;
        }
        if ( strcmp( str, "Coordinate3") == 0) {
            in >> c;
            if ( c == '{') {
                in >> str;
                if ( strcmp( str, "point") == 0) {
                    vout << "\"Coordinate3 { point\" block found." << endl;
                    blocks ++;
                    if ( ! block || n_block == blocks) {
                        if ( block)
                            vout << "    this block is read." << endl;
                        // point coordinate list starts here
                        offset = points.size();
                        in >> c;
                        CGAL_assertion( c == '[');
                        in >> c;
                        while ( in && ( c != ']')) {
                            in.putback( c);
                            Point p;
                            in >> p;
                            points.push_back( p);
                            in >> c;  // comma or ']' or digit
                            if ( c == ',')
                                in >> c;   // ignoring commas
                        }
                        vout << "    " << points.size() - offset
                             << " point coordinates read." << endl;
                    }
                }
            }
        }
        if ( block && n_block != blocks) {
            in >> str;
            continue;
        }

        if ( strcmp( str, "IndexedFaceSet") == 0) {
            in >> c;
            if ( c == '{') {
                in >> str;
                if ( strcmp( str, "coordIndex") == 0) {
                    vout << "\"IndexedFaceSet { coordIndex\" block found."
                         << endl;
                    // indices start here
                    std::size_t face_offset = facets.size();
                    in >> c;
                    CGAL_assertion( c == '[');
                    facets.push_back( Facet());
                    Facet* facet = &facets.back();
                    in >> c;
                    while ( in && ( c != ']')) {
                        in.putback( c);
                        int index;
                        in >> index;
                        if ( index >= 0) {
                            facet->push_back( index + offset);
                        } else {
                            facets.push_back( Facet());
                            facet = &facets.back();
                        }
                        in >> c;  // comma or ']' or digit
                        if ( c == ',')
                            in >> c;   // ignoring commas
                    }
                    facets.pop_back();
                    vout << "    " << facets.size() - face_offset
                         << " facets read." << endl;
                }
            }
        }
        if ( strcmp( str, "IndexedLineSet") == 0) {
            in >> c;
            if ( c == '{') {
                in >> str;
                if ( strcmp( str, "coordIndex") == 0) {
                    vout << "\"IndexedLineSet { coordIndex\" block found."
                         << endl;
                    // indices start here
                    std::size_t face_offset = facets.size();
                    in >> c;
                    CGAL_assertion( c == '[');
                    facets.push_back( Facet());
                    Facet* facet = &facets.back();
                    in >> c;
                    while ( in && ( c != ']')) {
                        in.putback( c);
                        int index;
                        in >> index;
                        if ( index >= 0) {
                            facet->push_back( index + offset);
                        } else {
                            facets.push_back( Facet());
                            facet = &facets.back();
                        }
                        in >> c;  // comma or ']' or digit
                        if ( c == ',')
                            in >> c;   // ignoring commas
                    }
                    facets.pop_back();
                    vout << "    " << facets.size() - face_offset
                         << " polylines read." << endl;
                }
            }
        }
        in >> str;
    }
    vout << endl;
    vout << "Total number of vertices: " << points.size() << endl;
    vout << "Total number of faces   : " << facets.size() << endl;
    vout << "Total number of blocks  : " << blocks << endl << endl;
}

// OFF File Writer
// ===============
void print_OFF( ostream& out) {
    CGAL::File_header_OFF  header( binary, noc, skel, verbose);
    CGAL::File_writer_OFF  writer( header);
    writer.write_header( out, points.size(), 0, facets.size());
    while( ! points.empty()) {
        writer.write_vertex( points.front().x(),
                             points.front().y(),
                             points.front().z());
        points.pop_front();
    }
    writer.write_facet_header();
    while( ! facets.empty()) {
        Facet& facet = facets.front();
        writer.write_facet_begin( facet.size());
        while( ! facet.empty()) {
            writer.write_facet_vertex_index( facet.front());
            facet.pop_front();
        }
        writer.write_facet_end();
        facets.pop_front();
    }
    writer.write_footer();
}

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
        else if ( strcmp( "-#", argv[i]) == 0) {
            block = true;
            i++;
            if ( i < argc) {
                n_block = atoi( argv[i]);
                if ( n_block < 0) {
                    cerr << argv[0] << ": error: <i> must be greater than "
                            " zero." << endl;
                    help = true;
                }
            } else {
                cerr << argv[0] << ": error: -# needs an integer parameter."
                     << endl;
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
        cerr << "       converts an Inventor file to OFF." << endl;
        cerr << "       -# <i>  convert only i-th block (i in 1..n)." << endl;
        cerr << "       -b      binary (default is ASCII)." << endl;
        cerr << "       -skel   Geomview SKEL format." << endl;
        cerr << "       -noc    no comments in file." << endl;
        cerr << "       -v      verbose." << endl;
        exit( ! help);
    }

    CGAL::Verbose_ostream vout( verbose);
    vout << argv[0] << ": verbosity on." << endl;

    const char* iname = "cin";
    istream*    p_in  = &cin;
    ifstream    in;
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

    const char* oname = "cout";
    ostream*    p_out = &cout;
    ofstream    out;
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

    CGAL::set_ascii_mode(* p_in);
    CGAL::set_ascii_mode(* p_out);

    vout << "Scanning Inventor file `" << iname << "' ....\n--------" << endl;
    iv_file_scanner( *p_in);
    vout << "--------\n    .... done." << endl;

    if ( !*p_in  && !p_in->eof()) {
        cerr << argv[0] << " read error: while reading file '"<< iname << "'."
             << endl;
        exit( 1);
    }

    vout << "print_OFF: `" << iname << "' ...." << endl;
    print_OFF( *p_out);
    vout << "    .... done." << endl;

    if ( !*p_out) {
        cerr << argv[0] << " write error: while writing file '"<< oname
             << "'." << endl;
        exit( 1);
    }
    return 0;
}
// EOF //
