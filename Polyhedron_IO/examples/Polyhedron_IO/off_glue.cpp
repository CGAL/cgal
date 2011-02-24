// Glue vertices of a polyhedron together that have equal coordinate values.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

typedef  CGAL::Simple_cartesian<double>  Kernel;
typedef  Kernel::Point_3                 Point;
typedef  Kernel::Vector_3                Vector;
struct   Vertex {
    Point   point;
    Vector  normal;
    int     index;
};
typedef  vector<Vertex>           Vertex_vector;
typedef  Vertex_vector::iterator  VIterator;

struct   VertexComp {
    bool operator()( const Vertex* v, const Vertex* w) const {
        return ( v->point.x() < w->point.x() ||
            (v->point.x() == w->point.x() && v->point.y() < w->point.y()) ||
            (v->point.x() == w->point.x() && (v->point.y() == w->point.y()
                                           && v->point.z() < w->point.z())));
    }
};

bool  verbose = false;
bool  binary  = false;
bool  skel    = false;
bool  noc     = false;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
    int n = 0; // number of filenames
    char *filename[2];
    bool help = false;
    for (int i=1 ; i < argc; i++) { // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
        else if ( strcmp( "-b", argv[i]) == 0)
            binary = true;
        else if ( strcmp( "-skel", argv[i]) == 0)
            skel = true;
        else if ( strcmp( "-noc", argv[i]) == 0)
            noc = true;
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
        cerr << "       glues vertices of equal coordinates together." << endl;
        cerr << "       -b      binary output (default is ASCII)." << endl;
        cerr << "       -skel   Geomview SKEL format." << endl;
        cerr << "       -noc    no comments in file." << endl;
        cerr << "       -v      verbose." << endl;
        exit( ! help);
    }

    CGAL::Verbose_ostream vout( verbose);
    vout << argv[0] << ": verbosity on." << endl;

    const char*  name = "<cin>";
    istream*     p_in = &cin;
    ifstream     in;
    if ( n > 0) {
        in.open( filename[0]);
        p_in = &in;
        name = filename[0];
    }
    if ( !in) {
        cerr << argv[0] << ": error: cannot open file '"<< name
         << "' for reading." <<endl;
        exit( 1);
    }

    CGAL::File_scanner_OFF  scanner( * p_in);
    if ( ! (*p_in)) {
        cerr << " " << endl;
        cerr << argv[0] << ": input error: file format is not in OFF." << endl;
        std::abort();
    }
    Vertex_vector            vertices;
    vector<Vertex*>          sorted_vertices;
    // Avoid any reallocation
    vertices.reserve(        scanner.size_of_vertices());
    sorted_vertices.reserve( scanner.size_of_vertices());

    float x, y, z;
    for (std::size_t i = 0; i < scanner.size_of_vertices(); i++) {
        Vertex  vertex;
        scanner.scan_vertex( x, y, z);
        vertex.point = Point( x, y, z);
        //scanner.scan_normal( x, y, z);
        vertex.normal = Vector( x, y, z);
        scanner.skip_to_next_vertex( i);
        vertex.index = -1;
        vertices.push_back( vertex);
        sorted_vertices.push_back( & vertices.back());
    }
    vout << scanner.size_of_vertices() << " vertices read." << endl;

    sort( sorted_vertices.begin(), sorted_vertices.end(), VertexComp());
    int current_index = 0;
    sorted_vertices[0]->index = 0;
    for (std::size_t i = 1; i < scanner.size_of_vertices(); i++) {
        if ( sorted_vertices[i]->point != sorted_vertices[i-1]->point)
            current_index++;
        sorted_vertices[i]->index = current_index;
    }
    current_index++;
    vout << "Merged to " << current_index << " vertices." << endl;

    const char*  oname = "<cout>";
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

    vout << "CGAL::File_writer_OFF( " << (binary ? ", binary" : ", ASCII")
         << ") ...." << endl;
    CGAL::File_header_OFF  header( binary, noc, skel, verbose);
    CGAL::File_writer_OFF  writer( header);
    writer.write_header(*p_out,
                        current_index,
                        0,
                        scanner.size_of_facets(),
                        scanner.has_normals());
    vector<Vertex*>::iterator v = sorted_vertices.begin();
    writer.write_vertex((*v)->point.x(), (*v)->point.y(), (*v)->point.z());
    if ( scanner.has_normals()) {
        writer.write_normal((*v)->normal.x(),
                            (*v)->normal.y(),
                            (*v)->normal.z());
    }
    ++v;
    for ( ; v != sorted_vertices.end(); ++v) {
        if ( (*v)->index != (*(v-1))->index) {
            writer.write_vertex( (*v)->point.x(),
                                 (*v)->point.y(),
                                 (*v)->point.z());
            if ( scanner.has_normals()) {
                writer.write_normal( (*v)->normal.x(),
                                     (*v)->normal.y(),
                                     (*v)->normal.z());
            }
        }
    }

    // Copy facets and translate vertex indices.
    writer.write_facet_header();
    for (std::size_t i = 0; i < scanner.size_of_facets(); i++) {
      std::size_t no;  // number of vertices of a facet.
        scanner.scan_facet( no, i);
        writer.write_facet_begin( no);
        for ( std::size_t j = 0; j < no; j++) {
          std::size_t index;
            scanner.scan_facet_vertex_index( index, i);
            writer.write_facet_vertex_index(  vertices[index].index);
        }
        scanner.skip_to_next_facet( i);
        writer.write_facet_end();
    }
    writer.write_footer();
    vout << "    .... done." << endl;

    if ( ! * p_in) {
        cerr << argv[0] << " read error: while reading file '"<< name << "'."
             << endl;
        exit( 1);
    }
    if ( !*p_out) {
        cerr << argv[0] << " write error: while writing file '"<< oname
             << "'."  << endl;
        exit( 1);
    }
    return 0;
}
