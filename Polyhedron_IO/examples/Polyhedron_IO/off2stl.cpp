// Convert from OFF format to StereoLithography StL format.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

bool  verbose  = false;
bool  binary    = false;

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef CGAL::Polyhedron_3<Kernel>      Polyhedron;
typedef Polyhedron::Vertex_iterator     Vertex_iterator;
typedef Polyhedron::Facet_iterator      Facet_iterator;
typedef Polyhedron::Halfedge_handle     Halfedge_handle;


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
        cerr << "       convert a CGAL object (OFF) to StereoLithography StL "
                "format." << endl;
        cerr << "       -v      verbose." << endl;
        cerr << "       -b      binary." << endl;
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

    vout << "Reading polyhedron ..." << endl;
    Polyhedron P;
    (*p_in) >> P;
    vout << "    .... done." << endl;

    if ( !*p_in) {
        cerr << argv[0] << " read error: while reading file '"<< iname << "'."
             << endl;
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

    vout << "Write file ...." << endl;
    *p_out << "solid " << oname << endl;

    // find "bottom/left/front" corner to translate into positive octant
    Vertex_iterator vi = P.vertices_begin();
    Point p = vi->point();
    double minx = p.x();
    double miny = p.y();
    double minz = p.z();
    for ( ; vi != P.vertices_end(); ++vi) {
        p = vi->point();
        if ( p.x() < minx)
            minx = p.x();
        if ( p.y() < miny)
            miny = p.y();
        if ( p.z() < minz)
            minz = p.z();
    }
    // translate into positive octant
    Vector trans( -minx, -miny, -minz);
    for ( Vertex_iterator i = P.vertices_begin(); i != P.vertices_end(); ++i) {
        i->point() = i->point() + trans;
    }
    // write triangles
    for ( Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_handle h = i->halfedge();
        if ( h->next()->next()->next() != h) {
            cerr << argv[0] << " format error: polyhedron in file '"<<
                iname << "' is not triangulated." << endl;
            exit( 1);
        }
        Point p = h->vertex()->point();
        Point q = h->next()->vertex()->point();
        Point r = h->next()->next()->vertex()->point();
        // compute normal
        Vector n = CGAL::cross_product( q-p, r-p);
        Vector norm = n / std::sqrt( n * n);
        *p_out << "    facet normal " << norm << endl;
        *p_out << "      outer loop " << endl;
        *p_out << "        vertex " << p << endl;
        *p_out << "        vertex " << q << endl;
        *p_out << "        vertex " << r << endl;
        *p_out << "      endloop " << endl;
        *p_out << "    endfacet " << endl;
    }

    *p_out << "endsolid " << oname << endl;
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
