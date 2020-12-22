#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

// #include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/convex_hull_3.h>

using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel     R;
//typedef CGAL::Cartesian<double>                               R;
typedef CGAL::Polyhedron_3<R>                                 Polyhedron;
typedef Polyhedron::Point_3                                   Point;

void read( const char* name, Polyhedron& poly) {
    ifstream in( name);
    if ( ! in) {
        cerr << "minkowsky_sum: error: cannot open file '"<< name
             << "' for reading." << endl;
        exit( 1);
    }
    in >> poly;
    if ( ! in) {
        cerr << "minkowsky_sum: error: reading from file '"<< name << "'."
             << endl;
        exit( 1);
    }
}

template <class InpIt, class ForIt, class OutIt, class AdBiFct>
OutIt fold( InpIt first1, InpIt beyond1,
            ForIt first2, ForIt beyond2,
            OutIt result,
            AdBiFct fct) {
    for ( ; first1 != beyond1; ++first1) {
        for ( ForIt i = first2; i != beyond2; ++i) {
            *result++ = fct( *first1, *i);
        }
    }
    return result;
}

struct Add_points {
    template <class Vertex>
    Point operator()( const Vertex& p, const Vertex& q) const {
        using CGAL::ORIGIN;
        return ORIGIN + (p.point()-ORIGIN) + (q.point()-ORIGIN);
    }
};

int main( int argc, char **argv) {
    if ( argc != 3) {
        cerr << "Usage: " << argv[0] << " <infile1> <infile2>" << endl;
        cerr << "       Minkowsky sum of two 3d polyhedra in OFF format."
             << endl;
        cerr << "       Output in OFF to stdout." << endl;
        exit(1);
    }
    Polyhedron P1;
    Polyhedron P2;
    read( argv[1], P1);
    read( argv[2], P2);

    CGAL::Timer t;
    t.start();

    vector<Point> points;
    Add_points add;
    fold( P1.vertices_begin(), P1.vertices_end(),
          P2.vertices_begin(), P2.vertices_end(),
          back_inserter( points),
          add);
    Polyhedron P3;
    convex_hull_3( points.begin(), points.end(), P3);

    t.stop();
    std::cout << "Runtime Minkowski Sum: " << t.time() << std::endl;

    //    cout << P3;
    return 0;
}
