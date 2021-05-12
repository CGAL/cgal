#include <CGAL/basic.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/leda_integer.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

using namespace std;
typedef leda_integer                                          NT;
typedef CGAL::Simple_homogeneous<NT>                          R;
typedef CGAL::Polyhedron_3<R>                                 Polyhedron;
typedef Polyhedron::Point_3                                   Point;
typedef CGAL::Nef_polyhedron_3<R>                             Nef_polyhedron;

void read( const char* name, Polyhedron_3& poly) {
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
OutIt fold( InpIt first1,
            ForIt first2, ForIt beyond2,
            OutIt result,
            AdBiFct fct) {
    InpIt beyond1(first1);
    CGAL_For_all(first1,beyond1) {
        for ( ForIt i = first2; i != beyond2; ++i) {
            *result++ = fct( *first1->vertex(), *i);
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
    Add_points add;
    list<Nef_polyhedron> nef;
    Polyhedron::Facet_iterator fit;
    for(fit = P1.facets_begin(); fit != P1.facets_end(); fit++) {
      vector<Point> points;
      std::cerr << "size = " << nef.size() << std::endl;
      fold( fit->facet_begin(),
            P2.vertices_begin(), P2.vertices_end(),
            back_inserter( points),
            add);
      Polyhedron P3;
      convex_hull_3( points.begin(), points.end(), P3);
      nef.push_back(Nef_polyhedron(P3));
    }
    while(nef.size() > 1) {
      std::cerr << "size = " << nef.size() << std::endl;
      Nef_polyhedron Temp(nef.front());
      nef.pop_front();
      nef.push_back(Temp.join(nef.front()));
      nef.pop_front();
    }
    std::cout << nef.front();
    return 0;
}
