#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <boost/shared_ptr.hpp>
#include <CGAL/Gmpq.h>

#include <fstream>

#ifndef CGAL_USE_GMP
#define CGAL_USE_GMP
#endif

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_2                                         Point_2;
typedef CGAL::Polygon_2<Kernel>                                 Polygon_2;
typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>        SsBuilderTraits;
typedef CGAL::Straight_skeleton_2<Kernel>                       Ss;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss>   SsBuilder;

void print_polygon_in_exact_form(const Polygon_2& poly);
void print_polygon_in_double(const Polygon_2& poly);
void read_polygon_from_file(Polygon_2& pgn);

int main() {
    Polygon_2 pgn;
    read_polygon_from_file(pgn);
    if(!pgn.is_counterclockwise_oriented()) {
        std::cerr << "Polygon is not CCW Oriented" << std::endl;
    }
    if(!pgn.is_simple()) {
        std::cerr << "Polygon is not simple" << std::endl;
    }

    // print_polygon_in_exact_form(pgn);  // Observe that the polygon is read from the input file correctly
    // print_polygon_in_double(pgn);

    SsBuilder ssb;
    ssb.enter_contour(pgn.vertices_begin(), pgn.vertices_end());
    boost::shared_ptr<Ss> straight_ske = ssb.construct_skeleton();
    if(!straight_ske->is_valid()) {
        std::cerr << "Straight skeleton is not valid" << std::endl;
    }
    return 0;
}

void print_polygon_in_exact_form(const Polygon_2& poly) {
    std::cout << "Polygon with " << poly.size() << " vertices" << std::endl ;
    for(Polygon_2::Vertex_const_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi) {
        CGAL::Gmpq x (CGAL::exact(vi->x()));
        CGAL::Gmpq y (CGAL::exact(vi->y()));
        std::cout << x.numerator() << "/" << x.denominator() << std::endl;
        std::cout << y.numerator() << "/" << y.denominator() << std::endl;
    }
}

void print_polygon_in_double(const Polygon_2& poly) {
    std::cout.precision(20);
    std::cout << "Polygon with " << poly.size() << " vertices" << std::endl ;
    for(Polygon_2::Vertex_const_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi) {
        std::cout << "(" << vi->x() << "," << vi->y() << ")" << std::endl;
    }
}

void read_polygon_from_file(Polygon_2& pgn) {

  std::ifstream in("data/issue3382.txt");

  CGAL::Gmpq gx, gy;
  while(in >> gx  >> gy) {
    FT x(gx);
    FT y(gy);
    Point_2 p(x, y);
    pgn.push_back(p);
  }
}
