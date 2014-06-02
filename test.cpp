#include "MyMink.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>

#include <iostream>
#include <vector>
#include <list>
using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Point_2<Kernel> Point_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

// Check if two polygons with holes are the same.
bool are_equal(const Polygon_with_holes_2& ph1, const Polygon_with_holes_2& ph2)
{
    list<Polygon_with_holes_2> sym_diff;

    CGAL::symmetric_difference(ph1, ph2, std::back_inserter(sym_diff));

    return (sym_diff.empty());
}


int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: <program> (0=old|1=new) polygon-file polygon-file [result-file]" << endl;
        return 0;
    }

    Polygon_2 a, b;
    Polygon_with_holes_2 c;

    Polygon_with_holes_2 (*algo)(const Polygon_2&, const Polygon_2&);
    switch (atoi(argv[1])) {
        case 0:
            algo = &CGAL::minkowski_sum_2;
            break;
        case 1:
            algo = &CGAL::minkowski_sum_2_;
            break;
        default:
            cerr << "Unknown algorithm" << endl;
            return 1;
    }

    fstream data1(argv[2]);
    data1 >> a;
    fstream data2(argv[3]);
    data2 >> b;

    bool verify = argc > 4;
    if (verify) {
        fstream data3(argv[4]);
        data3 >> c;
    }

    Polygon_with_holes_2 sum2 = algo(a, b);

    if (verify) {
        if (!are_equal(sum2, c)) {
            cerr << "Expected:" << endl;
            cerr << c;
            cerr << "Got:" << endl;
            cerr << sum2;
            return 1;
        }
    } else {
        return 2;
    }

    return 0;
}
