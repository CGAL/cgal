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
        cout << "Usage: <program> polygon-file polygon-file result-file" << endl;
        return 0;
    }

    Polygon_2 a, b;
    Polygon_with_holes_2 c;

    fstream data1(argv[1]);
    data1 >> a;
    fstream data2(argv[2]);
    data2 >> b;
    fstream data3(argv[3]);
    data3 >> c;

    Polygon_with_holes_2 sum2 = CGAL::minkowski_sum_2(a, b);
    if (!are_equal(sum2, c)) {
        cerr << "minkowski_sum_2 NOT OK" << endl;
        cerr << "Expected:" << endl;
        cerr << c;
        cerr << "Got:" << endl;
        cerr << sum2;
        return 1;
    }

    Polygon_with_holes_2 sum = minkowski_sum_2_(a, b);
    if (!are_equal(sum, c)) {
        cerr << "New minkowski_sum_2 NOT OK" << endl;
        cerr << "Expected:" << endl;
        cerr << c;
        cerr << "Got:" << endl;
        cerr << sum;
        return 1;
    }

    return 0;
}
