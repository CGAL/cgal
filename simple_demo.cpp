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
    if (argc < 3) {
        cout << "Usage: <program> file1 file2" << endl;
        return 0;
    }

    Polygon_2 polygon_a, polygon_b;

    fstream data1(argv[1]);
    data1 >> polygon_a;
    fstream data2(argv[2]);
    data2 >> polygon_b;

    Polygon_with_holes_2 sum = minkowski_sum_2_(polygon_a,polygon_b);
    cout << sum << endl << endl;

    Polygon_with_holes_2 sum2 = CGAL::minkowski_sum_2(polygon_a,polygon_b);
    cout << sum2 << endl;

    if (are_equal(sum, sum2)) {
        cout << "Polygons are equal\n";
    } else {
        cout << "Polygons are NOT equal\n";
    }
}
