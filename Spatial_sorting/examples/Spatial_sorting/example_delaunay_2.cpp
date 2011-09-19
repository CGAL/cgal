#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/spatial_sort.h>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;

int main ()
{
    std::vector<K::Point_2> v;

    std::cout << "Reading..." << std::endl;
    std::istream_iterator<K::Point_2> begin(std::cin), end;
    std::copy(begin, end, std::back_inserter(v));

    // Comment the following three lines to get a massive speed down.
    std::cout << "Sorting..." << std::endl;
    std::random_shuffle(v.begin(), v.end());
    CGAL::spatial_sort(v.begin(), v.end());

    std::cout << "Delaunay..." << std::endl;
    DT dt;

    DT::Face_handle f;
    for (std::vector<K::Point_2>::const_iterator p = v.begin(); p != v.end(); ++p)
        f = dt.insert(*p, f)->face();

    std::cout << "Delaunay is_valid()..." << std::endl;
    dt.is_valid();

    std::cout << "Ok." << std::endl;

    return 0;
}
