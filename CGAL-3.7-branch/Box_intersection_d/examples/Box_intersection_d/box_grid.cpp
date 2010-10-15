#include <CGAL/box_intersection_d.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>

typedef CGAL::Box_intersection_d::Box_d<int,2> Box;

// coordinates for 9 boxes of a grid
int p[9*4]   = { 0,0,1,1,  1,0,2,1,  2,0,3,1, // lower
                 0,1,1,2,  1,1,2,2,  2,1,3,2, // middle
                 0,2,1,3,  1,2,2,3,  2,2,3,3};// upper
// 9 boxes
Box boxes[9] = { Box( p,    p+ 2),  Box( p+ 4, p+ 6),  Box( p+ 8, p+10),
                 Box( p+12, p+14),  Box( p+16, p+18),  Box( p+20, p+22),
                 Box( p+24, p+26),  Box( p+28, p+30),  Box( p+32, p+34)};
// 2 selected boxes as query; center and upper right
Box query[2] = { Box( p+16, p+18),  Box( p+32, p+34)};

// callback function object writing results to an output iterator
template <class OutputIterator>
struct Report {
    OutputIterator it;
    Report( OutputIterator i) : it(i) {} // store iterator in object
    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const Box& a, const Box&) { *it++ = a.id(); }
};
template <class Iter> // helper function to create the function object
Report<Iter> report( Iter it) { return Report<Iter>(it); }

int main() {
    // run the intersection algorithm and store results in a vector
    std::vector<std::size_t> result;
    CGAL::box_intersection_d( boxes, boxes+9, query, query+2,
                              report( std::back_inserter( result)));
    // sort, check, and show result
    std::sort( result.begin(), result.end());
    std::size_t check1[13] = {0,1,2,3,4,4,5,5,6,7,7,8,8};
    assert(result.size() == 13 && std::equal(check1,check1+13,result.begin()));
    std::copy( result.begin(), result.end(),
               std::ostream_iterator<std::size_t>( std::cout, " "));
    std::cout << std::endl;

    // run it again but for different cutoff value and half-open boxes
    result.clear();
    CGAL::box_intersection_d( boxes, boxes+9, query, query+2,
                              report( std::back_inserter( result)),
                              std::ptrdiff_t(1),
                              CGAL::Box_intersection_d::HALF_OPEN);
    // sort, check, and show result
    std::sort( result.begin(), result.end());
    std::size_t check2[2]  = {4,8};
    assert(result.size() == 2 && std::equal(check2, check2+2, result.begin()));
    std::copy( result.begin(), result.end(),
               std::ostream_iterator<std::size_t>( std::cout, " "));
    std::cout << std::endl;
    return 0;
}
